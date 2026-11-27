"""Byte-range download scheduler.

Every file being fetched, across all datasets, becomes part of one pool of byte ranges, served by a shared set of curl connections whose number is tuned to the highest measured throughput.

Throughput on a path is often limited per TCP connection (a VPN tunnel, or a server pacing each connection) rather than in total.
With n connections the total rate R(n) rises until the link or the server saturates, and fetching B bytes takes B / R(n).
The scheduler therefore

  * keeps every connection busy until the last byte.
    Ranges come from one queue, ordered by the files' request order so whole files finish early and can be checksummed while the rest downloads.
    Once the queue is empty, an idle connection takes the far half of the largest range still in flight by work stealing, so no connection idles before the last byte.
  * picks n by measurement.
    It raises n while the measured total keeps rising by more than GAIN, falls back to the best n seen, and probes again every REPROBE_S seconds, since the link changes.
  * drops a connection that stalls (under 10 kB/s for STALL_S seconds) and requeues its unfetched bytes.
    Near the end, a range that was receiving bytes and then made no progress for TAKEOVER_S seconds is handed to an idle connection.
    A range still waiting for its first byte is not, since a queueing server would only put the new request at the back of its queue.
  * learns each host's connection limit.
    A host starts at HOST_START connections.
    When the measured total is still rising at the limit, the probe for more connections raises the limit of every host at it, so a server that takes more connections gets them, and the measurement then keeps them only if they pay.
    A first byte slower than SLOW_FIRST_S, the sign of a server queueing connections, lowers the host's limit below the connections it had open, and the host is not raised again for REPROBE_S.
    A 429 or 503 answer cuts the host's limit to three quarters of the connections it had open and pauses the host for its Retry-After, or for a growing randomized backoff without one.
    A rejection is not a failure of the range, which waits and retries, so a busy server slows the fetch without failing it.

Every answer must be a 206 whose Content-Range starts at the requested byte before any of it is written.
A 200 with the whole file, a proxy's error page or a range for another offset is discarded, so the file can never end up the right size with the wrong bytes.
A range fails its file only after RANGE_MAX_ATTEMPTS attempts in a row without a byte, and at once when the server's total no longer matches the expected size.

Bytes are written in place (pwrite) into a preallocated <dest>.part, and <part>.lock keeps a second fetch process from writing the same file.
The completed intervals are saved to <dest>.part.done only after an fdatasync, so a restart resumes from what has reached the disk.
An empty .done is written as soon as the .part is preallocated, so a full-size .part never exists without one.
A plain .part prefix left by a single-connection curl, and the .part.N/.part.ranges segments of the earlier fixed-split fetch, are taken over as completed intervals.

Only curl and the Python standard library are used, so it runs on any Linux."""

import contextlib
import fcntl
import json
import os
import random
import re
import subprocess
import threading
import time
import urllib.parse

from units import BYTES_PER_GB, BYTES_PER_MB, BYTES_PER_MIB

CHUNK = 32 << 20  # smallest initial range
# Initial ranges are a sixteenth of the file, up to CHUNK_MAX, so a fast link makes few requests
CHUNK_MAX = 256 << 20
MIN_SPLIT = (
    4 << 20
)  # a range in flight is split only when both halves get at least this
MAX_CONNS = 64  # hard ceiling (--max-connections)
START_CONNS = 4
WINDOW_S = 10.0  # throughput is measured over this window
GAIN = 0.05  # more connections are kept only for a total at least 5% higher
REPROBE_S = 120.0
LOG_S = 60
STALL_S = 30
TAKEOVER_S = 10
# Servers limit connections per client.
# On 2026-10-16 ENA served 16 concurrent ranges at once, each with its first byte within 2.6 s.
# At 24 it delayed some by 6.5 s, and at 32 it left 4 without a byte for 40 s.
# Past such a limit new connections only wait in line, so each host gets a cap.
HOST_START = 16
SLOW_FIRST_S = 5.0  # a first byte slower than this means the host queues connections
SLOW_ANSWER_MIN_CAP = START_CONNS  # slow answers alone never cut a host below this many connections
REJECT_CODES = ("429", "503")  # the server asks the client to slow down
REJECT_MAX_WAIT_S = 120.0  # longest pause of a host without a Retry-After
REJECT_KEEP_NUM, REJECT_KEEP_DEN = 3, 4  # a rejection keeps 3/4 of the open connections
REJECTS_MAX = 200  # rejections of one range before its file fails
TRACE = os.environ.get("RANGEFETCH_TRACE")  # file that gets one line per range attempt
RANGE_MAX_ATTEMPTS = 20  # failed attempts of one range before its file fails
RETRY_SLEEP_PER_ATTEMPT_S = 2  # a failed range waits this many seconds per attempt so far before it is retried ...
RETRY_SLEEP_MAX_S = 30  # ... up to this long
BACKOFF_BASE = 2.0  # a rejected host pauses BACKOFF_BASE ** (rejections in a row) seconds ...
BACKOFF_JITTER_MIN = 0.5  # ... times a random factor in [BACKOFF_JITTER_MIN, BACKOFF_JITTER_MIN + 1)
INITIAL_RANGES = 16  # a file is first cut into this many ranges, within [CHUNK, CHUNK_MAX]
CAP_RAISE_STEP = 2  # a connection cap grows by at least this many connections ...
CAP_RAISE_FACTOR = 1.5  # ... or by this factor, whichever is more
BEST_RATE_DECAY = 0.7  # a connection count's best rate decays by this factor per window, so stale peaks fade
DONE_RANGES_MERGE_AT = 64  # finished ranges of a file are merged once there are more than this many
SEGMENT_READ_BYTES = 16 * BYTES_PER_MIB  # read size when a segment file is copied into the part file
CURL_READ_BYTES = BYTES_PER_MIB  # read size of curl's output
SYNC_INTERVAL_S = 5  # files are fdatasynced and their progress recorded this often
NO_FIRST_BYTE = -1  # trace value of a range attempt that got no byte
PART_FILE_MODE = 0o644  # permissions of a new .part file
NO_HTTP_CODE = "000"  # curl's %{http_code} when no response arrived
CURL = ["curl", "-fsSL", "--speed-limit", "10000", "--speed-time", str(STALL_S)]


# A range is split in two, and only when both halves keep MIN_SPLIT bytes.
SPLIT_PARTS = 2


def _merge(intervals):
    """Sorted union of the [start, end) intervals `intervals`, with overlapping and touching ones joined."""
    out = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return out


def _gaps(done, size):
    """The [start, end) intervals of [0, size) not covered by `done`."""
    out, pos = [], 0
    for start, end in _merge(done):
        if start > pos:
            out.append([pos, start])
        pos = max(pos, end)
    if pos < size:
        out.append([pos, size])
    return out


def supports_ranges(url):
    """True when the server answers a one-byte range with 206, False for any other answer, and None when curl got no answer.

    A server that ignores Range would send the whole file to every connection.
    None keeps a dropped probe from being remembered as "no ranges", which would send the file down the single-connection path and discard a partial .part."""
    result = subprocess.run(
        [
            "curl",
            "-sSL",
            "-r",
            "0-0",
            "-o",
            os.devnull,
            "-w",
            "%{http_code}",
            "--retry",
            "4",
            "--retry-all-errors",
            url,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    code = result.stdout.strip()[-len(NO_HTTP_CODE):]
    if code in ("", NO_HTTP_CODE):
        return None
    return code == "206"


@contextlib.contextmanager
def part_lock(part):
    """Hold an exclusive lock on <part>.lock, so two fetch processes never write one .part with separate .done files."""
    lock_path = part + ".lock"
    fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, PART_FILE_MODE)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit(f"{part}: another fetch is writing this file") from None
        try:
            yield
        finally:
            os.remove(lock_path)
    finally:
        os.close(fd)


def content_length(url):
    """Content-Length of url after redirects, or None when the server gives none."""
    result = subprocess.run(
        ["curl", "-sSIL", "--retry", "4", "--retry-all-errors", url],
        capture_output=True,
        text=True,
        check=False,
    )
    length = None
    for line in result.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            length = int(
                line.split(":", 1)[1]
            )  # the last one belongs to the final redirect target
    return length


def _response(head):
    """(HTTP status, Retry-After seconds or None) of the last response in the header text `head`."""
    blocks = [block for block in re.split(r"\r?\n\r?\n", head) if block.strip()]
    if not blocks:
        return None, None
    last = blocks[-1]
    status_match = re.match(r"HTTP/\S+\s+(\d{3})", last)
    retry_match = re.search(r"(?im)^retry-after:\s*(\d+)\s*$", last)
    return (status_match.group(1) if status_match else None), (float(retry_match.group(1)) if retry_match else None)


# Under curl -L, redirects, 100 Continue and a proxy's CONNECT answer each print a header block before the real response.
_INTERIM_STATUS = re.compile(rb"HTTP/\S+\s+(1\d\d|3\d\d|200 connection established)", re.IGNORECASE)
_HEAD_END = re.compile(rb"\r?\n\r?\n")


def _read_head(stream):
    """Split curl -D - output into (final response's header text, body bytes read past it).

    The body is None when the stream ended before a final header block, as when curl gave up waiting for an answer."""
    buf = b""
    while True:
        end = _HEAD_END.search(buf)
        if end:
            if _INTERIM_STATUS.match(buf):
                buf = buf[end.end() :]
                continue
            return buf[: end.start()].decode("latin-1"), buf[end.end() :]
        chunk = stream.read(CURL_READ_BYTES)
        if not chunk:
            return buf.decode("latin-1"), None
        buf += chunk


def _range_problem(head, range_lo):
    """Why `head` is not a 206 for bytes starting at range_lo, or None when it is.

    A proxy, a VPN error page or a server ignoring Range answers with something else of the same length, which would fill the range with the wrong bytes."""
    status, _ = _response(head)
    if status != "206":
        return f"status {status}"
    served = re.search(r"(?im)^content-range:\s*bytes\s+(\d+)-", head)
    if not served:
        return "206 without Content-Range"
    if int(served.group(1)) != range_lo:
        return f"bytes from {served.group(1)} instead of {range_lo}"
    return None


def _served_size(head):
    """Total file size in the Content-Range of `head`, or None when it gives none."""
    total = re.search(r"(?im)^content-range:\s*bytes\s+\d+-\d+/(\d+)", head)
    return int(total.group(1)) if total else None


class _Host:
    """The learned connection limit of one host and the time it may be asked again."""

    def __init__(self):
        self.cap = HOST_START
        self.until = 0.0  # no new connection before this time
        self.backoff = 0  # consecutive rejection rounds
        self.cut = 0.0  # time the limit was last lowered


class _File:
    """One file being fetched, with its completed intervals and its error state."""

    def __init__(self, url, part, size, order):
        self.url, self.part, self.size, self.order = url, part, size, order
        self.lock = threading.Lock()
        self.done = []  # completed [start, end) intervals, in memory
        self.synced = []  # the part of done known to be on disk
        self.left = 0  # bytes still to write
        self.err = None
        self.event = threading.Event()
        self.part_fd = None
        self.closed = False


class _Range:
    """A byte range [start, end) of a file assigned to one connection."""

    def __init__(self, file_state, start, end):
        self.file, self.start, self.end = (
            file_state,
            start,
            end,
        )  # start advances as bytes arrive, end may shrink when stolen from
        self.attempts = 0  # attempts in a row that ended without a byte
        self.rejects = 0
        self.rejected = False
        self.receiving = False  # the current attempt has delivered body bytes
        self.seen = time.time()  # time of the last progress
        self.proc = None
        self.host = urllib.parse.urlsplit(file_state.url).netloc


class Scheduler:
    """Pool of curl connections that serves the byte ranges of every file being fetched."""

    def __init__(self, max_conns=MAX_CONNS, log=print, progress=None):
        self.max_conns, self.log, self.progress = max_conns, log, progress
        self.lock = threading.Condition()
        self.queue, self.flight, self.files = [], [], []
        self.target = START_CONNS
        self.active = 0
        self.bytes = 0
        self.order = 0
        self.workers = 0
        self.started = False
        self.steals = self.takeovers = self.rejections = 0
        self.hosts = {}  # host -> _Host
        self.wrong_answers = set()  # (host, problem kind) pairs already logged

    # Files
    def fetch(self, url, part, size):
        """Fill the size bytes of `part` from url, blocking until done and raising on failure."""
        with self.lock:
            file_state = _File(url, part, size, self.order)
            self.order += 1
        done = self._adopt(part, size)
        file_state.part_fd = os.open(part, os.O_RDWR | os.O_CREAT, PART_FILE_MODE)
        os.ftruncate(file_state.part_fd, size)
        # The .done file is written before any range starts, so a full-size .part is never left without one.
        # Otherwise a crash before the first sync would leave a preallocated .part that _adopt takes for a complete prefix.
        self._save(part, done)
        file_state.done, file_state.synced = [list(interval) for interval in done], [list(interval) for interval in done]
        gaps = _gaps(done, size)
        file_state.left = sum(end - start for start, end in gaps)
        if file_state.left < size:
            self.log(
                f"    resuming {os.path.basename(part)}: {size - file_state.left} of {size} bytes already on disk",
                flush=True,
            )
        if file_state.left == 0:
            os.close(file_state.part_fd)
            self._clear_state(part)
            return
        chunk = max(CHUNK, min(CHUNK_MAX, size // INITIAL_RANGES))
        with self.lock:
            for start, end in gaps:
                for c in range(start, end, chunk):
                    self.queue.append(_Range(file_state, c, min(c + chunk, end)))
            self.queue.sort(key=lambda byte_range: (byte_range.file.order, byte_range.start))
            self.files.append(file_state)
            self._start()
            self.lock.notify_all()
        file_state.event.wait()
        with self.lock:
            # A file failed, so its queued ranges are dropped and its in-flight ones stop.
            self.queue = [byte_range for byte_range in self.queue if byte_range.file is not file_state]
            while any(byte_range.file is file_state for byte_range in self.flight):
                self.lock.wait(1.0)
            self.files.remove(file_state)
        with file_state.lock:
            file_state.closed = True
            os.fsync(file_state.part_fd)
            os.close(file_state.part_fd)
        if file_state.err:
            raise SystemExit(f"{part}: {file_state.err}")
        self._clear_state(part)

    def _adopt(self, part, size):
        """Completed intervals of `part` that are already on disk."""
        state = part + ".done"
        if os.path.exists(state):
            with open(state) as handle:
                return json.load(handle)
        done = []
        if os.path.exists(part):  # a prefix from a single-connection download
            done.append([0, min(os.path.getsize(part), size)])
        plan = part + ".ranges"  # fixed-split segments left by the earlier fetch
        if os.path.exists(plan):
            with open(plan) as handle:
                bounds = json.load(handle)
            # Old segment files are copied into place in the full-size .part file. Their spans then count as done.
            out_fd = os.open(part, os.O_RDWR | os.O_CREAT, PART_FILE_MODE)
            os.ftruncate(out_fd, size)
            for n, (start, end) in enumerate(bounds):
                seg = f"{part}.{n}"
                if not os.path.exists(seg):
                    continue
                with open(seg, "rb") as seg_handle:
                    pos = start
                    for chunk in iter(lambda: seg_handle.read(SEGMENT_READ_BYTES), b""):
                        # Segment bounds were inclusive, so a segment may not write past end.
                        chunk = chunk[: end + 1 - pos]
                        os.pwrite(out_fd, chunk, pos)
                        pos += len(chunk)
                done.append([start, pos])
            os.fsync(out_fd)
            os.close(out_fd)
            self._save(part, done)
            for n in range(len(bounds)):
                if os.path.exists(f"{part}.{n}"):
                    os.remove(f"{part}.{n}")
            os.remove(plan)
        return _merge(done)

    def _save(self, part, done):
        """Atomically replace part's .done file with the merged intervals."""
        tmp = part + ".done.tmp"
        with open(tmp, "w") as handle:
            json.dump(_merge(done), handle)
        os.replace(tmp, part + ".done")

    def _clear_state(self, part):
        """Remove part's interval files once the download is complete."""
        for path in (part + ".done", part + ".done.tmp"):
            if os.path.exists(path):
                os.remove(path)

    # Workers
    def _start(self):
        """Start the controller and syncer threads once, and workers up to max_conns."""
        if not self.started:
            self.started = True
            threading.Thread(target=self._control, daemon=True).start()
            threading.Thread(target=self._syncer, daemon=True).start()
        while self.workers < self.max_conns:
            threading.Thread(
                target=self._worker, args=(self.workers,), daemon=True
            ).start()
            self.workers += 1

    def _next(self, slot):
        """The next range for worker `slot`, or None when it should wait."""
        if slot >= self.target:
            return None
        per = {}
        for byte_range in self.flight:
            per[byte_range.host] = per.get(byte_range.host, 0) + 1
        now = time.time()

        def room(host):
            """True when host `host` is not paused and is below its learned limit."""
            host_state = self._host(host)
            return host_state.until <= now and per.get(host, 0) < host_state.cap

        for n, byte_range in enumerate(self.queue):
            if room(byte_range.host):
                return self.queue.pop(n)
        open_ = {byte_range.host for byte_range in self.flight if room(byte_range.host)}
        # Work stealing takes the far half of the largest range in flight.
        big = max(
            (byte_range for byte_range in self.flight if byte_range.host in open_),
            key=lambda byte_range: byte_range.end - byte_range.start,
            default=None,
        )
        if big is not None and big.end - big.start >= SPLIT_PARTS * MIN_SPLIT:
            mid = big.start + (big.end - big.start) // SPLIT_PARTS
            byte_range = _Range(big.file, mid, big.end)
            big.end = mid
            self.steals += 1
            return byte_range
        # In the endgame, a range without progress for TAKEOVER_S moves whole to this idle connection, and the stalled one stops at its next read.
        # The stalled one is killed first, so the host's connection count does not grow.
        # A range still waiting for its first byte is left alone, because a queueing server (ENA waits up to 40 s) would put the new request at the back of its queue.
        # curl's own --speed-time aborts it if no byte ever comes.
        stuck = [
            byte_range
            for byte_range in self.flight
            if byte_range.start < byte_range.end
            and byte_range.receiving
            and now - byte_range.seen > TAKEOVER_S
            and self._host(byte_range.host).until <= now
        ]
        if stuck:
            old = max(stuck, key=lambda byte_range: now - byte_range.seen)
            byte_range = _Range(old.file, old.start, old.end)
            old.end = old.start
            if old.proc is not None and old.proc.poll() is None:
                old.proc.kill()
            self.takeovers += 1
            return byte_range
        return None

    def _host(self, host):
        """The state of host `host`, created on first use.
        The caller holds the lock."""
        if host not in self.hosts:
            self.hosts[host] = _Host()
        return self.hosts[host]

    def _rejected(self, host, retry_after):
        """Pause host `host` after a 429 or 503 and cut its connection limit once per pause."""
        with self.lock:
            host_state, now = self._host(host), time.time()
            busy = sum(1 for byte_range in self.flight if byte_range.host == host)
            self.rejections += 1
            if host_state.until > now:
                return  # the connections already open get the same answer
            host_state.backoff += 1
            wait = retry_after
            if wait is None:
                wait = min(REJECT_MAX_WAIT_S, BACKOFF_BASE**host_state.backoff) * (BACKOFF_JITTER_MIN + random.random())
            host_state.until = now + wait
            host_state.cap = max(1, min(host_state.cap, busy) * REJECT_KEEP_NUM // REJECT_KEEP_DEN)
            host_state.cut = now
            self.log(
                f"    fetch: {host} asked to slow down, pausing {wait:.0f} s at {host_state.cap} connections",
                flush=True,
            )

    def _wrong_answer(self, host, problem):
        """Log, once per host and kind of problem, a range answer that was discarded instead of written."""
        kind = problem.split()[0]
        with self.lock:
            if (host, kind) in self.wrong_answers:
                return
            self.wrong_answers.add((host, kind))
        self.log(f"    fetch: {host} answered a range with {problem}; discarded and retried", flush=True)

    def _answered(self, host, first, busy):
        """Lower host `host`'s limit when a range's first byte came slowly.

        `busy` is how many connections the host had open, counting this one."""
        with self.lock:
            host_state, now = self._host(host), time.time()
            if first > SLOW_FIRST_S:
                # The connections queued together report together, so one window cuts once.
                # A slow answer only suggests queueing, and on a high-latency path every answer is slow, so the cap stops at SLOW_ANSWER_MIN_CAP.
                # Only an explicit 429 or 503 may push a host lower.
                if now - host_state.cut > WINDOW_S and host_state.cap > SLOW_ANSWER_MIN_CAP:
                    host_state.cap = max(SLOW_ANSWER_MIN_CAP, min(host_state.cap, busy) - 1)
                    host_state.cut = now
            else:
                host_state.backoff = 0

    def _host_cap(self, hosts):
        """The most connections the per-host limits allow for `hosts`.  The caller holds the lock."""
        return min(self.max_conns, sum(self._host(host).cap for host in hosts))

    def _raise_hosts(self, hosts):
        """Raise the limit of each host in `hosts` that is at it and was not lowered lately.

        The caller holds the lock."""
        now = time.time()
        per = {}
        for byte_range in self.flight:
            per[byte_range.host] = per.get(byte_range.host, 0) + 1
        for host in hosts:
            host_state = self._host(host)
            if per.get(host, 0) >= host_state.cap and host_state.until <= now and now - host_state.cut > REPROBE_S:
                host_state.cap = min(self.max_conns, max(host_state.cap + CAP_RAISE_STEP, int(host_state.cap * CAP_RAISE_FACTOR)))

    def _worker(self, slot):
        """Worker loop that fetches ranges and requeues unfinished bytes, giving up after 20 failures."""
        while True:
            with self.lock:
                byte_range = self._next(slot)
                while byte_range is None:
                    self.lock.wait(1.0)
                    byte_range = self._next(slot)
                byte_range.seen = time.time()
                self.flight.append(byte_range)
                self.active += 1
                start_before = byte_range.start
            succeeded = self._run(byte_range)
            with self.lock:
                self.flight.remove(byte_range)
                self.active -= 1
                if (
                    byte_range.start < byte_range.end and not byte_range.file.err
                ):  # bytes left, so the range goes back to the queue in file order
                    # Only attempts in a row without a byte count, so a flaky link that keeps making progress never fails the file.
                    if byte_range.rejected:
                        byte_range.rejects += 1
                    elif byte_range.start > start_before:
                        byte_range.attempts = 0
                    elif not succeeded:
                        byte_range.attempts += 1
                    if byte_range.attempts > RANGE_MAX_ATTEMPTS or byte_range.rejects > REJECTS_MAX:
                        byte_range.file.err = (
                            f"range {byte_range.start}-{byte_range.end} failed {RANGE_MAX_ATTEMPTS} times in a row without a byte"
                            if byte_range.attempts > RANGE_MAX_ATTEMPTS
                            else f"range {byte_range.start}-{byte_range.end} rejected {REJECTS_MAX} times (429/503)"
                        )
                        byte_range.file.event.set()
                    else:
                        self.queue.append(byte_range)
                        self.queue.sort(key=lambda queued: (queued.file.order, queued.start))
                self.lock.notify_all()
            # A rejected range waits through its host's pause instead, in _next.
            if not succeeded and not byte_range.rejected:
                time.sleep(min(RETRY_SLEEP_MAX_S, RETRY_SLEEP_PER_ATTEMPT_S * byte_range.attempts))

    def _run(self, byte_range):
        """Stream `byte_range` into its file, returning False on a curl failure or a stall."""
        file_state = byte_range.file
        if file_state.err:
            byte_range.start = byte_range.end
            return True
        with self.lock:
            range_lo, range_hi = byte_range.start, byte_range.end
            busy = sum(1 for other in self.flight if other.host == byte_range.host)
            byte_range.receiving = False
        byte_range.rejected = False
        start_time, first = time.time(), None
        status = retry_after = None
        # -D - puts the headers ahead of the body on stdout, so the response is checked before a byte reaches the file.
        proc = subprocess.Popen(
            CURL + ["-D", "-", "-r", f"{range_lo}-{range_hi - 1}", file_state.url],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            bufsize=0,
        )
        with self.lock:
            byte_range.proc = proc
            if byte_range.start >= byte_range.end:  # taken over while starting
                proc.kill()
        try:
            head, chunk = _read_head(proc.stdout)
            status, retry_after = _response(head)
            if status in REJECT_CODES:
                byte_range.rejected = True
                return False
            if chunk is None:
                return False
            first = time.time() - start_time
            served = _served_size(head)
            if served is not None and served != file_state.size:
                # The server's copy changed size, so no number of retries can produce the published file.
                with self.lock:
                    file_state.err = f"server now serves {served} bytes, expected {file_state.size}"
                    byte_range.start = byte_range.end
                file_state.event.set()
                return True
            problem = _range_problem(head, range_lo)
            if problem:
                self._wrong_answer(byte_range.host, problem)
                return False
            while not file_state.err:
                # The read that ended the headers may hold no body yet.
                if not chunk:
                    chunk = proc.stdout.read(CURL_READ_BYTES)
                    if not chunk:
                        break
                # The bytes are claimed before they are written.
                # A steal moves end but never below start, so the two halves stay disjoint.
                with self.lock:
                    start = byte_range.start
                    count = min(len(chunk), byte_range.end - start)
                    byte_range.start += count
                    byte_range.seen = time.time()
                    byte_range.receiving = True
                    self.bytes += count
                # pwrite takes an explicit offset, so parallel workers never share a file position.
                os.pwrite(file_state.part_fd, chunk[:count], start)
                with file_state.lock:
                    file_state.done.append([start, start + count])
                    file_state.left -= count
                    if len(file_state.done) > DONE_RANGES_MERGE_AT:
                        file_state.done = _merge(file_state.done)
                    finished = file_state.left == 0
                if finished:
                    file_state.event.set()
                if byte_range.start >= byte_range.end:
                    break
                chunk = b""
        finally:
            if proc.poll() is None:
                proc.kill()
            proc.wait()
            # Without the close, each attempt would keep its pipe descriptor until garbage collection, and long fetches retry thousands of times.
            proc.stdout.close()
            # A 429 or 503 marks the host as throttling.
            # No answer at all until curl's speed limit gave up is how a queueing server looks, so it lowers the host's limit like a slow first byte.
            if byte_range.rejected:
                self._rejected(byte_range.host, retry_after)
            else:
                self._answered(byte_range.host, time.time() - start_time if first is None else first, busy)
            if TRACE:
                with open(TRACE, "a") as handle:
                    handle.write(
                        f"{start_time:.1f}\t{byte_range.host}\t{range_lo}\t{range_hi}\t{byte_range.start}\t{proc.returncode}\t"
                        f"{NO_FIRST_BYTE if first is None else round(first, 2)}\t{time.time() - start_time:.1f}\n"
                    )
        return byte_range.start >= byte_range.end

    def _syncer(self):
        """Every few seconds, fdatasync each file and then record what is on disk."""
        while True:
            time.sleep(SYNC_INTERVAL_S)
            with self.lock:
                files = list(self.files)
            for file_state in files:
                with file_state.lock:
                    if file_state.closed:
                        continue
                    snap = _merge(file_state.done)
                    os.fdatasync(file_state.part_fd)  # snap is on disk from here on
                    self._save(file_state.part, snap)

    def _write_progress(self, rate, busy):
        """Write the progress file, replaced atomically.

        It holds the time, the rate over the last window, the busy and targeted connection counts and every file in flight."""
        if not self.progress:
            return
        with self.lock:
            files = [
                {"part": file_state.part, "size": file_state.size, "left": file_state.left} for file_state in self.files
            ]
            target = self.target
        progress = {
            "time": time.time(),
            "pid": os.getpid(),
            "rate_bytes_s": rate,
            "connections": busy,
            "target": target,
            "files": files,
        }
        tmp = self.progress + ".tmp"
        try:
            with open(tmp, "w") as handle:
                json.dump(progress, handle)
            os.replace(tmp, self.progress)
        except OSError:
            pass  # progress is informational, so it never fails a fetch

    def _control(self):
        """Choose the connection count from measured totals.

        The choice is the smallest count within GAIN of the best rate seen, and when that is the largest count tried, a higher one is probed.
        Old measurements decay, and every REPROBE_S the neighbouring counts are measured again."""
        best = {}  # n -> rate
        last_probe = last_log = time.time()
        while True:
            with self.lock:
                bytes_before = self.bytes
            time.sleep(WINDOW_S)
            with self.lock:
                rate = (self.bytes - bytes_before) / WINDOW_S
                conns, busy = self.target, self.active
                pending = len(self.queue) + len(self.flight)
                left = sum(file_state.left for file_state in self.files)
                hosts = {byte_range.host for byte_range in self.queue + self.flight}
                cap = self._host_cap(hosts)
            self._write_progress(rate, busy)
            if pending == 0:
                continue
            conns = min(conns, cap)
            if busy >= conns:  # only a window with n busy connections measures n
                best[conns] = max(rate, BEST_RATE_DECAY * best.get(conns, 0.0))
                peak = max(best.values())
                good = min(conns_tried for conns_tried in best if best[conns_tried] >= peak / (1 + GAIN))
                if time.time() - last_probe > REPROBE_S:
                    last_probe = time.time()
                    best = {good: best[good]}
                if good == max(best) and good < self.max_conns:
                    if good >= cap:
                        # The total still rises at the host limits, so the probe raises them.
                        with self.lock:
                            self._raise_hosts(hosts)
                            cap = self._host_cap(hosts)
                    good = min(cap, max(good + CAP_RAISE_STEP, int(good * CAP_RAISE_FACTOR)))
                conns = good
            if conns != self.target:
                with self.lock:
                    self.target = conns
                    self.lock.notify_all()
            if time.time() - last_log >= LOG_S:
                last_log = time.time()
                self.log(
                    f"    fetch: {rate / BYTES_PER_MB:.1f} MB/s over {busy} connections "
                    f"(target {self.target}), {left / BYTES_PER_GB:.1f} GB left; "
                    f"{self.steals} splits, {self.takeovers} takeovers, "
                    f"{self.rejections} rejections",
                    flush=True,
                )


_SCHEDULER = None
_S_LOCK = threading.Lock()


def scheduler(max_conns=None, log=print, progress=None):
    """The process's single scheduler.

    `progress` names a JSON file rewritten every WINDOW_S with the measured rate, the connections and each file's bytes left, which `compare_tools_published_sets.py status` reads."""
    global _SCHEDULER
    with _S_LOCK:
        if _SCHEDULER is None:
            _SCHEDULER = Scheduler(max_conns or MAX_CONNS, log, progress)
        return _SCHEDULER
