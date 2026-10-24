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
    Near the end, a range without progress for TAKEOVER_S seconds is handed to an idle connection.
  * learns each host's connection limit.
    A host starts at HOST_START connections.
    When the measured total is still rising at the limit, the probe for more connections raises the limit of every host at it, so a server that takes more connections gets them, and the measurement then keeps them only if they pay.
    A first byte slower than SLOW_FIRST_S, the sign of a server queueing connections, lowers the host's limit below the connections it had open, and the host is not raised again for REPROBE_S.
    A 429 or 503 answer cuts the host's limit to three quarters of the connections it had open and pauses the host for its Retry-After, or for a growing randomized backoff without one.
    A rejection is not a failure of the range, which waits and retries, so a busy server slows the fetch without failing it.

Bytes are written in place (pwrite) into a preallocated <dest>.part.
The completed intervals are saved to <dest>.part.done only after an fdatasync, so a restart resumes from what has reached the disk.
A plain .part prefix left by a single-connection curl, and the .part.N/.part.ranges segments of the earlier fixed-split fetch, are taken over as completed intervals.

Only curl and the Python standard library are used, so it runs on any Linux."""

import json
import os
import random
import re
import subprocess
import threading
import time

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
REJECT_CODES = ("429", "503")  # the server asks the client to slow down
REJECT_MAX_WAIT_S = 120.0  # longest pause of a host without a Retry-After
REJECT_KEEP_NUM, REJECT_KEEP_DEN = 3, 4  # a rejection keeps 3/4 of the open connections
REJECTS_MAX = 200  # rejections of one range before its file fails
TRACE = os.environ.get("RANGEFETCH_TRACE")  # file that gets one line per range attempt
CURL = ["curl", "-fsSL", "--speed-limit", "10000", "--speed-time", str(STALL_S)]


def _merge(iv):
    """Sorted union of the [s, e) intervals iv, with overlapping and touching ones joined."""
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def _gaps(done, size):
    """The [s, e) intervals of [0, size) not covered by `done`."""
    out, pos = [], 0
    for s, e in _merge(done):
        if s > pos:
            out.append([pos, s])
        pos = max(pos, e)
    if pos < size:
        out.append([pos, size])
    return out


def supports_ranges(url):
    """True when the server answers a one-byte range with 206.
    A server that ignores Range would send the whole file to every connection."""
    r = subprocess.run(
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
    return r.stdout.strip().endswith("206")


def content_length(url):
    """Content-Length of url after redirects, or None when the server gives none."""
    r = subprocess.run(
        ["curl", "-sSIL", "--retry", "4", "--retry-all-errors", url],
        capture_output=True,
        text=True,
        check=False,
    )
    n = None
    for line in r.stdout.splitlines():
        if line.lower().startswith("content-length:"):
            n = int(
                line.split(":", 1)[1]
            )  # the last one belongs to the final redirect target
    return n


def _response(hdr):
    """(HTTP status, Retry-After seconds or None) of the last response in a curl header dump."""
    try:
        with open(hdr, errors="replace") as fh:
            text = fh.read()
    except OSError:
        return None, None
    finally:
        try:
            os.remove(hdr)
        except OSError:
            pass
    blocks = [b for b in re.split(r"\r?\n\r?\n", text) if b.strip()]
    if not blocks:
        return None, None
    last = blocks[-1]
    m = re.match(r"HTTP/\S+\s+(\d{3})", last)
    ra = re.search(r"(?im)^retry-after:\s*(\d+)\s*$", last)
    return (m.group(1) if m else None), (float(ra.group(1)) if ra else None)


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
        self.done = []  # completed [s, e) intervals, in memory
        self.synced = []  # the part of done known to be on disk
        self.left = 0  # bytes still to write
        self.err = None
        self.event = threading.Event()
        self.fd = None
        self.closed = False


class _Range:
    """A byte range [s, e) of a file assigned to one connection."""

    def __init__(self, f, s, e):
        self.f, self.s, self.e = (
            f,
            s,
            e,
        )  # s advances as bytes arrive, e may shrink when stolen from
        self.attempts = 0
        self.rejects = 0
        self.rejected = False
        self.seen = time.time()  # time of the last progress
        self.proc = None
        self.host = f.url.split("/")[2]


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

    # Files
    def fetch(self, url, part, size):
        """Fill the size bytes of `part` from url, blocking until done and raising on failure."""
        with self.lock:
            f = _File(url, part, size, self.order)
            self.order += 1
        done = self._adopt(part, size)
        f.fd = os.open(part, os.O_RDWR | os.O_CREAT, 0o644)
        os.ftruncate(f.fd, size)
        f.done, f.synced = [list(x) for x in done], [list(x) for x in done]
        gaps = _gaps(done, size)
        f.left = sum(e - s for s, e in gaps)
        if f.left < size:
            self.log(
                f"    resuming {os.path.basename(part)}: {size - f.left} of {size} bytes already on disk",
                flush=True,
            )
        if f.left == 0:
            os.close(f.fd)
            self._clear_state(part)
            return
        chunk = max(CHUNK, min(CHUNK_MAX, size // 16))
        with self.lock:
            for s, e in gaps:
                for c in range(s, e, chunk):
                    self.queue.append(_Range(f, c, min(c + chunk, e)))
            self.queue.sort(key=lambda r: (r.f.order, r.s))
            self.files.append(f)
            self._start()
            self.lock.notify_all()
        f.event.wait()
        with self.lock:
            # A file failed, so its queued ranges are dropped and its in-flight ones stop.
            self.queue = [r for r in self.queue if r.f is not f]
            while any(r.f is f for r in self.flight):
                self.lock.wait(1.0)
            self.files.remove(f)
        with f.lock:
            f.closed = True
            os.fsync(f.fd)
            os.close(f.fd)
        if f.err:
            raise SystemExit(f"{part}: {f.err}")
        self._clear_state(part)

    def _adopt(self, part, size):
        """Completed intervals of `part` that are already on disk."""
        state = part + ".done"
        if os.path.exists(state):
            with open(state) as fh:
                return json.load(fh)
        done = []
        if os.path.exists(part):  # a prefix from a single-connection download
            done.append([0, min(os.path.getsize(part), size)])
        plan = part + ".ranges"  # fixed-split segments left by the earlier fetch
        if os.path.exists(plan):
            with open(plan) as fh:
                bounds = json.load(fh)
            fd = os.open(part, os.O_RDWR | os.O_CREAT, 0o644)
            os.ftruncate(fd, size)
            for n, (s, e) in enumerate(bounds):
                seg = f"{part}.{n}"
                if not os.path.exists(seg):
                    continue
                with open(seg, "rb") as sf:
                    pos = s
                    for b in iter(lambda: sf.read(1 << 24), b""):
                        b = b[: e + 1 - pos]
                        os.pwrite(fd, b, pos)
                        pos += len(b)
                done.append([s, pos])
            os.fsync(fd)
            os.close(fd)
            self._save(part, done)
            for n in range(len(bounds)):
                if os.path.exists(f"{part}.{n}"):
                    os.remove(f"{part}.{n}")
            os.remove(plan)
        return _merge(done)

    def _save(self, part, done):
        """Atomically replace part's .done file with the merged intervals."""
        tmp = part + ".done.tmp"
        with open(tmp, "w") as fh:
            json.dump(_merge(done), fh)
        os.replace(tmp, part + ".done")

    def _clear_state(self, part):
        """Remove part's interval files once the download is complete."""
        for p in (part + ".done", part + ".done.tmp"):
            if os.path.exists(p):
                os.remove(p)

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
        for r in self.flight:
            per[r.host] = per.get(r.host, 0) + 1
        now = time.time()

        def room(h):
            """True when host h is not paused and is below its learned limit."""
            hs = self._host(h)
            return hs.until <= now and per.get(h, 0) < hs.cap

        for n, r in enumerate(self.queue):
            if room(r.host):
                return self.queue.pop(n)
        open_ = {r.host for r in self.flight if room(r.host)}
        # Work stealing takes the far half of the largest range in flight.
        big = max(
            (r for r in self.flight if r.host in open_),
            key=lambda r: r.e - r.s,
            default=None,
        )
        if big is not None and big.e - big.s >= 2 * MIN_SPLIT:
            mid = big.s + (big.e - big.s) // 2
            r = _Range(big.f, mid, big.e)
            big.e = mid
            self.steals += 1
            return r
        # In the endgame, a range without progress for TAKEOVER_S moves whole to this idle connection, and the stalled one stops at its next read.
        # The stalled one is killed first, so the host's connection count does not grow.
        stuck = [
            r
            for r in self.flight
            if r.s < r.e
            and now - r.seen > TAKEOVER_S
            and self._host(r.host).until <= now
        ]
        if stuck:
            old = max(stuck, key=lambda r: now - r.seen)
            r = _Range(old.f, old.s, old.e)
            old.e = old.s
            if old.proc is not None and old.proc.poll() is None:
                old.proc.kill()
            self.takeovers += 1
            return r
        return None

    def _host(self, h):
        """The state of host h, created on first use.
        The caller holds the lock."""
        if h not in self.hosts:
            self.hosts[h] = _Host()
        return self.hosts[h]

    def _rejected(self, h, retry_after):
        """Pause host h after a 429 or 503 and cut its connection limit once per pause."""
        with self.lock:
            hs, now = self._host(h), time.time()
            busy = sum(1 for r in self.flight if r.host == h)
            self.rejections += 1
            if hs.until > now:
                return  # the connections already open get the same answer
            hs.backoff += 1
            wait = retry_after
            if wait is None:
                wait = min(REJECT_MAX_WAIT_S, 2.0**hs.backoff) * (0.5 + random.random())
            hs.until = now + wait
            hs.cap = max(1, min(hs.cap, busy) * REJECT_KEEP_NUM // REJECT_KEEP_DEN)
            hs.cut = now
            self.log(
                f"    fetch: {h} asked to slow down, pausing {wait:.0f} s at {hs.cap} connections",
                flush=True,
            )

    def _answered(self, h, first, busy):
        """Lower host h's limit when a range's first byte came slowly.

        `busy` is how many connections the host had open, counting this one."""
        with self.lock:
            hs, now = self._host(h), time.time()
            if first > SLOW_FIRST_S:
                # The connections queued together report together, so one window cuts once.
                if now - hs.cut > WINDOW_S:
                    hs.cap = max(1, min(hs.cap, busy) - 1)
                    hs.cut = now
            else:
                hs.backoff = 0

    def _host_cap(self, hosts):
        """The most connections the per-host limits allow for `hosts`.  The caller holds the lock."""
        return min(self.max_conns, sum(self._host(h).cap for h in hosts))

    def _raise_hosts(self, hosts):
        """Raise the limit of each host in `hosts` that is at it and was not lowered lately.

        The caller holds the lock."""
        now = time.time()
        per = {}
        for r in self.flight:
            per[r.host] = per.get(r.host, 0) + 1
        for h in hosts:
            hs = self._host(h)
            if per.get(h, 0) >= hs.cap and hs.until <= now and now - hs.cut > REPROBE_S:
                hs.cap = min(self.max_conns, max(hs.cap + 2, int(hs.cap * 1.5)))

    def _worker(self, slot):
        """Worker loop that fetches ranges and requeues unfinished bytes, giving up after 20 failures."""
        while True:
            with self.lock:
                r = self._next(slot)
                while r is None:
                    self.lock.wait(1.0)
                    r = self._next(slot)
                r.seen = time.time()
                self.flight.append(r)
                self.active += 1
            ok = self._run(r)
            with self.lock:
                self.flight.remove(r)
                self.active -= 1
                if (
                    r.s < r.e
                ):  # bytes left, so the range goes back to the queue in file order
                    if r.rejected:
                        r.rejects += 1
                    else:
                        r.attempts = r.attempts + (0 if ok else 1)
                    if r.attempts > 20 or r.rejects > REJECTS_MAX:
                        r.f.err = (
                            f"range {r.s}-{r.e} failed 20 times"
                            if r.attempts > 20
                            else f"range {r.s}-{r.e} rejected {REJECTS_MAX} times (429/503)"
                        )
                        r.f.event.set()
                    else:
                        self.queue.append(r)
                        self.queue.sort(key=lambda x: (x.f.order, x.s))
                self.lock.notify_all()
            # A rejected range waits through its host's pause instead, in _next.
            if not ok and not r.rejected:
                time.sleep(min(30, 2 * r.attempts))

    def _run(self, r):
        """Stream range r into its file, returning False on a curl failure or a stall."""
        f = r.f
        if f.err:
            r.s = r.e
            return True
        with self.lock:
            lo, hi = r.s, r.e
            busy = sum(1 for x in self.flight if x.host == r.host)
        r.rejected = False
        t0, first = time.time(), None
        hdr = f"{f.part}.hdr.{threading.get_ident()}"
        p = subprocess.Popen(
            CURL + ["-D", hdr, "-r", f"{lo}-{hi - 1}", f.url],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            bufsize=0,
        )
        with self.lock:
            r.proc = p
            if r.s >= r.e:  # taken over while starting
                p.kill()
        try:
            while not f.err:
                b = p.stdout.read(1 << 20)
                if not b:
                    break
                if first is None:
                    first = time.time() - t0
                # The bytes are claimed before they are written.
                # A steal moves r.e but never below r.s, so the two halves stay disjoint.
                with self.lock:
                    start = r.s
                    n = min(len(b), r.e - start)
                    r.s += n
                    r.seen = time.time()
                    self.bytes += n
                os.pwrite(f.fd, b[:n], start)
                with f.lock:
                    f.done.append([start, start + n])
                    f.left -= n
                    if len(f.done) > 64:
                        f.done = _merge(f.done)
                    finished = f.left == 0
                if finished:
                    f.event.set()
                if r.s >= r.e:
                    break
        finally:
            if p.poll() is None:
                p.kill()
            p.wait()
            status, retry_after = _response(hdr)
            if first is None and status in REJECT_CODES:
                r.rejected = True
                self._rejected(r.host, retry_after)
            elif first is not None:
                self._answered(r.host, first, busy)
            if TRACE:
                with open(TRACE, "a") as fh:
                    fh.write(
                        f"{t0:.1f}\t{r.host}\t{lo}\t{hi}\t{r.s}\t{p.returncode}\t"
                        f"{-1 if first is None else round(first, 2)}\t{time.time() - t0:.1f}\n"
                    )
        return r.s >= r.e

    def _syncer(self):
        """Every few seconds, fdatasync each file and then record what is on disk."""
        while True:
            time.sleep(5)
            with self.lock:
                files = list(self.files)
            for f in files:
                with f.lock:
                    if f.closed:
                        continue
                    snap = _merge(f.done)
                    os.fdatasync(f.fd)  # snap is on disk from here on
                    self._save(f.part, snap)

    def _write_progress(self, rate, busy):
        """Write the progress file, replaced atomically.

        It holds the time, the rate over the last window, the busy and targeted connection counts and every file in flight."""
        if not self.progress:
            return
        with self.lock:
            files = [
                {"part": f.part, "size": f.size, "left": f.left} for f in self.files
            ]
            target = self.target
        d = {
            "time": time.time(),
            "pid": os.getpid(),
            "rate_bytes_s": rate,
            "connections": busy,
            "target": target,
            "files": files,
        }
        tmp = self.progress + ".tmp"
        try:
            with open(tmp, "w") as o:
                json.dump(d, o)
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
                b0 = self.bytes
            time.sleep(WINDOW_S)
            with self.lock:
                rate = (self.bytes - b0) / WINDOW_S
                n, busy = self.target, self.active
                pending = len(self.queue) + len(self.flight)
                left = sum(f.left for f in self.files)
                hosts = {r.host for r in self.queue + self.flight}
                cap = self._host_cap(hosts)
            self._write_progress(rate, busy)
            if pending == 0:
                continue
            n = min(n, cap)
            if busy >= n:  # only a window with n busy connections measures n
                best[n] = max(rate, 0.7 * best.get(n, 0.0))
                peak = max(best.values())
                good = min(m for m in best if best[m] >= peak / (1 + GAIN))
                if time.time() - last_probe > REPROBE_S:
                    last_probe = time.time()
                    best = {good: best[good]}
                if good == max(best) and good < self.max_conns:
                    if good >= cap:
                        # The total still rises at the host limits, so the probe raises them.
                        with self.lock:
                            self._raise_hosts(hosts)
                            cap = self._host_cap(hosts)
                    good = min(cap, max(good + 2, int(good * 1.5)))
                n = good
            if n != self.target:
                with self.lock:
                    self.target = n
                    self.lock.notify_all()
            if time.time() - last_log >= LOG_S:
                last_log = time.time()
                self.log(
                    f"    fetch: {rate / 1e6:.1f} MB/s over {busy} connections "
                    f"(target {self.target}), {left / 1e9:.1f} GB left; "
                    f"{self.steals} splits, {self.takeovers} takeovers, "
                    f"{self.rejections} rejections",
                    flush=True,
                )


_S = None
_S_LOCK = threading.Lock()


def scheduler(max_conns=None, log=print, progress=None):
    """The process's single scheduler.

    `progress` names a JSON file rewritten every WINDOW_S with the measured rate, the connections and each file's bytes left, which `compare_tools_published_sets.py status` reads."""
    global _S
    with _S_LOCK:
        if _S is None:
            _S = Scheduler(max_conns or MAX_CONNS, log, progress)
        return _S
