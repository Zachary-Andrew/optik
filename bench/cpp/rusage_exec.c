/* rusage_exec.c - run a command and report its own peak RSS.
 *
 *   rusage_exec [-a] [-o OUT] [-e ERR] AS_LIMIT_BYTES RESULT_FILE cmd [args...]
 *
 * Forks from this tiny process (so the child's ru_maxrss is not inflated by a large parent, since Linux folds the pre-exec address space into the exec'd task's maxrss), applies an address-space limit (0 = none), execs the command and writes "exit_code wall_seconds maxrss_kb written_bytes read_bytes" to RESULT_FILE.
 * The I/O fields are the block I/O of the command and the children it reaped (ru_oublock / ru_inblock in 512-byte units, i.e. the bytes it sent to and took from storage, as Cuttlefish 2 and 3 report disk use).
 * Exits with the command's status.
 *
 *   -a      append the line instead of overwriting.
 *           Under hyperfine every run, warmups first, adds one line, and a bench keeps the last --runs lines and drops the warmups'
 *   -o OUT  the command's stdout to OUT, -e ERR its stderr to ERR (each
 *           truncated per run, so after hyperfine they hold the last run's) */
#define _GNU_SOURCE
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

static void redirect(const char* path, int fd) {
    const int f = open(path, O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
    if (f < 0) { perror(path); _exit(127); }
    if (dup2(f, fd) < 0) { perror("dup2"); _exit(127); }
}

int main(int argc, char** argv) {
    int append = 0, a = 1;
    const char* out = NULL;
    const char* err = NULL;
    while (a < argc && argv[a][0] == '-' && argv[a][1] != '\0') {
        if (strcmp(argv[a], "-a") == 0) { append = 1; a += 1; }
        else if (strcmp(argv[a], "-o") == 0 && a + 1 < argc) { out = argv[a + 1]; a += 2; }
        else if (strcmp(argv[a], "-e") == 0 && a + 1 < argc) { err = argv[a + 1]; a += 2; }
        else break;
    }
    if (argc - a < 3) {
        fprintf(stderr, "usage: rusage_exec [-a] [-o OUT] [-e ERR] AS_LIMIT RESULT cmd [args]\n");
        return 2;
    }
    const unsigned long long lim = strtoull(argv[a], NULL, 10);
    const char* result = argv[a + 1];
    char** cmd = argv + a + 2;
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    pid_t pid = fork();
    if (pid < 0) { perror("fork"); return 2; }
    if (pid == 0) {
        if (out) redirect(out, 1);
        if (err) redirect(err, 2);
        if (lim) { struct rlimit r = {lim, lim}; setrlimit(RLIMIT_AS, &r); }
        execvp(cmd[0], cmd);
        perror("execvp");
        _exit(127);
    }
    int st = 0;
    struct rusage ru;
    if (wait4(pid, &st, 0, &ru) < 0) { perror("wait4"); return 2; }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    const int code = WIFEXITED(st) ? WEXITSTATUS(st) : 128 + WTERMSIG(st);
    FILE* f = fopen(result, append ? "a" : "w");
    if (f) {
        fprintf(f, "%d %.6f %ld %lld %lld\n", code, (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9,
                ru.ru_maxrss, (long long)ru.ru_oublock * 512LL, (long long)ru.ru_inblock * 512LL);
        fclose(f);
    }
    return code;
}
