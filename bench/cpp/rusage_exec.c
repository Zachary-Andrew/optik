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
/* Mode of the -o / -e output files. */
enum { OUTPUT_FILE_MODE = 0644 };
/* Exit status when the command cannot be started, and the base added to a terminating signal, as shells report them. */
enum { EXIT_EXEC_FAILED = 127, EXIT_SIGNAL_BASE = 128 };
/* Positional arguments: the address-space limit, the result file and the command. */
enum { MIN_POSITIONAL_ARGS = 3 };
/* ru_oublock / ru_inblock count 512-byte blocks. */
enum { BLOCK_BYTES = 512 };
/* Nanoseconds per second. */
#define NS_PER_SEC 1e9


static void redirect(const char* path, int target_fd) {
    const int src_fd = open(path, O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, OUTPUT_FILE_MODE);
    if (src_fd < 0) { perror(path); _exit(EXIT_EXEC_FAILED); }
    if (dup2(src_fd, target_fd) < 0) { perror("dup2"); _exit(EXIT_EXEC_FAILED); }
}

int main(int argc, char** argv) {
    int append = 0, arg_idx = 1;
    const char* out = NULL;
    const char* err = NULL;
    while (arg_idx < argc && argv[arg_idx][0] == '-' && argv[arg_idx][1] != '\0') {
        if (strcmp(argv[arg_idx], "-a") == 0) { append = 1; arg_idx += 1; }
        else if (strcmp(argv[arg_idx], "-o") == 0 && arg_idx + 1 < argc) { out = argv[arg_idx + 1]; arg_idx += 2; }
        else if (strcmp(argv[arg_idx], "-e") == 0 && arg_idx + 1 < argc) { err = argv[arg_idx + 1]; arg_idx += 2; }
        else break;
    }
    if (argc - arg_idx < MIN_POSITIONAL_ARGS) {
        fprintf(stderr, "usage: rusage_exec [-a] [-o OUT] [-e ERR] AS_LIMIT RESULT cmd [args]\n");
        return 2;
    }
    const unsigned long long lim = strtoull(argv[arg_idx], NULL, 10);
    const char* result = argv[arg_idx + 1];
    char** cmd = argv + arg_idx + 2;
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    pid_t pid = fork();
    if (pid < 0) { perror("fork"); return 2; }
    if (pid == 0) {
        if (out) redirect(out, 1);
        if (err) redirect(err, 2);
        if (lim) { struct rlimit limit_spec = {lim, lim}; setrlimit(RLIMIT_AS, &limit_spec); }
        execvp(cmd[0], cmd);
        perror("execvp");
        _exit(EXIT_EXEC_FAILED);
    }
    int status = 0;
    struct rusage usage;
    if (wait4(pid, &status, 0, &usage) < 0) { perror("wait4"); return 2; }
    clock_gettime(CLOCK_MONOTONIC, &end);
    const int code = WIFEXITED(status) ? WEXITSTATUS(status) : EXIT_SIGNAL_BASE + WTERMSIG(status);
    FILE* result_file = fopen(result, append ? "a" : "w");
    if (result_file) {
        fprintf(result_file, "%d %.6f %ld %lld %lld\n", code, (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / NS_PER_SEC,
                usage.ru_maxrss, (long long)usage.ru_oublock * (long long)BLOCK_BYTES, (long long)usage.ru_inblock * (long long)BLOCK_BYTES);
        fclose(result_file);
    }
    return code;
}
