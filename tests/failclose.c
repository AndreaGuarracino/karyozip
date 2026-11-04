/* LD_PRELOAD shim for run.sh and diff.sh: a write error reported late, as NFS, a quota or a failing disk may report
   one only at fsync or close. A descriptor whose path (read from /proc/self/fd) contains $KZ_FAIL_FD_MATCH fails
   fsync, fdatasync and close with EIO; KZ_FAIL_FD_OP=sync or close fails only that one. close still closes, as Linux
   does when it reports an error. Only calls through the symbols are seen: stdio's fclose closes inside libc.
     cc -shared -fPIC -o failclose.so failclose.c -ldl */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int hit(int fd, const char *op) {
    const char *m = getenv("KZ_FAIL_FD_MATCH"), *o = getenv("KZ_FAIL_FD_OP");
    if (!m || !*m || (o && *o && strcmp(o, op))) return 0;
    char l[64], p[4096];
    snprintf(l, sizeof l, "/proc/self/fd/%d", fd);
    ssize_t n = readlink(l, p, sizeof p - 1);
    if (n <= 0) return 0;
    p[n] = 0;
    return strstr(p, m) != NULL;
}

int close(int fd) {
    static int (*real)(int);
    if (!real) real = (int (*)(int))dlsym(RTLD_NEXT, "close");
    int h = hit(fd, "close"), r = real(fd);
    if (h) { errno = EIO; return -1; }
    return r;
}
int fsync(int fd) {
    static int (*real)(int);
    if (!real) real = (int (*)(int))dlsym(RTLD_NEXT, "fsync");
    if (hit(fd, "sync")) { errno = EIO; return -1; }
    return real(fd);
}
int fdatasync(int fd) {
    static int (*real)(int);
    if (!real) real = (int (*)(int))dlsym(RTLD_NEXT, "fdatasync");
    if (hit(fd, "sync")) { errno = EIO; return -1; }
    return real(fd);
}
