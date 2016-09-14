
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"

lulog *dbg = NULL;

int unpack_c(const char *pattern, int **offsets, int *length) {
    return LU_OK;
}

int unpack_b(const char *pattern, int **offsets, int *length) {
    return LU_OK;
}

int unpack_a(const char *pattern, int **offsets, int *length) {
    return LU_OK;
}

int unpack(const char *pattern, int **offsets, int *length) {

    LU_STATUS

    if (strchr(pattern, 'A')) {
        LU_CHECK(unpack_a(pattern, offsets, length));
    } else if (strchr(pattern, 'B')) {
        LU_CHECK(unpack_b(pattern, offsets, length));
    } else if (strchr(pattern, 'C')) {
        LU_CHECK(unpack_c(pattern, offsets, length));
    } else {
        luerror(dbg, "Did not find group type (A, B, C) in %s", pattern);
        status = LU_ERR_ARG;
    }

    LU_NO_CLEANUP
}

int plot(const char *pattern) {

    LU_STATUS;
    int *offsets = NULL, length = 0;

    luinfo(dbg, "Pattern '%s'", pattern);
    LU_CHECK(unpack(pattern, &offsets, &length));

LU_CLEANUP
    free(offsets);
    LU_RETURN
}

void usage(const char *progname) {
    luinfo(dbg, "Plot the given spoke pattern");
    luinfo(dbg, "%s -h        display this message", progname);
    luinfo(dbg, "%s pattern   plot pattern to pattern.png", progname);
}

// error handling is for lulib routines; don't bother elsewhere.
int main(int argc, char** argv) {

    LU_STATUS

    lulog_mkstderr(&dbg, lulog_level_debug);
    if (argc != 2 || !strcmp("-h", argv[1])) {
        usage(argv[0]);
    } else {
        plot(argv[1]);
    }

LU_CLEANUP
    if (dbg) status = dbg->free(&dbg, status);
    return status;
}
