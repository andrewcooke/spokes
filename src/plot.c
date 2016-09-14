
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

    LU_STATUS

    const char *p = pattern;
    int sign = 1;

    while (*p != 'A') {
        if (*p == ',') sign = 1;
        else if (*p == '-') sign = -sign;
        else {
            int offset = *p - '0';
            if (!(*offsets = realloc(*offsets, (1 + *length) * sizeof(**offsets)))) return LU_ERR_MEM;
            (*offsets)[*length] = sign * offset;
            (*length)++;
            sign = 1;
        }
        p++;
    }
    p++;   // drop A
    if ((*offsets)[*length-1]) luwarn(dbg, "Central offset for group A is non-zero");
    if (!(*offsets = realloc(*offsets, (2 * *length - 1) * sizeof(**offsets)))) return LU_ERR_MEM;
    for (int i = 0; i < *length; ++i) (*offsets)[*length + i - 1] = (*offsets)[*length - i - 1];
    *length = 2 * *length - 1;
    if (*p) {
        int padding = *p - '0';
        if (!(*offsets = realloc(*offsets, (padding + *length) * sizeof(**offsets)))) return LU_ERR_MEM;
        for (int i = 0; i < padding; ++i) (*offsets)[*length + i] = 0;
        *length = padding + *length;
    }

    LU_NO_CLEANUP
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

int dump_pattern(int *offsets, int length) {

    LU_STATUS
    char *buffer = NULL, *p;

    LU_ALLOC(dbg, buffer, length * 3 + 4);
    p = buffer;
    for (int i = 0; i < length; ++i) {
        if (i) p += sprintf(p, ",");
        p += sprintf(p, "%d", offsets[i]);
    }
    *p = '\0';

    luinfo(dbg, "Pattern: %s (length %d)", buffer, length);

LU_CLEANUP
    free(buffer);
    LU_RETURN
}

int plot(const char *pattern) {

    LU_STATUS;
    int *offsets = NULL, length = 0;

    luinfo(dbg, "Pattern '%s'", pattern);
    LU_CHECK(unpack(pattern, &offsets, &length));
    LU_CHECK(dump_pattern(offsets, length));

LU_CLEANUP
    free(offsets);
    LU_RETURN
}

void usage(const char *progname) {
    luinfo(dbg, "Plot the given spoke pattern");
    luinfo(dbg, "%s -h        display this message", progname);
    luinfo(dbg, "%s pattern   plot pattern to pattern.png", progname);
    luinfo(dbg, "(file name has commas removed)", progname);
}

// error handling is for lulib routines; don't bother elsewhere.
int main(int argc, char** argv) {

    LU_STATUS

    lulog_mkstderr(&dbg, lulog_level_debug);
    if (argc != 2 || !strcmp("-h", argv[1])) {
        usage(argv[0]);
    } else {
        LU_CHECK(plot(argv[1]));
    }

LU_CLEANUP
    if (dbg) status = dbg->free(&dbg, status);
    return status;
}
