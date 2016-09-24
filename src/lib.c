
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <cairo/cairo.h>

#include "lu/status.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"

#include "lib.h"


void draw_circle(cairo_t *cr, float r) {
    cairo_move_to(cr, r, 0);
    cairo_arc(cr, 0, 0, r, 0, 2*M_PI);
    cairo_stroke(cr);
}

void draw_line(cairo_t *cr, float x0, float y0, float x1, float y1) {
    cairo_move_to(cr, x0, y0);
    cairo_line_to(cr, x1, y1);
    cairo_stroke(cr);
}

static int unpack_generic(lulog *dbg, const char *pattern, int **offsets, int *length, const char stop, int *padding) {

    LU_STATUS

    const char *p = pattern;
    int sign = 1;

    while (*p != stop) {
        if (*p == ',') sign = 1;
        else if (*p == '-') sign = -sign;
        else if (*p >= '0' && *p <= '9'){
            int offset = *p - '0';
            if (!(*offsets = realloc(*offsets, (1 + *length) * sizeof(**offsets)))) return LU_ERR_MEM;
            (*offsets)[*length] = sign * offset;
            (*length)++;
            sign = 1;
        } else {
            luerror(dbg, "Bad character '%c'", *p);
            status = LU_ERR;
            goto exit;
        }
        p++;
    }
    p++;   // drop stop character
    if (*p) {
        *padding = *p - '0';
    } else {
        *padding = 0;
    }

    LU_NO_CLEANUP
}

static int apply_padding(int **offsets, int *length, int padding) {

    LU_STATUS

    if (padding) {
        if (!(*offsets = realloc(*offsets, (padding + *length) * sizeof(**offsets)))) return LU_ERR_MEM;
        for (int i = 0; i < padding; ++i) (*offsets)[*length + i] = 0;
        *length = padding + *length;
    }

    LU_NO_CLEANUP
}

static int unpack_c(lulog *dbg, const char *pattern, int **offsets, int *length) {

    LU_STATUS
    int padding;

    LU_CHECK(unpack_generic(dbg, pattern, offsets, length, 'C', &padding))
    LU_CHECK(apply_padding(offsets, length, padding))

    LU_NO_CLEANUP
}

static int unpack_b(lulog *dbg, const char *pattern, int **offsets, int *length, int *padding) {

    LU_STATUS

    LU_CHECK(unpack_generic(dbg, pattern, offsets, length, 'B', padding))

    if (!(*offsets = realloc(*offsets, (2 * *length) * sizeof(**offsets)))) return LU_ERR_MEM;
    for (int i = 0; i < *length; ++i) (*offsets)[*length + i] = -(*offsets)[*length - i - 1];
    *length = 2 * *length;

    LU_CHECK(apply_padding(offsets, length, *padding))

    LU_NO_CLEANUP
}

static int unpack_a(lulog *dbg, const char *pattern, int **offsets, int *length, int *padding) {

    LU_STATUS

    LU_CHECK(unpack_generic(dbg, pattern, offsets, length, 'A', padding))

    if ((*offsets)[*length-1]) luwarn(dbg, "Central offset for group A is non-zero");
    if (!(*offsets = realloc(*offsets, (2 * *length - 1) * sizeof(**offsets)))) return LU_ERR_MEM;
    for (int i = 1; i < *length; ++i) (*offsets)[*length + i - 1] = -(*offsets)[*length - i - 1];
    *length = 2 * *length - 1;

    LU_CHECK(apply_padding(offsets, length, *padding))

    LU_NO_CLEANUP
}

int unpack(lulog *dbg, const char *pattern, int **offsets, int *length, char *type, int *padding) {

    LU_STATUS

    if (strchr(pattern, 'A')) {
        LU_CHECK(unpack_a(dbg, pattern, offsets, length, padding));
        if (type) *type = 'A';
    } else if (strchr(pattern, 'B')) {
        LU_CHECK(unpack_b(dbg, pattern, offsets, length, padding));
        if (type) *type = 'B';
    } else if (strchr(pattern, 'C')) {
        LU_CHECK(unpack_c(dbg, pattern, offsets, length));
        if (type) *type = 'C';
        *padding = 0;
    } else {
        luerror(dbg, "Did not find group type (A, B, C) in %s", pattern);
        status = LU_ERR_ARG;
    }

    LU_NO_CLEANUP
}

int dump_pattern(lulog *dbg, int *offsets, int length) {

    LU_STATUS
    char *buffer = NULL, *p;

    LU_ALLOC(dbg, buffer, length * 3 + 4)
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

int rim_size(lulog *dbg, int length, int *holes) {
    LU_STATUS
    switch(length) {
    case 1:
    case 2:
    case 4:
    case 8:
        *holes = 32;
        break;
    case 3:
    case 6:
    case 9:
        *holes = 36;
        break;
    case 7:
        *holes = 28;
        break;
    case 5:
    case 10:
        *holes = 20;
        break;
    default:
        luerror(dbg, "Cannot infer number of holes for a pattern of length %d", length);
        status = LU_ERR_ARG;
    }
    if (!status) luinfo(dbg, "Will use a rim with %d holes", *holes);
    LU_RETURN
}

int make_path(lulog *dbg, const char *pattern, char **path) {

    LU_STATUS
    const char *p1;
    char *p2;

    LU_ALLOC(dbg, *path, strlen(pattern) + 5)
    p1 = pattern; p2 = *path;
    while(*p1) {
        if (*p1 != ',') *(p2++) = *p1;
        p1++;
    }
    p2 += sprintf(p2, ".png");
    *p2 = '\0';

    luinfo(dbg, "File path: %s", *path);

LU_CLEANUP
    LU_RETURN
}


