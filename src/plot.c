
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <cairo/cairo.h>

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"


lulog *dbg = NULL;


int unpack_c(const char *pattern, int **offsets, int *length) {
    return LU_OK;
}

int unpack_b(const char *pattern, int **offsets, int *length) {

    LU_STATUS

    const char *p = pattern;
    int sign = 1;

    while (*p != 'B') {
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
    p++;   // drop B
    if (!(*offsets = realloc(*offsets, (2 * *length) * sizeof(**offsets)))) return LU_ERR_MEM;
    for (int i = 0; i < *length; ++i) (*offsets)[*length + i] = -(*offsets)[*length - i - 1];
    *length = 2 * *length;
    if (*p) {
        int padding = *p - '0';
        if (!(*offsets = realloc(*offsets, (padding + *length) * sizeof(**offsets)))) return LU_ERR_MEM;
        for (int i = 0; i < padding; ++i) (*offsets)[*length + i] = 0;
        *length = padding + *length;
    }

    LU_NO_CLEANUP
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
    for (int i = 1; i < *length; ++i) (*offsets)[*length + i - 1] = -(*offsets)[*length - i - 1];
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

int rim_size(int length, int *holes) {
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
    case 5:
    case 10:
        *holes = 20;
        break;
    default:
        luerror(dbg, "Cannot infer number of holes for a pattern of length %d", length);
        status = LU_ERR_ARG;
    }
    if (!status) luinfo(dbg, "Will plot a rim with %d holes", *holes);
    LU_RETURN
}

void draw_circle(cairo_t *cr, float r) {
    cairo_move_to(cr, r, 0);
    cairo_arc(cr, 0, 0, r, 0, 2*M_PI);
    cairo_stroke(cr);
}

void draw_spoke(cairo_t *cr, int hub, int offset, float r_hub, float r_rim, int holes) {
    float fudge = -M_PI / 2;  // rotate so red is in a nice place
    float t_hub = 2 * M_PI * hub / holes + fudge;
    float t_rim = 2 * M_PI * (hub + 2 * offset) / holes + fudge;
    cairo_move_to(cr, r_hub * cos(t_hub), r_hub * sin(t_hub));
    cairo_line_to(cr, r_rim * cos(t_rim), r_rim * sin(t_rim));
    cairo_stroke(cr);
}

void draw_pattern(cairo_t *cr, int *offsets, int length, float r_hub, float r_rim, int holes, int start, int direction) {
    for (int i = 0; i < length; ++i) {
        draw_spoke(cr, start + 2 * i * direction, offsets[i] * direction, r_hub, r_rim, holes);
    }
}

int draw(int *offsets, int length, int holes, const char *path) {

    LU_STATUS;
    cairo_surface_t *surface;
    cairo_t *cr;
    int nx = 200, ny = 200;
    float r_hub = 0.15, r_rim = 0.9, wheel_width = 0.03, wheel_grey = 0.5;
    float spoke_width = 0.015, red = 0.3, spoke_grey = 0.7;

    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, nx, ny);
    cr = cairo_create(surface);

    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_paint(cr);

    // centre at (0,0) with corner at (1,1)
    cairo_translate(cr, nx/2, ny/2);
    cairo_scale(cr, nx/2, ny/2);

    cairo_set_line_width(cr, wheel_width);
    cairo_set_source_rgb(cr, wheel_grey, wheel_grey, wheel_grey);
    draw_circle(cr, r_hub);
    draw_circle(cr, r_rim);

    cairo_set_line_width(cr, spoke_width);
    cairo_set_source_rgb(cr, spoke_grey, spoke_grey, spoke_grey);
    for (int i = 0; i < holes / (2 * length); ++i) {
        draw_pattern(cr, offsets, length, r_hub, r_rim, holes, 2 * i * length + 1, -1);
    }
    cairo_set_source_rgb(cr, 0, 0, 0);
    for (int i = 1; i < holes / (2 * length); ++i) {
        draw_pattern(cr, offsets, length, r_hub, r_rim, holes, 2 * i * length, 1);
    }
    cairo_set_source_rgb(cr, red, 0, 0);
    draw_pattern(cr, offsets, length, r_hub, r_rim, holes, 0, 1);

    cairo_surface_write_to_png(surface, path);

LU_CLEANUP
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    LU_RETURN
}

int make_path(const char *pattern, char **path) {

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

int plot(const char *pattern) {

    LU_STATUS;
    int *offsets = NULL, length = 0, holes = 0;
    char *path = NULL;

    luinfo(dbg, "Pattern '%s'", pattern);
    LU_CHECK(unpack(pattern, &offsets, &length));
    LU_CHECK(dump_pattern(offsets, length));
    LU_CHECK(rim_size(length, &holes));
    LU_CHECK(make_path(pattern, &path));
    LU_CHECK(draw(offsets, length, holes, path));

LU_CLEANUP
    free(path);
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
