
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <cairo/cairo.h>

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"

#include "lib.h"

lulog *dbg = NULL;


typedef struct {
    double x;
    double y;
} xy;

typedef struct {
    char type;
    int n_offsets;
    int *offset;
    int align;        // extra offset of "other side"
    int n_holes;
    xy *nipple;
    double r_hub;
    double windup;    // angle of rotation of hub
    double r_rim;     // uncompressed
    double *l_spoke;  // slack
    double e_rim;
    double e_spoke;
} wheel;


// hole 0 at 12 o'clock
xy hole_xy(double r, int hole, int n_holes, double delta) {
    xy xy;
    double theta = delta + 2 * M_PI * hole / (float)n_holes;
    xy.x = r * sin(theta);
    xy.y = r * cos(theta);
    return xy;
}

void plot_wheel(wheel *wheel, const char *path) {

    int nx = 500, ny = 500;
    double line_width = wheel->r_rim / 100;

    cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, ny, ny);
    cairo_t *cr = cairo_create(surface);

    cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
    cairo_paint(cr);

    cairo_translate(cr, nx/2, ny/2);
    cairo_scale(cr, nx/(2.2 * wheel->r_rim), ny/(2.2 * wheel->r_rim));

    cairo_set_line_width(cr, line_width);
    cairo_set_source_rgb(cr, 0, 0, 0);
    draw_circle(cr, wheel->r_hub);
    draw_circle(cr, wheel->r_rim);

    for (int i = 0; i < wheel->n_holes; ++i) {
        xy hub = hole_xy(wheel->r_hub, i, wheel->n_holes, wheel->windup);
        xy rim = wheel->nipple[i];
        draw_line(cr, hub.x, hub.y, rim.x, rim.y);
    }


    cairo_surface_write_to_png(surface, path);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);

}

void widths_for_tension(wheel *wheel) {

}

void tension_for_circle(wheel *wheel) {

}

void lace(wheel *wheel) {
    for (int i = 0; i < wheel->n_holes; i += 2) {
        int offset = wheel->offset[(i / 2) % wheel->n_offsets];
        wheel->nipple[i] = hole_xy(wheel->r_rim, i + 2 * offset, wheel->n_holes, 0);
    }
    for (int i = 0; i < wheel->n_holes; i += 2) {
        int offset = wheel->offset[(i / 2) % wheel->n_offsets];
        wheel->nipple[(i + 2 * wheel->align + 1) % wheel->n_holes] = hole_xy(wheel->r_rim, i + 2 * (offset + wheel->align) + 1, wheel->n_holes, 0);
    }
    for (int repeat = 0; repeat < 2; ++repeat) {
        tension_for_circle(wheel);
        widths_for_tension(wheel);
    }
}

int make_wheel(int *offsets, int length, int holes, int padding, char type, wheel **wheel) {
    LU_STATUS
    LU_ALLOC(dbg, *wheel, 1);
    (*wheel)->type = type;
    (*wheel)->n_offsets = length;
    LU_ALLOC(dbg, (*wheel)->offset, length);
    for (int i = 0; i < length; ++i) (*wheel)->offset[i] = offsets[i];
    (*wheel)->align = padding;
    (*wheel)->n_holes = holes;
    LU_ALLOC(dbg, (*wheel)->nipple, holes);
    (*wheel)->r_hub = 25;
    (*wheel)->windup = 0;
    (*wheel)->r_rim = 280;
    LU_ALLOC(dbg, (*wheel)->l_spoke, holes);
    (*wheel)->e_spoke = 1;
    (*wheel)->e_rim = 100 * (*wheel)->e_spoke;
    LU_NO_CLEANUP
}

void free_wheel(wheel *wheel) {
    if (wheel) {
        free(wheel->offset);
        free(wheel->nipple);
        free(wheel->l_spoke);
        free(wheel);
    }
}

int stress(const char *pattern) {

    LU_STATUS;
    int *offsets = NULL, length = 0, holes = 0, padding;
    char type, *path = NULL;
    wheel *wheel = NULL;

    luinfo(dbg, "Pattern '%s'", pattern);
    LU_CHECK(unpack(dbg, pattern, &offsets, &length, &type, &padding));
    LU_ASSERT(strchr("AB", type), LU_ERR, dbg, "Only symmetric types supported")
    LU_CHECK(dump_pattern(dbg, offsets, length));
    LU_CHECK(rim_size(dbg, length, &holes));
    LU_CHECK(make_wheel(offsets, length, holes, padding, type, &wheel));
    lace(wheel);
    LU_CHECK(make_path(dbg, pattern, &path));
    plot_wheel(wheel, path);

LU_CLEANUP
    free_wheel(wheel);
    free(path);
    free(offsets);
    LU_RETURN
}

void usage(const char *progname) {
    luinfo(dbg, "Plot the stresses for a given spoke pattern");
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
        LU_CHECK(stress(argv[1]));
    }

LU_CLEANUP
    if (dbg) status = dbg->free(&dbg, status);
    return status;
}
