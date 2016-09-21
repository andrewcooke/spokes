
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

int draw(int *offsets, int length, int holes, int nx, int ny, int align, const char *path) {

    LU_STATUS;
    float r_hub = 0.085, r_rim = 0.9, wheel_width = 0.03, wheel_grey = 0.5;
    float spoke_width = 0.015, red = 0.5, spoke_grey = 0.7;

    cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, nx, ny);
    cairo_t *cr = cairo_create(surface);

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
        draw_pattern(cr, offsets, length, r_hub, r_rim, holes, 2 * i * length - 1 - 2 * align, -1);
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

int plot_size(char type, int *nx, int *ny) {
    switch(type) {
    case 'A':
    case 'B':
        *nx = *ny = 200;
        return LU_OK;
    case 'C':
        *nx = *ny = 100;
        return LU_OK;
    default:
        luerror(dbg, "Unexpected type %c", type);
        return LU_ERR;
    }
}

int plot(const char *pattern) {

    LU_STATUS;
    int *offsets = NULL, length = 0, holes = 0, nx = 0, ny = 0, padding;
    char *path = NULL, type;

    luinfo(dbg, "Pattern '%s'", pattern);
    LU_CHECK(unpack(dbg, pattern, &offsets, &length, &type, &padding));
    LU_CHECK(dump_pattern(dbg, offsets, length));
    LU_CHECK(rim_size(dbg, length, &holes));
    LU_CHECK(make_path(dbg, pattern, &path));
    LU_CHECK(plot_size(type, &nx, &ny));
    LU_CHECK(draw(offsets, length, holes, nx, ny, padding, path));

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
