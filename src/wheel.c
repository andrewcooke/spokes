
#include <math.h>

#include "cairo/cairo.h"

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/strings.h"
#include "lu/dynamic_memory.h"

#include "wheel.h"
#include "lib.h"


int make_wheel(lulog *dbg, int *offsets, int length, int holes, int padding, char type, wheel **wheel) {
    LU_STATUS
    LU_ALLOC(dbg, *wheel, 1);
    (*wheel)->type = type;
    (*wheel)->n_offsets = length;
    LU_ALLOC(dbg, (*wheel)->offset, length);
    for (int i = 0; i < length; ++i) (*wheel)->offset[i] = offsets[i];
    (*wheel)->align = padding;
    (*wheel)->n_holes = holes;
    LU_ALLOC(dbg, (*wheel)->rim, holes);
    LU_ALLOC(dbg, (*wheel)->hub, holes);
    LU_ALLOC(dbg, (*wheel)->hub_to_rim, holes);
    LU_ALLOC(dbg, (*wheel)->rim_to_hub, holes);
    (*wheel)->r_hub = 25;  // mm
    (*wheel)->r_rim = 280;  // mm
    (*wheel)->l_chord = 2 * (*wheel)->r_rim * sin(M_PI / (*wheel)->n_holes);
    LU_ALLOC(dbg, (*wheel)->l_spoke, holes);
    (*wheel)->tension = 1000;  // 100 kgf = 1000 N
    (*wheel)->e_spoke = 200000 * 3;  // 200000 N/mm^2 for steel approx; 2mm diameter spoke
    (*wheel)->e_rim = 100 * (*wheel)->e_spoke;
//    (*wheel)->e_rim = 10 * (*wheel)->e_spoke;
    LU_NO_CLEANUP
}

void free_wheel(wheel *wheel) {
    if (wheel) {
        free(wheel->offset);
        free(wheel->rim);
        free(wheel->hub);
        free(wheel->hub_to_rim);
        free(wheel->rim_to_hub);
        free(wheel->l_spoke);
        free(wheel);
    }
}

int copy_wheel(lulog *dbg, wheel *w, wheel **c) {
    LU_STATUS
    LU_CHECK(make_wheel(dbg, w->offset, w->n_offsets, w->n_holes, w->align, w->type, c))
    for (int i = 0; i < w->n_holes; ++i) {
        (*c)->hub[i] = w->hub[i];
        (*c)->rim[i] = w->rim[i];
        (*c)->hub_to_rim[i] = w->hub_to_rim[i];
        (*c)->rim_to_hub[i] = w->rim_to_hub[i];
    }
    LU_NO_CLEANUP
}

xy add(xy a, xy b) {
    xy c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return c;
}

xy sub(xy a, xy b) {
    xy c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    return c;
}

double dot(xy a, xy b) {
    return a.x * b.x + a.y * b.y;
}

double length(xy xy) {
    return sqrt(xy.x * xy.x + xy.y * xy.y);
}

xy scalar_mult(double k, xy a) {
    xy c;
    c.x = k * a.x;
    c.y = k * a.y;
    return c;
}

xy norm(xy a) {
    return scalar_mult(1.0 / length(a), a);
}


// hole 0 at 12 o'clock
xy xy_on_circle(double r, int hole, int n_holes) {
    xy xy;
    double theta = 2 * M_PI * hole / (float)n_holes;
    xy.x = r * sin(theta);
    xy.y = r * cos(theta);
    return xy;
}

void draw_line_xy(cairo_t *cr, xy a, xy b) {
    draw_line(cr, a.x, a.y, b.x, b.y);
}

void open_plot(wheel * wheel, int nx, int ny, double line_width, cairo_t **cr, cairo_surface_t **surface) {

    *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, ny, ny);
    *cr = cairo_create(*surface);

    cairo_set_source_rgb(*cr, 1.0, 1.0, 1.0);
    cairo_paint(*cr);

    cairo_translate(*cr, nx/2, ny/2);
    cairo_scale(*cr, nx/(2.2 * wheel->r_rim), -ny/(2.2 * wheel->r_rim));

    cairo_set_line_width(*cr, line_width);
    cairo_set_source_rgb(*cr, 0, 0, 0);
}

void close_plot(cairo_t *cr, cairo_surface_t *surface, const char *path) {
    cairo_surface_write_to_png(surface, path);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}

void draw_wheel(cairo_t *cr, wheel *wheel) {
    for (int hub = 0; hub < wheel->n_holes; ++hub) {
        int rim = wheel->hub_to_rim[hub];
        draw_line_xy(cr, wheel->hub[hub], wheel->rim[rim]);
        draw_line_xy(cr, wheel->hub[hub], wheel->hub[(hub+1) % wheel->n_holes]);
        draw_line_xy(cr, wheel->rim[rim], wheel->rim[(rim+1) % wheel->n_holes]);
    }
}

void plot_wheel(wheel *wheel, const char *path) {
    cairo_surface_t *surface = NULL;
    cairo_t *cr = NULL;
    open_plot(wheel, 500, 500, wheel->r_rim / 100, &cr, &surface);
    draw_wheel(cr, wheel);
    close_plot(cr, surface, path);
}

xy zoom(xy a, xy b, double scale) {
    xy d = sub(b, a);
    double r = sqrt(d.x * d.x + d.y * d.y), theta = atan2(d.y, d.x);
    r *= scale;
    d.x = r * cos(theta); d.y = r * sin(theta);
    return add(a, d);
}

void draw_deform(cairo_t *cr, wheel *a, wheel *b, double scale) {
    for (int i = 0; i < a->n_holes; ++i) {
        int j = a->hub_to_rim[i];
        xy hub_a = a->hub[i], hub_b = b->hub[i];
        xy rim1_a = a->rim[j], rim1_b = b->rim[j];
        int k = (j + 1) % a->n_holes;
        xy rim2_a = a->rim[k], rim2_b = b->rim[k];
        draw_line_xy(cr, zoom(hub_a, hub_b, scale), zoom(rim1_a, rim1_b, scale));
        draw_line_xy(cr, zoom(rim1_a, rim1_b, scale), zoom(rim2_a, rim2_b, scale));
    }
}

void plot_deform(wheel *original, wheel *deformed, load *l, const char *path, double scale) {
    cairo_surface_t *surface = NULL;
    cairo_t *cr = NULL;
    open_plot(original, 500, 500, original->r_rim / 100, &cr, &surface);
    cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
    draw_wheel(cr, original);
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    draw_deform(cr, original, deformed, scale);
    cairo_set_source_rgb(cr, 0.5, 0.0, 0.0);
    xy p = zoom(original->rim[l->i_rim], deformed->rim[l->i_rim], scale);
    draw_line_xy(cr, p, add(p, scalar_mult(100, l->g_norm)));
    close_plot(cr, surface, path);
}

int plot_multi_deform(lulog *dbg, wheel *original, wheel *deformed, load *l, const char *pattern) {
    LU_STATUS
    lustr path = {0};
    for (int i = 0; i < 6; ++i) {
        LU_CHECK(lustr_printf(dbg, &path, "%s-%d.png", pattern, i))
        double scale = pow(10, i);
        plot_deform(original, deformed, l, path.c, scale);
        luinfo(dbg, "Scale %g plot: %s", scale, path.c);
        LU_CHECK(lustr_clear(dbg, &path))
    }
LU_CLEANUP
    status = lustr_free(&path, status);
    LU_RETURN
}

