
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "cairo/cairo.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_roots.h"

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"

#include "lib.h"

lulog *dbg = NULL;
gsl_rng *rng = NULL;

#define X 0
#define Y 1
#define G 9.8


typedef struct {
    double x;
    double y;
} xy;

typedef struct {
    double a11;
    double a12;
    double a21;
    double a22;
} m;

typedef struct {
    char type;
    int n_offsets;
    int *offset;
    int align;        // extra offset of "other side"
    int n_holes;
    double r_hub;
    double r_rim;     // uncompressed
    double l_chord;
    xy *rim;          // vary
    xy *hub;          // fixed
    int *hub_to_rim;  // spoke from hub[n] to rim[hub_to_rim[n]] has length[n]
    int *rim_to_hub;  // inverse of above
    double *l_spoke;  // unloaded, indexed by hub
    double tension;   // target tension for spokes
    double e_rim;
    double e_spoke;
} wheel;

typedef struct {
    int i_rim;
    double mass;
    xy g_norm;
    xy start;
    xy end;
} load;

typedef struct {
    int i_rim;
    m *to_xy;
    xy resolve;
    xy zero;
    wheel *wheel;
    load *load;
} data;

int make_wheel(int *offsets, int length, int holes, int padding, char type, wheel **wheel) {
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

int copy_wheel(wheel *w, wheel **c) {
    LU_STATUS
    LU_CHECK(make_wheel(w->offset, w->n_offsets, w->n_holes, w->align, w->type, c))
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

xy mul(m m, xy v) {
    xy c;
    c.x = m.a11 * v.x + m.a12 * v.y;
    c.y = m.a21 * v.x + m.a22 * v.y;
    return c;
}

double dot(xy a, xy b) {
    return a.x * b.x + a.y * b.y;
}

double cross(xy a, xy b) {
    return a.x * b.y - a.y * b.x;
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
    cairo_scale(*cr, nx/(2.2 * wheel->r_rim), ny/(2.2 * wheel->r_rim));

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

xy zoom(xy a, xy b, double r_scale, double t_scale) {
    xy d = sub(b, a);
    double r = sqrt(d.x * d.x + d.y * d.y), theta = atan2(d.y, d.x);
    r *= r_scale; theta *= t_scale;
    d.x = r * cos(theta); d.y = r * sin(theta);
    return add(a, d);
}

void draw_deform(cairo_t *cr, wheel *a, wheel *b, double r_scale, double t_scale) {
    for (int i = 0; i < a->n_holes; ++i) {
        int j = a->hub_to_rim[i];
        xy hub_a = a->hub[i], hub_b = b->hub[i];
        xy rim1_a = a->rim[j], rim1_b = b->rim[j];
        int k = (j + 1) % a->n_holes;
        xy rim2_a = a->rim[k], rim2_b = b->rim[k];
        draw_line_xy(cr, zoom(hub_a, hub_b, r_scale, t_scale), zoom(rim1_a, rim1_b, r_scale, t_scale));
        draw_line_xy(cr, zoom(rim1_a, rim1_b, r_scale, t_scale), zoom(rim2_a, rim2_b, r_scale, t_scale));
    }
}

void plot_deform(wheel *a, wheel *b, const char *path) {
    cairo_surface_t *surface = NULL;
    cairo_t *cr = NULL;
    open_plot(a, 500, 500, a->r_rim / 100, &cr, &surface);
    cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
    draw_wheel(cr, a);
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    draw_deform(cr, a, b, 1e2, 1e2);
    close_plot(cr, surface, path);
}

double angle_to_rim(wheel *wheel, int i_hub) {
    int i_rim = wheel->hub_to_rim[i_hub];
    xy hub = wheel->hub[i_hub];
    xy rim = wheel->rim[i_rim];
    // angle of spoke relative to horizontal
    double theta = atan2(rim.y - hub.y, rim.x - hub.x);
    // angle of radius to horizontal
    double psi = atan2(rim.y, rim.x);
    // angle of spoke relative to tangent, where 0 = along rim anticlock,
    // 90 = towards axle, 180 = along rim clock (but in radians...)
    return M_PI / 2 - psi + theta;
}

xy force(wheel *w, int i_rim, xy rim) {

    int i_before = (i_rim - 1 + w->n_holes) % w->n_holes;
    int i_after = (i_rim + 1 + w->n_holes) % w->n_holes;
    int i_hub = w->rim_to_hub[i_rim];

    xy spoke = sub(rim, w->hub[i_hub]);
    double l = length(spoke);
    double l0 = w->l_spoke[i_hub];
    xy f = {-w->e_spoke * (l - l0) * spoke.x / (l * l0), -w->e_spoke * (l - l0) * spoke.y / (l * l0)};

    xy chord = sub(rim, w->rim[i_before]);
    l = length(chord);
    l0 = w->l_chord;
    f.x -= w->e_rim * (l - l0) * chord.x / (l * l0);
    f.y -= w->e_rim * (l - l0) * chord.y / (l * l0);

    chord = sub(rim, w->rim[i_after]);
    l = length(chord);
    l0 = w->l_chord;
    f.x -= w->e_rim * (l - l0) * chord.x / (l * l0);
    f.y -= w->e_rim * (l - l0) * chord.y / (l * l0);

    return f;
}

double f_vertical(double v, void *params) {

    data *d = (data*)params;
    wheel *w = d->wheel;

    xy uv = {0, v};
    xy rim = add(d->zero, mul(*d->to_xy, uv));
//    ludebug(dbg, "Vertical %g: (%g,%g)", v, rim.x, rim.y);

    xy fxy = force(w, d->i_rim, rim);
    double f = dot(d->resolve, fxy);
//    ludebug(dbg, "Force (%g,%g) at %g: %g", fxy.x, fxy.y, v, f);
    return f;
}

double f_horizontal(double u, void *params) {

    data *d = (data*)params;
    wheel *w = d->wheel;

    xy uv = {u, 0};
    xy rim = add(d->zero, mul(*d->to_xy, uv));
//    ludebug(dbg, "Horizontal %g: (%g,%g)", u, rim.x, rim.y);

    xy fxy = force(w, d->i_rim, rim);
    double f = dot(d->resolve, fxy);
//    ludebug(dbg, "Force (%g,%g) at %g: %g", fxy.x, fxy.y, u, f);
    return f;
}

#define ITER_MAX 100
#define EPSABS 1e-3
#define EPSREL 1e-8

int relax_hole(wheel *w, int i_rim, double delta) {

    LU_STATUS
    data params;
    params.i_rim = i_rim;
    params.wheel = w;
    gsl_root_fsolver *s = NULL;
    LU_ASSERT(s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent), LU_ERR, dbg, "Cannot create solver")

    int i_before = (i_rim - 1 + w->n_holes) % w->n_holes;
    int i_after = (i_rim + 1 + w->n_holes) % w->n_holes;
    xy before = w->rim[i_before], after = w->rim[i_after];
    ludebug(dbg, "Before: (%g,%g); After: (%g,%g)", before.x, before.y, after.x, after.y);
    xy u = norm(sub(after, before));       // x axis
    xy v = {-u.y, u.x};                    // y axis, pointing away from hub
    ludebug(dbg, "X axis: (%g,%g); Y axis: (%g,%g)", u.x, u.y, v.x, v.y);
    m to_uv = {u.x, u.y, v.x, v.y};
    m to_xy = {v.y, -u.y, -v.x, u.x};      // inverse of to_uv
    params.to_xy = &to_xy;

    params.zero = w->rim[i_rim];
    params.resolve = u;

    gsl_function f;
    f.function = f_horizontal;
    f.params = &params;
    LU_ASSERT(!gsl_root_fsolver_set(s, &f, -0.5, 0.5), LU_ERR, dbg, "Cannot set solver")

    for (int iter = 0; iter < ITER_MAX; ++iter) {
        LU_ASSERT(!gsl_root_fsolver_iterate(s), LU_ERR, dbg, "Solver failed");
        double x_lo = gsl_root_fsolver_x_lower(s);
        double x_hi = gsl_root_fsolver_x_upper(s);
        double x = gsl_root_fsolver_root(s);
//        ludebug(dbg, "Tangential: %d; %g (%g - %g)", iter, x, x_lo, x_hi);
//        int gsl_status = gsl_root_test_interval(x_lo, x_hi, EPSABS, EPSREL);
        double residual = f.function(x, f.params);
        int gsl_status = gsl_root_test_residual(residual, EPSABS);
        if (gsl_status == GSL_SUCCESS) {ludebug(dbg, "Horizontal residual: %g", residual); break;}
        LU_ASSERT(gsl_status == GSL_CONTINUE, LU_ERR, dbg, "Test failed");
    }
    xy uv = {gsl_root_fsolver_root(s), 0};
    luinfo(dbg, "Tangential: %g", uv.x);
    params.zero = add(params.zero, mul(*params.to_xy, uv));

    xy q = sub(params.zero, before);
    double d = cross(u, q);                // v coord of hole above u axis
    ludebug(dbg, "Distance from (%g,%g) to horizontal: %g", params.zero.x, params.zero.y, d);
    LU_ASSERT(d > 0, LU_ERR, dbg, "Rim inflected");
    params.resolve = v;
    f.function = f_vertical;
    LU_ASSERT(!gsl_root_fsolver_set(s, &f, -d, 0.5), LU_ERR, dbg, "Cannot set solver")

    for (int iter = 0; iter < ITER_MAX; ++iter) {
        LU_ASSERT(!gsl_root_fsolver_iterate(s), LU_ERR, dbg, "Solver failed");
        double x_lo = gsl_root_fsolver_x_lower(s);
        double x_hi = gsl_root_fsolver_x_upper(s);
        double x = gsl_root_fsolver_root(s);
//        ludebug(dbg, "Radial: %d; %g (%g - %g)", iter, x, x_lo, x_hi);
//        int gsl_status = gsl_root_test_interval (x_lo, x_hi, EPSABS, EPSREL);
        double residual = f.function(x, f.params);
        int gsl_status = gsl_root_test_residual(residual, EPSABS);
        if (gsl_status == GSL_SUCCESS) {ludebug(dbg, "Vertical residual: %g", residual); break;}
        LU_ASSERT(gsl_status == GSL_CONTINUE, LU_ERR, dbg, "Test failed");
    }
    double root = gsl_root_fsolver_root(s);
    uv = (xy){0, gsl_root_fsolver_root(s)};
    luinfo(dbg, "Radial: %g", uv.y);
    params.zero = add(params.zero, mul(*params.to_xy, uv));

    xy full = sub(params.zero, w->rim[params.i_rim]);
    xy frac = scalar_mult(delta, full);
    ludebug(dbg, "Shift: (%g,%g) (Delta: %g)", frac.x, frac.y, delta);
    w->rim[params.i_rim] = add(w->rim[params.i_rim], frac);

LU_CLEANUP
    if (s) gsl_root_fsolver_free(s);
    LU_RETURN
}

float residual(wheel *wheel) {
    double residual = 0;
    for (int i = 0; i < wheel->n_holes; ++i) {
        xy f = force(wheel, i, wheel->rim[i]);
        double delta = length(f);
        ludebug(dbg, "Force at %d: %g", i, delta);
        residual += delta;
    }
    return residual;
}

int relax(wheel *wheel, load *load) {

    LU_STATUS
    double force = 0;
    int iter = 0;

    do {
        force = 0; iter++;
        int start = gsl_rng_uniform_int(rng, wheel->n_holes);
        for (int i = 0; i < wheel->n_holes; ++i) {
            int i_rim = (i + start) % wheel->n_holes;
            ludebug(dbg, "Relaxing hole %d at (%g,%g)", i_rim, wheel->rim[i_rim].x, wheel->rim[i_rim].y);
            LU_CHECK(relax_hole(wheel, i_rim, 1 - pow(0.9, iter)))
        }
        force = residual(wheel);
        ludebug(dbg, "Iteration: %d; Residual: %g", iter, force);
    } while (force > 10);

LU_CLEANUP
    LU_RETURN
}

#define TARGET_WOBBLE 1e-6
#define DAMP_TRUE 0.5
#define DAMP_TENSION 0.5

int true(wheel *w) {

    LU_STATUS

    LU_CHECK(relax(w, NULL))

    while(1) {

        double r_target = 0;
        for (int i = 0; i < w->n_holes; ++i) r_target += length(w->rim[i]);
        r_target = r_target / w->n_holes;
        double wobble = 0, shift = 0;
        for (int i = 0; i < w->n_holes; ++i) wobble += fabs(length(w->rim[i]) - r_target);
        wobble /= w->n_holes;
        luinfo(dbg, "Average radial wobble %gmm", wobble);
        if (wobble <= TARGET_WOBBLE) break;
        for (int i = 0; i < w->n_holes; ++i) {
            double l = length(w->rim[i]);
            if (l > r_target) {
                double excess = DAMP_TRUE * (l - r_target);
                w->l_spoke[w->rim_to_hub[i]] -= excess;
                shift += excess;
            }
        }
        luinfo(dbg, "Total correction %gmm", shift);
        LU_CHECK(relax(w, NULL));

        double tension = 0;
        for (int i = 0; i < w->n_holes; ++i) {
            int j = w->rim_to_hub[i];
            double l = length(sub(w->rim[i], w->hub[j]));
            double l0 = w->l_spoke[j];
            tension += w->e_spoke * (l - l0) / l0;
        }
        tension /= w->n_holes;
        // df/dl = e/l0 so dl = df l0 / e
        double error = tension - w->tension;
        luinfo(dbg, "Tension excess is %g", error);
        for (int i = 0; i < w->n_holes; ++i) {
            // if tension is too large, correction is positive and spoke is extended
            w->l_spoke[i] += DAMP_TENSION * w->l_spoke[i] * error / w->e_spoke;
        }
        LU_CHECK(relax(w, NULL));
    }

    luinfo(dbg, "True!");

LU_CLEANUP
    LU_RETURN
}

int lace(wheel *wheel) {

    LU_STATUS
    double *tension = NULL;

    LU_ALLOC(dbg, tension, wheel->n_holes);

    for (int i_hub = 0; i_hub < wheel->n_holes; i_hub += 2) {
        int offset = wheel->offset[(i_hub / 2) % wheel->n_offsets];
        int i_rim = (i_hub + 2 * offset + wheel->n_holes) % wheel->n_holes;
        wheel->hub_to_rim[i_hub] = i_rim;
        wheel->rim_to_hub[i_rim] = i_hub;
        wheel->hub[i_hub] = xy_on_circle(wheel->r_hub, i_hub, wheel->n_holes);
        wheel->rim[i_rim] = xy_on_circle(wheel->r_rim, i_rim, wheel->n_holes);
    }
    for (int i = 0; i < wheel->n_holes; i += 2) {
        int offset = wheel->offset[(i / 2) % wheel->n_offsets];
        int i_hub = (i + 2 * wheel->align + 1 + wheel->n_holes) % wheel->n_holes;
        int i_rim = (i + 2 * (offset + wheel->align) + 1 + wheel->n_holes) % wheel->n_holes;
        wheel->hub_to_rim[i_hub] = i_rim;
        wheel->rim_to_hub[i_rim] = i_hub;
        wheel->hub[i_hub] = xy_on_circle(wheel->r_hub, i_hub, wheel->n_holes);
        wheel->rim[i_rim] = xy_on_circle(wheel->r_rim, i_rim, wheel->n_holes);
    }

    // set lengths so that all have desired radial tension
    for (int i_hub = 0; i_hub < wheel->n_holes; ++i_hub) {
        int i_rim = wheel->hub_to_rim[i_hub];
        double l = length(sub(wheel->hub[i_hub], wheel->rim[i_rim]));
        double theta = angle_to_rim(wheel, i_hub);
        double strain = wheel->tension / (sin(theta) * wheel->e_spoke);
        wheel->l_spoke[i_hub] = l * (1 - strain);
    }

LU_CLEANUP
    free(tension);
    LU_RETURN
}

int deform(wheel *wheel) {

    LU_STATUS
    load *l = NULL;

    LU_ALLOC(dbg, l, 1)
    l->g_norm.x = -1;
    l->g_norm.y = 0;
    l->i_rim = 4;
    l->start = wheel->rim[l->i_rim];

    for (int i = -3; i < 2; ++i) {
        l->mass = pow(10, i+1);
        ludebug(dbg,"Mass %gkg", l->mass);
        LU_CHECK(relax(wheel, l))
    }

LU_CLEANUP
    free(l);
    LU_RETURN
}

int stress(const char *pattern) {

    LU_STATUS;
    int *offsets = NULL, length = 0, holes = 0, padding;
    char type, *path = NULL;
    wheel *wheel = NULL, *original = NULL;

    luinfo(dbg, "Pattern '%s'", pattern);
    LU_CHECK(unpack(dbg, pattern, &offsets, &length, &type, &padding))
    LU_ASSERT(strchr("AB", type), LU_ERR, dbg, "Only symmetric types supported")
    LU_CHECK(dump_pattern(dbg, offsets, length))
    LU_CHECK(rim_size(dbg, length, &holes))
    LU_CHECK(make_wheel(offsets, length, holes, padding, type, &wheel))
    LU_CHECK(lace(wheel))
    LU_CHECK(true(wheel))
    LU_CHECK(copy_wheel(wheel, &original))
//    LU_CHECK(deform(wheel))
    LU_CHECK(make_path(dbg, pattern, &path))
//    plot_deform(original, wheel, path);
    plot_wheel(original, path);

LU_CLEANUP
    free_wheel(original);
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

int main(int argc, char** argv) {

    LU_STATUS

    lulog_mkstderr(&dbg, lulog_level_debug);
    if (argc != 2 || !strcmp("-h", argv[1])) {
        usage(argv[0]);
    } else {
        LU_ASSERT(rng = gsl_rng_alloc(gsl_rng_mt19937), LU_ERR, dbg, "Could not create PRNG")
        LU_CHECK(stress(argv[1]));
    }

LU_CLEANUP
    if (dbg) status = dbg->free(&dbg, status);
    if (rng) gsl_rng_free(rng);
    return status;
}
