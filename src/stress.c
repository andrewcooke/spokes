
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>

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
volatile int sig_exit = 0;

#define X 0
#define Y 1
#define G 9.8

#define ZOOM 100

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
} load;

typedef struct {
    int i_rim;
    xy descent;
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
    draw_deform(cr, a, b, ZOOM, ZOOM);
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

void describe(wheel *w) {
    for (int i = 0; i < w->n_holes; ++i) {
        double ls = length(sub(w->rim[i], w->hub[w->rim_to_hub[i]]));
        double l0s = w->l_spoke[w->rim_to_hub[i]];
        double xs = (ls - l0s) / l0s;
        double fs = w->e_spoke * xs;
        int j = (i - 1 + w->n_holes) % w->n_holes;
        double lr = length(sub(w->rim[i], w->rim[j]));
        double l0r = w->l_chord;
        double xr = (lr - l0r) / l0r;
        double fr = w->e_rim * xr;
        luinfo(dbg, "%d: Spoke %g%% %gN; Rim %g%% %gN", i, 100 * xs, fs, 100 * xr, fr);
    }
}

xy force(wheel *w, int i_rim, xy rim, load *load) {

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

    if (load && i_rim == load->i_rim) {
        xy extra = scalar_mult(load->mass * G, load->g_norm);
        f = add(extra, f);
    }

    return f;
}

double f_descent(double x, void *params) {

    data *d = (data*)params;
    wheel *w = d->wheel;

    xy rim = add(d->zero, scalar_mult(x, d->descent));
//    ludebug(dbg, "Descent %g: (%g,%g)", x, rim.x, rim.y);

    xy fxy = force(w, d->i_rim, rim, d->load);
    double f = dot(d->descent, fxy);
//    ludebug(dbg, "Force (%g,%g) at %g: %g", fxy.x, fxy.y, x, f);
    return f;
}

#define ITER_MAX 100
#define EPSABS 1e-3
#define EPSREL 1e-8

int relax_hole(wheel *w, load *load, int i_rim, double delta, int verbose) {

    LU_STATUS

    data params;
    params.load = load;
    params.i_rim = i_rim;
    params.wheel = w;
    params.zero = w->rim[i_rim];
    // steepest descent direction is the direction of the force
    params.descent = norm(force(w, i_rim, params.zero, load));
    if (verbose) {
        ludebug(dbg, "Relaxing at %d (%g,%g) along (%g,%g)",
                i_rim, params.zero.x, params.zero.y, params.descent.x, params.descent.y);
    }

    gsl_root_fsolver *s = NULL;
    LU_ASSERT(s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent), LU_ERR, dbg, "Cannot create solver")

    // calculate limit for range
    int i_before = (i_rim - 1 + w->n_holes) % w->n_holes;
    int i_after = (i_rim + 1 + w->n_holes) % w->n_holes;
    xy before = w->rim[i_before], after = w->rim[i_after];
    if (verbose) ludebug(dbg, "Before: (%g,%g); After: (%g,%g)", before.x, before.y, after.x, after.y);
    xy u = norm(sub(after, before));         // x axis
    xy v = {-u.y, u.x}; // y axis, pointing away from hub
    xy q = sub(params.zero, before);
    // distance above line in uv frame divided by sin(angle)
    double du = cross(u, q);
    double d = -du / dot(v, params.descent); // distance along descent line
    if (verbose) ludebug(dbg, "Limit %g in uv; %g along descent", du, d);
    double lo = -fabs(d), hi = fabs(d);
    lo = lo < -10 ? -10 : lo; hi = hi > 10 ? 10 : hi;
    if (verbose) ludebug(dbg, "Range: %g - %g along (%g,%g)", lo, hi, params.descent.x, params.descent.y);

    gsl_function f;
    f.function = f_descent;
    f.params = &params;
    LU_ASSERT(!gsl_root_fsolver_set(s, &f, lo, hi), LU_ERR, dbg, "Cannot set solver")

    for (int iter = 0; iter < ITER_MAX; ++iter) {
        LU_ASSERT(!gsl_root_fsolver_iterate(s), LU_ERR, dbg, "Solver failed");
        double x_lo = gsl_root_fsolver_x_lower(s);
        double x_hi = gsl_root_fsolver_x_upper(s);
        double x = gsl_root_fsolver_root(s);
        if (verbose) ludebug(dbg, "Descent: %d; %g (%g - %g)", iter, x, x_lo, x_hi);
//        int gsl_status = gsl_root_test_interval(x_lo, x_hi, EPSABS, EPSREL);
        double residual = f.function(x, f.params);
        int gsl_status = gsl_root_test_residual(residual, EPSABS);
        if (gsl_status == GSL_SUCCESS) break;
        LU_ASSERT(gsl_status == GSL_CONTINUE, LU_ERR, dbg, "Test failed");
    }

    double root = gsl_root_fsolver_root(s);
    xy end = add(params.zero, scalar_mult(root, params.descent));
    xy f_end = force(w, i_rim, end, load);
    if (verbose) ludebug(dbg, "End point: (%g,%g); Force: (%g,%g)", end.x, end.y, f_end.x, f_end.y);

    xy full = sub(end, params.zero);
    xy frac = scalar_mult(delta, full);
    w->rim[params.i_rim] = add(w->rim[params.i_rim], frac);
    if (verbose) {
        luinfo(dbg, "Shift: (%g,%g) (Delta: %g) to (%g,%g)",
                frac.x, frac.y, delta, w->rim[params.i_rim].x, w->rim[params.i_rim].y);
    }

LU_CLEANUP
    if (s) gsl_root_fsolver_free(s);
    LU_RETURN
}

float residual(wheel *wheel, load *load) {
    double residual = 0;
    for (int i = 0; i < wheel->n_holes; ++i) {
        xy f = force(wheel, i, wheel->rim[i], load);
        double delta = length(f);
//        ludebug(dbg, "Force at %d: %g", i, delta);
        residual += delta;
    }
    return residual;
}

// damp should be 0.9 for an initial step of 0.1, 0.99 for 0.01, etc
int relax(wheel *wheel, load *load, double step, double damp, int verbose) {

    LU_STATUS
    double force = 0;
    int iter = 0;

    do {
        force = 0; iter++;
        int start = gsl_rng_uniform_int(rng, wheel->n_holes);
        int sign = (2 * gsl_rng_uniform_int(rng, 2)) - 1;
        if (iter == 1 && load) start = load->i_rim;
        for (int i = 0; i < wheel->n_holes; ++i) {
            int i_rim = (sign * i + start + 2 * wheel->n_holes) % wheel->n_holes;
//            ludebug(dbg, "Relaxing hole %d at (%g,%g)", i_rim, wheel->rim[i_rim].x, wheel->rim[i_rim].y);
            LU_CHECK(relax_hole(wheel, load, i_rim, step * (1 - pow(damp, iter)), verbose))
        }
        force = residual(wheel, load);
        if (sig_exit || force <= 10 || iter < 100 || iter % 1000 == 0) luinfo(dbg, "Iteration: %d; Residual: %g", iter, force);
    } while (sig_exit || force > 10);
    if (sig_exit) luwarn(dbg, "SIGINT aborted relax");

LU_CLEANUP
    LU_RETURN
}

#define TARGET_WOBBLE 1e-5
#define DAMP_TRUE 0.5
#define DAMP_TENSION 0.5

int true(wheel *w) {

    LU_STATUS

    LU_CHECK(relax(w, NULL, 0.5, 0.9, 0))

    while(!sig_exit) {

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
        LU_CHECK(relax(w, NULL, 0.5, 0, 0));

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
        LU_CHECK(relax(w, NULL, 0.5, 0, 0));
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
    l->g_norm.x = 0;
    l->g_norm.y = -1;
    l->i_rim = 0;

//    for (int i = 0; i < 10; ++i) {
//        l->mass = 10 * (i+1);
//        luinfo(dbg,"Mass %gkg", l->mass);
//        LU_CHECK(relax(wheel, l))
//    }
    l->mass = 20;
    LU_CHECK(relax(wheel, l, 0.1, 0.99, 0))

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
    describe(wheel);
    LU_CHECK(copy_wheel(wheel, &original))
    LU_CHECK(deform(wheel))
    describe(wheel);
    LU_CHECK(make_path(dbg, pattern, &path))
    plot_deform(original, wheel, path);
//    plot_wheel(original, path);

LU_CLEANUP
    free_wheel(original);
    free_wheel(wheel);
    free(path);
    free(offsets);
    LU_RETURN
}

void new_handler(int sig) {
//    luwarn(dbg, "Handler called with %d", sig);
    sig_exit = 1;
}

int set_handler() {
    LU_STATUS
    struct sigaction action = {0};
    action.sa_handler = new_handler;
    LU_ASSERT(!sigaction(SIGINT, &action, NULL), LU_ERR, dbg, "Could not set handler")
    luinfo(dbg, "Handler set");
    LU_NO_CLEANUP
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
        LU_CHECK(set_handler())
        while (!sig_exit) sleep(1);
//        gsl_set_error_handler_off();
        LU_ASSERT(rng = gsl_rng_alloc(gsl_rng_mt19937), LU_ERR, dbg, "Could not create PRNG")
        LU_CHECK(stress(argv[1]));
    }

LU_CLEANUP
    if (dbg) status = dbg->free(&dbg, status);
    if (rng) gsl_rng_free(rng);
    return status;
}
