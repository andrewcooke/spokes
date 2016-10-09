
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "cairo/cairo.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_vector.h"

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"

#include "lib.h"

lulog *dbg = NULL;

#define X 0
#define Y 1
#define G 9.8


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
    double r_hub;
    double r_rim;     // uncompressed
    double l_chord;
    xy *rim;          // vary
    xy *hub;          // fixed
    int *hub_to_rim;  // spoke from hub[n] to rim[hub_to_rim[n]] has length[n]
    int *rim_to_hub;  // inverse of above
    double *l_spoke;  // unloaded
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
//    (*wheel)->e_rim = 100 * (*wheel)->e_spoke;
    (*wheel)->e_rim = 10 * (*wheel)->e_spoke;
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
    for (int i = 0; i < w->n_holes; ++i) (*c)->rim[i] = w->rim[i];
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

typedef struct {
    wheel *wheel;
    double *spoke_strain;
    xy *spoke;
    double *chord_strain;
    xy *chord;
    load *load;
} data;

void calculate_data(const gsl_vector *rim, data *d) {

    wheel *w = d->wheel;

    memset(d->spoke_strain, 0, w->n_holes * sizeof(*d->spoke_strain));
    memset(d->spoke, 0, w->n_holes * sizeof(*d->spoke));
    memset(d->chord_strain, 0, w->n_holes * sizeof(*d->chord_strain));
    memset(d->chord, 0, w->n_holes * sizeof(*d->chord));

    for (int i = 0; i < w->n_holes; ++i) {
        xy *hub = &w->hub[w->rim_to_hub[i]];
        xy *spoke = &d->spoke[i];
        spoke->x = gsl_vector_get(rim, 2*i+X) - hub->x;
        spoke->y = gsl_vector_get(rim, 2*i+Y) - hub->y;
        double l = w->l_spoke[w->rim_to_hub[i]];
        d->spoke_strain[i] = (length(*spoke) - l) / l;;
//        ludebug(dbg, "Spoke %d extended by %g%%", i, 100 * d->spoke_extension[i]);
    }

    for (int after = 0; after < w ->n_holes; ++after) {
        int before = (after - 1 + w->n_holes) % w->n_holes;
        xy *chord = &d->chord[after];
        chord->x = gsl_vector_get(rim, 2*after+X) - gsl_vector_get(rim, 2*before+X);
        chord->y = gsl_vector_get(rim, 2*after+Y) - gsl_vector_get(rim, 2*before+Y);
        d->chord_strain[after] = (length(*chord) - w->l_chord) / w->l_chord;
//        ludebug(dbg, "Chord %d compressed by %g%%", after, 100 * d->chord_compression[after]);
    }

    if (d->load) {
        d->load->end.x = gsl_vector_get(rim, 2*d->load->i_rim+X);
        d->load->end.y = gsl_vector_get(rim, 2*d->load->i_rim+Y);
    }
}

void calculate_energy(data *d, double *energy) {

    wheel *w = d->wheel;
    load *l = d->load;
    *energy = 0;

    for (int i = 0; i < w->n_holes; ++i) {
        double strain = d->spoke_strain[i];
        *energy += w->e_spoke * strain * strain / 2;
    }
    for (int i = 0; i < w ->n_holes; ++i) {
        double strain = d->chord_strain[i];
        *energy += w->e_rim * strain * strain / 2;
    }
    if (l) {
        xy disp = sub(l->end, l->start);
        double s = dot(disp, l->g_norm);
        *energy += s * l->mass * G;
    }
}

void calculate_neg_force(data *d, gsl_vector *neg_force) {

    wheel *w = d->wheel;
    load *l = d->load;
    double force = 0;
    gsl_vector_set_zero(neg_force);

    for (int i = 0; i < w->n_holes; ++i) {
        xy *spoke = &d->spoke[i];
        double l = length(*spoke);
        // there's a missing - sign because derivative is -force
        double fx = w->e_spoke * d->spoke_strain[i] * spoke->x / l;
        double fy = w->e_spoke * d->spoke_strain[i] * spoke->y / l;
        gsl_vector_set(neg_force, 2*i+X, gsl_vector_get(neg_force, 2*i+X) + fx);
        gsl_vector_set(neg_force, 2*i+Y, gsl_vector_get(neg_force, 2*i+Y) + fy);
//        ludebug(dbg, "Forces due to spoke %d: %g, %g", i, fx, fy);
    }

    for (int i = 0; i < w ->n_holes; ++i) {
        xy *chord = &d->chord[i];
        double l = length(*chord);
        double fx = w->e_rim * d->chord_strain[i] * chord->x / l;
        double fy = w->e_rim * d->chord_strain[i] * chord->y / l;
        gsl_vector_set(neg_force, 2*i+X, gsl_vector_get(neg_force, 2*i+X) - fx);
        gsl_vector_set(neg_force, 2*i+Y, gsl_vector_get(neg_force, 2*i+Y) - fy);
        int j = (i - 1 + w->n_holes) % w->n_holes;
        gsl_vector_set(neg_force, 2*j+X, gsl_vector_get(neg_force, 2*j+X) + fx);
        gsl_vector_set(neg_force, 2*j+Y, gsl_vector_get(neg_force, 2*j+Y) + fy);
//        ludebug(dbg, "Forces due to chord %d: %g, %g", i, fx, fy);
    }

    if (l) {
        gsl_vector_set(neg_force, 2*l->i_rim+X, gsl_vector_get(neg_force, 2*l->i_rim+X) - l->mass * G * l->g_norm.x);
        gsl_vector_set(neg_force, 2*l->i_rim+Y, gsl_vector_get(neg_force, 2*l->i_rim+Y) - l->mass * G * l->g_norm.y);
    }

    for (int i = 0; i < w->n_holes; ++i) {
        double fx = gsl_vector_get(neg_force, 2*i+X), fy = gsl_vector_get(neg_force, 2*i+Y);
        force += sqrt(fx * fx + fy * fy);
    }
    ludebug(dbg, "Force: %g", force);
}

double energy(const gsl_vector *rim, void *params) {
    data *d = (data*)params;
    double energy;
    calculate_data(rim, d);
    calculate_energy(d, &energy);
    ludebug(dbg, "Energy: %g", energy);
    return energy;
}

void neg_force(const gsl_vector *rim, void *params, gsl_vector *neg_force) {
    data *d = (data*)params;
    calculate_data(rim, d);
    calculate_neg_force(d, neg_force);
}

void energy_and_neg_force(const gsl_vector *rim, void *params, double *energy, gsl_vector *neg_force) {
    data *d = (data*)params;
    calculate_data(rim, d);
    calculate_energy(d, energy);
    calculate_neg_force(d, neg_force);
}

int relax(wheel *wheel, load *load) {

    LU_STATUS
    const gsl_multimin_fdfminimizer_type *T;
    gsl_vector *rim = NULL;
    gsl_multimin_fdfminimizer *s = NULL;
    gsl_multimin_function_fdf callbacks;
    data *d = NULL;
    int gsl_status = GSL_CONTINUE;

    LU_ALLOC(dbg, d, 1);
    d->wheel = wheel;
    d->load = load;
    LU_ALLOC(dbg, d->spoke, wheel->n_holes);
    LU_ALLOC(dbg, d->spoke_strain, wheel->n_holes);
    LU_ALLOC(dbg, d->chord, wheel->n_holes);
    LU_ALLOC(dbg, d->chord_strain, wheel->n_holes);

    callbacks.n = 2 * wheel->n_holes;
    callbacks.params = d;
    callbacks.f = energy;
    callbacks.df = neg_force;
    callbacks.fdf = energy_and_neg_force;

    rim = gsl_vector_alloc(2 * wheel->n_holes);
    for (int i = 0; i < wheel->n_holes; ++i) {
        gsl_vector_set(rim, 2*i+X, wheel->rim[i].x);
        gsl_vector_set(rim, 2*i+Y, wheel->rim[i].y);
    }

    double e = energy(rim, d);
    luinfo(dbg, "Initial energy: %g", e);

//    T = gsl_multimin_fdfminimizer_steepest_descent;
    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc(T, 2 * wheel->n_holes);
    luinfo(dbg, "Method %s", gsl_multimin_fdfminimizer_name(s));
    gsl_multimin_fdfminimizer_set(s, &callbacks, rim, 1e-1, 0.1);

    for (int iter = 0; iter < 1000 && gsl_status == GSL_CONTINUE; ++iter) {
        ludebug(dbg, "Iteration %d", iter);
        gsl_status = gsl_multimin_fdfminimizer_iterate(s);
        if (gsl_status == GSL_ENOPROG) {
            luwarn(dbg, "Cannot progress");
        // below not necessary - could consider restart
//        } else if (!gsl_status && iter < 1) {
//            ludebug(dbg, "Forcing new iteration");
//            gsl_status = GSL_CONTINUE;
        } else if (!gsl_status) {
            gsl_status = gsl_multimin_test_gradient(s->gradient, 1e4);
            if (gsl_status == GSL_SUCCESS) luinfo(dbg, "Minimum energy: %g", s->f);
        }
    }
    LU_ASSERT(gsl_status == GSL_SUCCESS, LU_ERR, dbg, "Equilibrium not found")

    double shift = 0.0;
    for (int i = 0; i < wheel->n_holes; ++i) {
        double x = gsl_vector_get(s->x, 2*i+X), y = gsl_vector_get(s->x, 2*i+Y);
        double dx = wheel->rim[i].x - x, dy = wheel->rim[i].y - y;
        shift += sqrt(dx * dx + dy * dy);
        wheel->rim[i].x = x;
        wheel->rim[i].y = y;
    }
    shift = shift / wheel->n_holes;
    ludebug(dbg, "Average shift %gmm", shift);

LU_CLEANUP
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(rim);
    if (d) {
        free(d->spoke);
        free(d->spoke_strain);
        free(d->chord);
        free(d->chord_strain);
        free(d);
    }
    LU_RETURN
}

#define TARGET_WOBBLE 1e-6

int true(wheel *w) {

    LU_STATUS

    LU_CHECK(relax(w, NULL))

    double r_target = 0, wobble = 2 * TARGET_WOBBLE;
    for (int i = 0; i < w->n_holes; ++i) r_target += length(w->rim[i]);
    r_target = r_target / w->n_holes;

    while(1) {
        wobble = 0;
        for (int i = 0; i < w->n_holes; ++i) wobble += fabs(length(w->rim[i]) - r_target);
        wobble = wobble / w->n_holes;
        luinfo(dbg, "Average radial wobble %gmm", wobble);
        if (wobble <= TARGET_WOBBLE) break;
        // this correction never used, so may be wrong...
        for (int i = 0; i < w->n_holes; ++i) {
            double delta = (length(w->rim[i]) - r_target) / r_target;
            w->l_spoke[i] = w->l_spoke[i] * (1 - 0.5 * delta);
        }
        LU_CHECK(relax(w, NULL));
    }

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
    l->g_norm.x = 1/sqrt(2);
    l->g_norm.y = 1/sqrt(2);
    l->i_rim = 4;
    l->mass = 100;
    l->start = wheel->rim[l->i_rim];

    LU_CHECK(relax(wheel, l))

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
    LU_CHECK(deform(wheel))
    LU_CHECK(make_path(dbg, pattern, &path))
    plot_wheel(wheel, path);

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
        LU_CHECK(stress(argv[1]));
    }

LU_CLEANUP
    if (dbg) status = dbg->free(&dbg, status);
    return status;
}
