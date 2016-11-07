
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>

#include "cairo/cairo.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_rng.h"

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

#define ZOOM 100000
#define MAX_ITER 100


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
    double *l_spoke;  // unloaded, indexed by hub
    double tension;   // target tension for spokes
    double e_rim;
    double e_spoke;
} wheel;

typedef struct {
    int i_rim;    // location of load
    double mass;  // load mass
    xy g_norm;    // direction in which gravity acts
    xy start;     // initial location of load (rim)
    xy end;       // final location of load (rim)
} load;

typedef struct {
    wheel *wheel;        // wheel in initial state
    xy *offset;          // radial and tangential offsets from coeffs
    xy *rim;             // rim locations (wheel->rim + offset)
    double *spoke_extn;  // extension of spoke (mm) (rim index)
    xy *spoke;           // vector along spoke (rim - hub) (rim index)
    double *chord_extn;  // extension of rim segment (mm)
    xy *chord;           // vector along rim segment
    load *load;          // additional load
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

xy zoom(xy a, xy b, double r_scale) {
    xy d = sub(b, a);
    double r = sqrt(d.x * d.x + d.y * d.y), theta = atan2(d.y, d.x);
    r *= r_scale;
    d.x = r * cos(theta); d.y = r * sin(theta);
    return add(a, d);
}

void draw_deform(cairo_t *cr, wheel *a, wheel *b, double r_scale) {
    for (int i = 0; i < a->n_holes; ++i) {
        int j = a->hub_to_rim[i];
        xy hub_a = a->hub[i], hub_b = b->hub[i];
        xy rim1_a = a->rim[j], rim1_b = b->rim[j];
        int k = (j + 1) % a->n_holes;
        xy rim2_a = a->rim[k], rim2_b = b->rim[k];
        draw_line_xy(cr, zoom(hub_a, hub_b, r_scale), zoom(rim1_a, rim1_b, r_scale));
        draw_line_xy(cr, zoom(rim1_a, rim1_b, r_scale), zoom(rim2_a, rim2_b, r_scale));
    }
}

void plot_deform(wheel *a, wheel *b, const char *path) {
    cairo_surface_t *surface = NULL;
    cairo_t *cr = NULL;
    open_plot(a, 500, 500, a->r_rim / 100, &cr, &surface);
    cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
    draw_wheel(cr, a);
    cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
    draw_deform(cr, a, b, ZOOM);
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

double eval_coeff(int i, int j, int n, double a) {
    double shift = a;
    if (i) {
        int n_cycles = (i + 1) / 2;  // i = 1,2 -> 1; i = n-1 -> n/2
        double x = j / (float)n;  // 0-1
        double y = 2 * M_PI * x * n_cycles;
        if (i % 2) {
            // when i = n-1 we are here with n_cycles = n/2
            // and get 1,-1,1,-1... for j = 0...n-1
            // eg j = 1, x = 1/n, y = pi
            shift = a * cos(y);
        } else {
            shift = a * sin(y);
        }
    }
//    ludebug(dbg, "Shift %g for coeff %d at rim %d/%d, coeff %g", shift, i, j, n, a);
    return shift;
}

void calculate_data(const gsl_vector *coeff, data *d) {

    wheel *w = d->wheel;

    memset(d->offset, 0, w->n_holes * sizeof(*d->offset));
    memset(d->rim, 0, w->n_holes * sizeof(*d->rim));
    memset(d->spoke_extn, 0, w->n_holes * sizeof(*d->spoke_extn));
    memset(d->spoke, 0, w->n_holes * sizeof(*d->spoke));
    memset(d->chord_extn, 0, w->n_holes * sizeof(*d->chord_extn));
    memset(d->chord, 0, w->n_holes * sizeof(*d->chord));

    for (int i = 0; i < w->n_holes; ++i) {
        for (int j = 0; j < w->n_holes; ++j) {
            d->offset[j].x += eval_coeff(i, j, w->n_holes, gsl_vector_get(coeff, 2*i+X));
            d->offset[j].y += eval_coeff(i, j, w->n_holes, gsl_vector_get(coeff, 2*i+Y));
        }
    }

    for (int i = 0; i < w->n_holes; ++i) {
        double dr = d->offset[i].x, dt = d->offset[i].y;
        double theta = atan2(w->rim[i].y, w->rim[i].x);
        d->rim[i].x = w->rim[i].x + dr * cos(theta) + dt * sin(theta);
        d->rim[i].y = w->rim[i].y + dr * sin(theta) - dt * cos(theta);
    }

    for (int i = 0; i < w->n_holes; ++i) {
        d->spoke[i] = sub(d->rim[i], w->hub[w->rim_to_hub[i]]);
        d->spoke_extn[i] = length(d->spoke[i]) - w->l_spoke[w->rim_to_hub[i]];
        ludebug(dbg, "Spoke %d extended by %gmm", i, d->spoke_extn[i]);
    }

    for (int after = 0; after < w ->n_holes; ++after) {
        int before = (after - 1 + w->n_holes) % w->n_holes;
        d->chord[after] = sub(d->rim[after], d->rim[before]);
        d->chord_extn[after] = length(d->chord[after]) - w->l_chord;
        ludebug(dbg, "Chord %d extended by %gmm", after, d->chord_extn[after]);
    }

    if (d->load) {
        d->load->end = d->rim[d->load->i_rim];
    }
}

void calculate_energy(data *d, double *energy) {

    wheel *w = d->wheel;
    load *l = d->load;
    *energy = 0;

    for (int i = 0; i < w->n_holes; ++i) {
        double extn = d->spoke_extn[i];
        // 1e-3 since mm
        *energy += 1e-3 * w->e_spoke * extn * extn / (2 * w->l_spoke[w->rim_to_hub[i]]);
    }
    for (int i = 0; i < w ->n_holes; ++i) {
        double extn = d->chord_extn[i];
        *energy += 1e-3 * w->e_rim * extn * extn / (2 * w->l_chord);
    }
    if (l) {
        xy disp = sub(l->end, l->start);
        double s = dot(disp, l->g_norm);
        ludebug(dbg, "Load moved by %gmm (%g,%g)", s, disp.x, disp.y);
        *energy -= s * l->mass * G * 1e-3;  // mm -> m
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
        double l0 = w->l_spoke[w->rim_to_hub[i]];
        // there's a missing - sign because derivative is -force
        double fx = w->e_spoke * d->spoke_extn[i] * spoke->x / (l * l0);
        double fy = w->e_spoke * d->spoke_extn[i] * spoke->y / (l * l0);
        gsl_vector_set(neg_force, 2*i+X, gsl_vector_get(neg_force, 2*i+X) + fx);
        gsl_vector_set(neg_force, 2*i+Y, gsl_vector_get(neg_force, 2*i+Y) + fy);
//        ludebug(dbg, "Forces due to spoke %d: %g, %g", i, fx, fy);
    }

    for (int i = 0; i < w ->n_holes; ++i) {
        xy *chord = &d->chord[i];
        double l = length(*chord);
        double l0 = w->l_chord;
        double fx = w->e_rim * d->chord_extn[i] * chord->x / (l * l0);
        double fy = w->e_rim * d->chord_extn[i] * chord->y / (l * l0);
        gsl_vector_set(neg_force, 2*i+X, gsl_vector_get(neg_force, 2*i+X) + fx);
        gsl_vector_set(neg_force, 2*i+Y, gsl_vector_get(neg_force, 2*i+Y) + fy);
        int j = (i - 1 + w->n_holes) % w->n_holes;
        gsl_vector_set(neg_force, 2*j+X, gsl_vector_get(neg_force, 2*j+X) - fx);
        gsl_vector_set(neg_force, 2*j+Y, gsl_vector_get(neg_force, 2*j+Y) - fy);
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

double energy(const gsl_vector *coeff, void *params) {
    data *d = (data*)params;
    double energy;
    calculate_data(coeff, d);
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

int alloc_data(data **d, wheel *w, load *l) {
    LU_STATUS
    LU_ALLOC(dbg, *d, 1);
    (*d)->wheel = w;
    (*d)->load = l;
    LU_ALLOC(dbg, (*d)->offset, w->n_holes)
    LU_ALLOC(dbg, (*d)->rim, w->n_holes)
    LU_ALLOC(dbg, (*d)->spoke, w->n_holes);
    LU_ALLOC(dbg, (*d)->spoke_extn, w->n_holes);
    LU_ALLOC(dbg, (*d)->chord, w->n_holes);
    LU_ALLOC(dbg, (*d)->chord_extn, w->n_holes);
    LU_NO_CLEANUP
}

void free_data(data *d) {
    if (d) {
        free(d->offset);
        free(d->rim);
        free(d->spoke);
        free(d->spoke_extn);
        free(d->chord);
        free(d->chord_extn);
        free(d);
    }
}

void alloc_coeff(gsl_vector **coeff, wheel *w) {
    *coeff = gsl_vector_calloc(2 * w->n_holes);
}

void log_energy(wheel *w, load *l, const char *msg) {
    gsl_vector *coeff = NULL;
    data *d = NULL;
    alloc_coeff(&coeff, w);
    alloc_data(&d, w, l);
    double e = energy(coeff, d);
    luinfo(dbg, "%s: %g", msg, e);
    free_data(d);
    gsl_vector_free(coeff);
}

void update_rim(gsl_vector *coeff, data *d, wheel *w) {
    double shift = 0.0;
    calculate_data(coeff, d);
    for (int i = 0; i < w->n_holes; ++i) {
        xy delta = sub(d->rim[i], w->rim[i]);
        shift += length(delta);
        w->rim[i] = d->rim[i];
    }
    shift = shift / w->n_holes;
    ludebug(dbg, "Average shift %gmm", shift);
}

int relax_f(wheel *wheel, load *load) {

    LU_STATUS
    gsl_vector *coeff = NULL, *step = NULL;
    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s = NULL;
    gsl_multimin_function callbacks;
    data *d = NULL;
    int gsl_status = GSL_CONTINUE;

    LU_CHECK(alloc_data(&d, wheel, load))
    alloc_coeff(&coeff, wheel);

    callbacks.n = 2 * wheel->n_holes;
    callbacks.params = d;
    callbacks.f = energy;

    T = gsl_multimin_fminimizer_nmsimplex;
    s = gsl_multimin_fminimizer_alloc(T, 2 * wheel->n_holes);
//    luinfo(dbg, "Method %s", gsl_multimin_fminimizer_name(s));
    step = gsl_vector_alloc(2 * wheel->n_holes);
    gsl_vector_set_all(step, 1e-4);
    gsl_multimin_fminimizer_set(s, &callbacks, coeff, step);

    for (int iter = 0; iter < MAX_ITER && gsl_status == GSL_CONTINUE && !sig_exit; ++iter) {
        ludebug(dbg, "Iteration %d", iter);
        gsl_status = gsl_multimin_fminimizer_iterate(s);
        if (gsl_status == GSL_ENOPROG) {
            luwarn(dbg, "Cannot progress");
        } else if (!gsl_status) {
            double size = gsl_multimin_fminimizer_size(s);
            ludebug(dbg, "Size %g", size);
            gsl_status = gsl_multimin_test_size(size, 1e-10);
            if (gsl_status == GSL_SUCCESS) luinfo(dbg, "Minimum energy: %g", s->fval);
        }
    }

    update_rim(s->x, d, wheel);

LU_CLEANUP
    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(step);
    gsl_vector_free(coeff);
    free_data(d);
    LU_RETURN
}

int relax_fdf(wheel *wheel, load *load) {

    LU_STATUS
    gsl_vector *rim = NULL, *step = NULL;
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s = NULL;
    gsl_multimin_function_fdf callbacks;
    data *d = NULL;
    int gsl_status = GSL_CONTINUE;

    LU_CHECK(alloc_data(&d, wheel, load))
    alloc_coeff(&rim, wheel);

    callbacks.n = 2 * wheel->n_holes;
    callbacks.params = d;
    callbacks.f = energy;
    callbacks.df = neg_force;
    callbacks.fdf = energy_and_neg_force;

    T = gsl_multimin_fdfminimizer_steepest_descent;
//    T = gsl_multimin_fdfminimizer_conjugate_fr;
//    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc(T, 2 * wheel->n_holes);
//    luinfo(dbg, "Method %s", gsl_multimin_fdfminimizer_name(s));
    gsl_multimin_fdfminimizer_set(s, &callbacks, rim, 1e-4, 1e-3);

    for (int iter = 0; iter < MAX_ITER && gsl_status == GSL_CONTINUE && !sig_exit; ++iter) {
        ludebug(dbg, "Iteration %d", iter);
        gsl_status = gsl_multimin_fdfminimizer_iterate(s);
        if (gsl_status == GSL_ENOPROG) {
            luwarn(dbg, "Cannot progress");
        } else if (!gsl_status) {
            gsl_status = gsl_multimin_test_gradient(s->gradient, 1e-4);
            if (gsl_status == GSL_SUCCESS) luinfo(dbg, "Minimum energy: %g", s->f);
        }
    }

    update_rim(s->x, d, wheel);

LU_CLEANUP
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(rim);
    free_data(d);
    LU_RETURN
}

int relax(wheel *w, load *l, int n) {
    LU_STATUS
    for (int i = 0; i < n; ++i) {
        if (i) log_energy(w, l, "Intermediate energy");
//        LU_CHECK(relax_fdf(w, l))
        LU_CHECK(relax_f(w, l))
    }
    LU_NO_CLEANUP
}

#define TARGET_WOBBLE 1e-3
#define DAMP_TRUE 0.5
#define DAMP_TENSION 0.5

int true(wheel *w) {

    LU_STATUS

    log_energy(w, NULL, "Energy before truing");
    LU_CHECK(relax(w, NULL, 1))

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
        LU_CHECK(relax(w, NULL, 1));

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
        LU_CHECK(relax(w, NULL, 1));
    }

    luinfo(dbg, "True!");
    log_energy(w, NULL, "Energy after truing");

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
    l->g_norm.y = -1/sqrt(2);
//    l->g_norm.x = 0;
//    l->g_norm.y = -1;
    l->i_rim = 0;
    l->start = wheel->rim[l->i_rim];
    l->mass = 10;
    ludebug(dbg,"Mass %gkg", l->mass);
    log_energy(wheel, l, "Energy before deforming");
    LU_CHECK(relax(wheel, l, 2))
    log_energy(wheel, l, "Energy after deforming");

LU_CLEANUP
    free(l);
    LU_RETURN
}

int stress(const char *pattern) {

    LU_STATUS
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
    luwarn(dbg, "Handler called with %d", sig);
    sig_exit = 1;
}

int set_handler() {
    LU_STATUS
    LU_ASSERT(!signal(SIGINT, &new_handler), LU_ERR, dbg, "Could not set handler")
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
        LU_ASSERT(rng = gsl_rng_alloc(gsl_rng_mt19937), LU_ERR, dbg, "Could not create PRNG")
        LU_CHECK(stress(argv[1]))
    }

LU_CLEANUP
    if (dbg) status = dbg->free(&dbg, status);
    if (rng) gsl_rng_free(rng);
    return status;
}
