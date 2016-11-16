
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>

#include "gsl/gsl_multimin.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_rng.h"

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/strings.h"
#include "lu/dynamic_memory.h"

#include "lib.h"
#include "wheel.h"

lulog *dbg = NULL;
gsl_rng *rng = NULL;
volatile int sig_exit = 0;

#define X 0
#define Y 1
#define G 9.8

#define MAX_ITER_OUTER 10000
#define MAX_ITER_INNER 10000
#define N_DEFORM 1
#define MAX_SIZE 1e-8
#define MAX_FORCE 1


/*
 * this coded does finite element analysis of wheels using a very simple
 * physical model (described somewhere in an email or online post by jobst)
 * where the rim is "hinged" at each spoke hole.
 *
 * unfortunately, for two leading two following, that model leads to a
 * non-physical solution, where the rim folds back on itself.  also, the
 * code below is extremely slow to converge - perhaps because i don't
 * know how best to propagate forces around the rim when not in
 * equilibrium.
 *
 * so i am not developing this further.
 */


struct data;

typedef void coeff_to_rim(const gsl_vector *coeff, struct data *d);

typedef struct data {
    wheel *wheel;          // wheel in initial state
    coeff_to_rim *to_rim;  // mapping from coeff to rim
    xy *offset;            // radial and tangential offsets from coeffs
    xy *rim;               // rim locations (wheel->rim + offset)
    double *spoke_extn;    // extension of spoke (mm) (rim index)
    xy *spoke;             // vector along spoke (rim - hub) (rim index)
    double *chord_extn;    // extension of rim segment (mm)
    xy *chord;             // vector along rim segment
    load *load;            // additional load
} data;

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

double eval_fourier_coeff(int i, int j, int n, double a) {
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

void fourier_coeff_to_rim(const gsl_vector *coeff, data *d) {

    wheel *w = d->wheel;

    for (int i = 0; i < w->n_holes; ++i) {
        for (int j = 0; j < w->n_holes; ++j) {
            d->offset[j].x += eval_fourier_coeff(i, j, w->n_holes, gsl_vector_get(coeff, 2*i+X));
            d->offset[j].y += eval_fourier_coeff(i, j, w->n_holes, gsl_vector_get(coeff, 2*i+Y));
        }
    }

    for (int i = 0; i < w->n_holes; ++i) {
        double dr = d->offset[i].x, dt = d->offset[i].y;
        double theta = atan2(w->rim[i].y, w->rim[i].x);
        d->rim[i].x = w->rim[i].x + dr * cos(theta) + dt * sin(theta);
        d->rim[i].y = w->rim[i].y + dr * sin(theta) - dt * cos(theta);
    }

}

void xy_coeff_to_rim(const gsl_vector *coeff, data *d) {

    wheel *w = d->wheel;

    for (int i = 0; i < w->n_holes; ++i) {
        d->offset[i].x = gsl_vector_get(coeff, 2*i+X);
        d->offset[i].y = gsl_vector_get(coeff, 2*i+Y);
    }

    for (int i = 0; i < w->n_holes; ++i) {
        d->rim[i] = add(w->rim[i], d->offset[i]);
    }

}


void calculate_data(const gsl_vector *coeff, data *d) {

    wheel *w = d->wheel;

    memset(d->offset, 0, w->n_holes * sizeof(*d->offset));
    memset(d->rim, 0, w->n_holes * sizeof(*d->rim));
    memset(d->spoke_extn, 0, w->n_holes * sizeof(*d->spoke_extn));
    memset(d->spoke, 0, w->n_holes * sizeof(*d->spoke));
    memset(d->chord_extn, 0, w->n_holes * sizeof(*d->chord_extn));
    memset(d->chord, 0, w->n_holes * sizeof(*d->chord));

    d->to_rim(coeff, d);

    for (int i = 0; i < w->n_holes; ++i) {
        d->spoke[i] = sub(d->rim[i], w->hub[w->rim_to_hub[i]]);
        d->spoke_extn[i] = length(d->spoke[i]) - w->l_spoke[w->rim_to_hub[i]];
//        ludebug(dbg, "Spoke %d extended by %gmm", i, d->spoke_extn[i]);
    }

    for (int after = 0; after < w ->n_holes; ++after) {
        int before = (after - 1 + w->n_holes) % w->n_holes;
        d->chord[after] = sub(d->rim[after], d->rim[before]);
        d->chord_extn[after] = length(d->chord[after]) - w->l_chord;
//        ludebug(dbg, "Chord %d extended by %gmm", after, d->chord_extn[after]);
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
//        ludebug(dbg, "Load moved by %gmm (%g,%g)", s, disp.x, disp.y);
        *energy -= s * l->mass * G * 1e-3;  // mm -> m
    }
}

double sum_force(wheel *w, gsl_vector *neg_force) {
    double total = 0;
    for (int i = 0; i < w->n_holes; ++i) {
        double fx = gsl_vector_get(neg_force, 2*i+X), fy = gsl_vector_get(neg_force, 2*i+Y);
        total += sqrt(fx * fx + fy * fy);
    }
    return total;
}

void calculate_neg_force(data *d, gsl_vector *neg_force) {

    wheel *w = d->wheel;
    load *l = d->load;
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
}

double energy(const gsl_vector *coeff, void *params) {
    data *d = (data*)params;
    double energy;
    calculate_data(coeff, d);
    calculate_energy(d, &energy);
//    ludebug(dbg, "Energy: %g", energy);
    return energy;
}

void neg_force(const gsl_vector *coeff, void *params, gsl_vector *neg_force) {
    data *d = (data*)params;
    calculate_data(coeff, d);
    calculate_neg_force(d, neg_force);
}

void energy_and_neg_force(const gsl_vector *coeff, void *params, double *energy, gsl_vector *neg_force) {
    data *d = (data*)params;
    calculate_data(coeff, d);
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

void log_energy(data *d, const char *msg, double *final_energy, double *final_force) {
    double energy, total;
    gsl_vector *neg_force = gsl_vector_alloc(d->wheel->n_holes * 2);
    calculate_energy(d, &energy);
    calculate_neg_force(d, neg_force);
    total = sum_force(d->wheel, neg_force);
    luinfo(dbg, "%s: Energy %g, Force %g", msg, energy, total);
    if (final_energy) *final_energy = energy;
    if (final_force) *final_force = total;
    free(neg_force);
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

int relax_f_fourier(wheel *wheel, load *load, double *final_energy, double *final_force) {

    LU_STATUS
    gsl_vector *coeff = NULL, *step = NULL;
    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s = NULL;
    gsl_multimin_function callbacks;
    data *d = NULL;
    int gsl_status = GSL_CONTINUE;

    LU_CHECK(alloc_data(&d, wheel, load))
    d->to_rim = &fourier_coeff_to_rim;
    alloc_coeff(&coeff, wheel);

    calculate_data(coeff, d);
    log_energy(d, "Before relax f Fourier", NULL, NULL);

    callbacks.n = 2 * wheel->n_holes;
    callbacks.params = d;
    callbacks.f = energy;

    T = gsl_multimin_fminimizer_nmsimplex;
    s = gsl_multimin_fminimizer_alloc(T, 2 * wheel->n_holes);
//    luinfo(dbg, "Method %s", gsl_multimin_fminimizer_name(s));
    step = gsl_vector_alloc(2 * wheel->n_holes);
    gsl_vector_set_all(step, 1e-4);
    gsl_multimin_fminimizer_set(s, &callbacks, coeff, step);

    for (int iter = 0; iter < MAX_ITER_OUTER && gsl_status == GSL_CONTINUE && !sig_exit; ++iter) {
//        ludebug(dbg, "Iteration %d", iter);
        gsl_status = gsl_multimin_fminimizer_iterate(s);
        if (gsl_status == GSL_ENOPROG) {
            luwarn(dbg, "Cannot progress");
        } else if (!gsl_status) {
            double size = gsl_multimin_fminimizer_size(s);
            if (!iter || !(iter & (iter - 1))) ludebug(dbg, "Size %g (%d)", size, iter);
            gsl_status = gsl_multimin_test_size(size, MAX_SIZE);
            if (gsl_status == GSL_SUCCESS) luinfo(dbg, "Minimum energy: %g", s->fval);
        }
    }

    log_energy(d, "After relax f Fourier", final_energy, final_force);
    update_rim(s->x, d, wheel);

LU_CLEANUP
    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(step);
    gsl_vector_free(coeff);
    free_data(d);
    LU_RETURN
}

double vec_len(gsl_vector *v) {
    double l = 0;
    for (int i = 0; i < v->size; ++i) l = l + gsl_vector_get(v, i) * gsl_vector_get(v, i);
    return sqrt(l);
}

int relax_fdf_xy(wheel *wheel, load *load, double *final_energy, double *final_force) {

    LU_STATUS
    gsl_vector *coeff = NULL, *step = NULL;
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s = NULL;
    gsl_multimin_function_fdf callbacks;
    data *d = NULL;
    int gsl_status = GSL_CONTINUE;

    LU_CHECK(alloc_data(&d, wheel, load))
    d->to_rim = &xy_coeff_to_rim;
    alloc_coeff(&coeff, wheel);

    calculate_data(coeff, d);
    log_energy(d, "Before relax fdf xy", NULL, NULL);

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
    gsl_multimin_fdfminimizer_set(s, &callbacks, coeff, 1e-4, 1e-3);

    for (int iter = 0; iter < MAX_ITER_OUTER && gsl_status == GSL_CONTINUE && !sig_exit; ++iter) {
//        ludebug(dbg, "Iteration %d", iter);
        gsl_status = gsl_multimin_fdfminimizer_iterate(s);
        if (gsl_status == GSL_ENOPROG) {
            luwarn(dbg, "Cannot progress");
        } else if (!gsl_status) {
            gsl_status = gsl_multimin_test_gradient(s->gradient, 1e-4);
            if (!iter || !(iter & (iter - 1))) ludebug(dbg, "Gradient %g (%d)", vec_len(s->gradient), iter);
            if (gsl_status == GSL_SUCCESS) luinfo(dbg, "Minimum energy: %g", s->f);
        }
    }

    log_energy(d, "After relax fdf xy", final_energy, final_force);
    update_rim(s->x, d, wheel);

LU_CLEANUP
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(coeff);
    free_data(d);
    LU_RETURN
}

int relax(wheel *w, load *l, int n) {
    LU_STATUS
    double final_energy, final_force;
    for (int i = 0; i < n; ++i) {
        ludebug(dbg, "Relax %d/%d", i , n);
        LU_CHECK(relax_fdf_xy(w, l, NULL, NULL))
        LU_CHECK(relax_f_fourier(w, l, &final_energy, &final_force))
        if (final_force <= MAX_FORCE) break;
    }
    LU_NO_CLEANUP
}

#define TARGET_WOBBLE 3e-3
#define DAMP_TRUE 0.5
#define DAMP_TENSION 0.5

int true(wheel *w) {

    LU_STATUS

    LU_CHECK(relax(w, NULL, MAX_ITER_INNER))
    LU_CHECK(dump_wheel(dbg, w, "untrue"))

    while(1) {

        double r_target = 0;
        for (int i = 0; i < w->n_holes; ++i) r_target += length(w->rim[i]);
        r_target = r_target / w->n_holes;
        double wobble = 0, shift = 0;
        for (int i = 0; i < w->n_holes; ++i) wobble += fabs(length(w->rim[i]) - r_target);
        wobble /= w->n_holes;
        luinfo(dbg, "Average radial wobble %gmm", wobble);
        if (wobble <= TARGET_WOBBLE || sig_exit) break;
        for (int i = 0; i < w->n_holes; ++i) {
            double l = length(w->rim[i]);
            if (l > r_target) {
                double excess = DAMP_TRUE * (l - r_target);
                w->l_spoke[w->rim_to_hub[i]] -= excess;
                shift += excess;
            }
        }
        luinfo(dbg, "Total correction %gmm", shift);
        LU_CHECK(relax(w, NULL, MAX_ITER_INNER));

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
        LU_CHECK(relax(w, NULL, MAX_ITER_INNER));
    }

    LU_CHECK(dump_wheel(dbg, w, "true"))
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

    LU_CHECK(dump_wheel(dbg, wheel, "laced"))

LU_CLEANUP
    free(tension);
    LU_RETURN
}

int deform(wheel *wheel, load *l) {

    LU_STATUS

    for (int i = 0; i < N_DEFORM; ++i) {
        l->mass = 10 * (float)(i + 1) / N_DEFORM;
        ludebug(dbg,"Mass %gkg", l->mass);
        l->start = wheel->rim[l->i_rim];
        LU_CHECK(relax(wheel, l, MAX_ITER_INNER))
    }

    LU_NO_CLEANUP
}

int alloc_load(load **l) {

    LU_STATUS

    LU_ALLOC(dbg, *l, 1)
    (*l)->g_norm.x = 1/sqrt(2);
    (*l)->g_norm.y = -1/sqrt(2);
//    (*l)->g_norm.x = 0;
//    (*l)->g_norm.y = -1;
    (*l)->i_rim = 0;

    LU_NO_CLEANUP
}

int stress(const char *pattern) {

    LU_STATUS
    int *offsets = NULL, length = 0, holes = 0, padding;
    char type, *path = NULL;
    wheel *wheel = NULL, *original = NULL;
    load *load = NULL;

    luinfo(dbg, "Pattern '%s'", pattern);
    LU_CHECK(unpack(dbg, pattern, &offsets, &length, &type, &padding))
    LU_ASSERT(strchr("AB", type), LU_ERR, dbg, "Only symmetric types supported")
    LU_CHECK(dump_pattern(dbg, offsets, length))
    LU_CHECK(rim_size(dbg, length, &holes))
    LU_CHECK(make_wheel(dbg, offsets, length, holes, padding, type, pattern, &wheel))
    LU_CHECK(lace(wheel))
    LU_CHECK(true(wheel))
    LU_CHECK(copy_wheel(dbg, wheel, &original))
    LU_CHECK(alloc_load(&load))
    LU_CHECK(deform(wheel, load))
    plot_multi_deform(dbg, original, wheel, load, pattern);
//    plot_wheel(original, path);

LU_CLEANUP
    free_wheel(original);
    free_wheel(wheel);
    free(path);
    free(offsets);
    free(load);
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
