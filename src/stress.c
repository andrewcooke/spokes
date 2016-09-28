
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

    for (int hub = 0; hub < wheel->n_holes; ++hub) {
        int rim = wheel->hub_to_rim[hub];
        draw_line_xy(cr, wheel->hub[hub], wheel->rim[rim]);
        draw_line_xy(cr, wheel->hub[hub], wheel->hub[(hub+1) % wheel->n_holes]);
        draw_line_xy(cr, wheel->rim[rim], wheel->rim[(rim+1) % wheel->n_holes]);
    }

    cairo_surface_write_to_png(surface, path);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);

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

xy polar_scale(xy components, xy location, double radial, double tangential, double *acc_r, double *acc_t) {
    xy n = norm(location);
    double r = radial * (components.x * n.x + components.y * n.y);
    double t = tangential * (components.x * n.y - components.x * n.x);
    *acc_r = *acc_r + fabs(r);
    *acc_t = *acc_t + fabs(t);
    double theta = atan2(location.y, location.x);
    xy c;
    c.x = r * cos(theta) + t * sin(theta);
    c.y = r * sin(theta) - t * cos(theta);
    return c;
}

int relax(wheel *wheel, double damping) {

    LU_STATUS
    double acc_r = 0, acc_t = 0, acc_f = 0;

    xy *forces = NULL;  // force at hole[index]
    LU_ALLOC(dbg, forces, wheel->n_holes)

    // accumulate total force

    // force on hole i_rim from spoke at hub i_hub
    for (int i_rim = 0; i_rim < wheel->n_holes; ++i_rim) {
        int i_hub = wheel->rim_to_hub[i_rim];
        xy hub = wheel->hub[i_hub];
        xy rim = wheel->rim[i_rim];
        xy spoke = sub(hub, rim);  // points inwards towards hub
        double l = length(spoke);
        // 1/l for strain, 1/l for normalization of spoke direction
        // stretched here, so l - l_spoke
//        forces[i_rim] = scalar_mult((l - wheel->l_spoke[i_hub]) * wheel->e_spoke / (l * l), spoke);
        forces[i_rim] = scalar_mult((l - wheel->l_spoke[i_hub]) / l, spoke);
    }

    // forces on hole i from rim chord to either side
    // we are looking at the chord from i_after-1 (i_before) to i_after
    xy before = wheel->rim[wheel->n_holes-1];
    for (int i_after = 0; i_after < wheel->n_holes; ++i_after) {
        xy after = wheel->rim[i_after];
        xy chord = sub(after, before);  // points from i_before towards i_after
        double l = length(chord);
        // compressive here, so l_chord - l
//        xy force = scalar_mult((wheel->l_chord - l) * wheel->e_rim / (l * l), chord);  // points towards i_after
        xy force = scalar_mult((wheel->l_chord - l) / l, chord);  // points towards i_after
        forces[i_after] = add(forces[i_after], force);
        int i_before = (i_after-1+wheel->n_holes) % wheel->n_holes;
        forces[i_before] = sub(forces[i_before], force);  // sub because pointing other way
        before = after;
    }

    // apply force
    for (int i_rim = 0; i_rim < wheel->n_holes; ++i_rim) {
        acc_f += length(forces[i_rim]);
//        xy displacement = polar_scale(forces[i_rim], wheel->rim[i_rim],
//                damping * wheel->l_spoke[wheel->rim_to_hub[i_rim]] / wheel->e_spoke,
//                damping * wheel->l_chord / wheel->e_rim,
//                &acc_r, &acc_t);
        xy displacement = scalar_mult(damping, forces[i_rim]);
        wheel->rim[i_rim] = add(wheel->rim[i_rim], displacement);
    }
    ludebug(dbg, "Damping %5.3f; total r: %7.2emm, t: %7.2emm, f: %7.2eN", damping, acc_r, acc_t, acc_f);

LU_CLEANUP
    free(forces);
    LU_RETURN
}

int true(wheel *wheel) {
    LU_STATUS
    int n = 1000;
    for (int i = 0; i < n; ++i) {
//        LU_CHECK(relax(wheel, pow(10, -1 + (i-n)/(n/3.0))))
        LU_CHECK(relax(wheel, 0.1))
    }
    LU_NO_CLEANUP
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
    LU_CHECK(lace(wheel));
    LU_CHECK(true(wheel));
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
