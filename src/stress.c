
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "cairo/cairo.h"
#include "lapacke_utils.h"

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"

#include "lib.h"

lulog *dbg = NULL;

#define LA_CHECK(log, stmt) if ((status = stmt)) {luerror(log, "LAPACKE error %d at %s:%d", status, __FILE__, __LINE__); goto exit;}


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

#define X 0
#define Y 1

void inc_a(double *a, int n, int xy_force, int i_force, int xy_posn, int i_posn, double value) {
    a[2*n*(2*i_posn+xy_posn) + 2*i_force+xy_force] = value;
}

void inc_b(double *b, int n, int xy_force, int i_force, double value) {
    b[2*i_force+xy_force] = value;
}

void inc_spoke(wheel *wheel, double *a, double *b, int n, int i_rim,
        double fxx, double fxy, double fx, double fyx, double fyy, double fy) {

    inc_a(a, n, X, i_rim, X, i_rim, -fxx);
    inc_a(a, n, X, i_rim, Y, i_rim, -fxy);
    inc_b(b, n, X, i_rim,            fx);
    inc_a(a, n, Y, i_rim, X, i_rim, -fyx);
    inc_a(a, n, Y, i_rim, Y, i_rim, -fyy);
    inc_b(b, n, Y, i_rim,            fy);

    if (i_rim == 31) {
        int i_hub = wheel->rim_to_hub[i_rim];
        xy rim = wheel->rim[i_rim], hub = wheel->hub[i_hub];
        double l0 = wheel->l_spoke[i_hub], l = length(sub(rim, hub));
        ludebug(dbg, "Force at hole %d (%.0f,%.0f) from spoke extended %4.2f%% with modulus %.0f",
                i_rim, wheel->rim[i_rim].x, wheel->rim[i_rim].y, 100 * (l - l0) / l0, wheel->e_spoke);
        ludebug(dbg, "is (%.0f dx %+.0f dy %+.0f, %.0f dx %+.0f dy %+.0f)",
                -fxx, -fxy, -fx, -fyx, -fyy, -fy);
    }
}

void inc_chord(wheel *wheel, double *a, double *b, int n, int i_before, int i_after,
        double fxx, double fxy, double fx, double fyx, double fyy, double fy) {

    // before is x1, after is x2
    inc_a(a, n, X, i_after,  X, i_before,  fxx);
    inc_a(a, n, X, i_after,  Y, i_before,  fxy);
    inc_a(a, n, X, i_after,  X, i_after,  -fxx);
    inc_a(a, n, X, i_after,  Y, i_after,  -fxy);
    inc_b(b, n, X, i_after,               -fx);
    inc_a(a, n, Y, i_after,  X, i_before,  fyx);
    inc_a(a, n, Y, i_after,  Y, i_before,  fyy);
    inc_a(a, n, Y, i_after,  X, i_after,  -fyx);
    inc_a(a, n, Y, i_after,  Y, i_after,  -fyy);
    inc_b(b, n, Y, i_after,               -fy);

    if (i_after == 31) {
        xy after = wheel->rim[i_after], before = wheel->rim[i_before];
        double l0 = wheel->l_chord, l = length(sub(after, before));
        ludebug(dbg, "Force at hole %d xy2=(%.0f,%.0f) from chord to hole %d xy1=(%.0f,%.0f) compressed %4.2f%% with modulus %.0f",
                i_after, after.x, after.y, i_before, before.x, before.y, 100 * (l0 - l) / l0, wheel->e_rim);
        ludebug(dbg, "is (%.0f dx1 %+.0f dx2 %+.0f dy1 %+.0f dy2 %+.0f, %.0f dx1 %+.0f dx2 %+.0f dy1 %+.0f dy2 %+.0f)",
                fxx, -fxx, fxy, -fxy, -fx, fyx, -fyx, fyy, -fyy, -fy);
    }

    inc_a(a, n, X, i_before, X, i_before, -fxx);
    inc_a(a, n, X, i_before, Y, i_before, -fxy);
    inc_a(a, n, X, i_before, X, i_after,   fxx);
    inc_a(a, n, X, i_before, Y, i_after,   fxy);
    inc_b(b, n, X, i_before,              -fx);
    inc_a(a, n, Y, i_before, X, i_before, -fyx);
    inc_a(a, n, Y, i_before, Y, i_before, -fyy);
    inc_a(a, n, Y, i_before, X, i_after,   fyx);
    inc_a(a, n, Y, i_before, Y, i_after,   fyy);
    inc_b(b, n, Y, i_before,              -fy);

    if (i_before == 31) {
        xy after = wheel->rim[i_after], before = wheel->rim[i_before];
        double l0 = wheel->l_chord, l = length(sub(after, before));
        ludebug(dbg, "Force at hole %d xy1=(%.0f,%.0f) from chord to hole %d xy2=(%.0f,%.0f) compressed %4.2f%% with modulus %.0f",
                i_before, before.x, before.y, i_after, after.x, after.y, 100 * (l0 - l) / l0, wheel->e_rim);
        ludebug(dbg, "is (%.0f dx1 %+.0f dx2 %+.0f dy1 %+.0f dy2 %+.0f, %.0f dx1 %+.0f dx2 %+.0f dy1 %+.0f dy2 %+.0f)",
                -fxx, fxx, -fxy, fxy, -fx, -fyx, fyx, -fyy, fyy, -fy);
    }
}

int validate(int n, double *a, double *x, double *b) {

    LU_STATUS

    double *b_test = NULL;

    LU_ALLOC(dbg, b_test, n * n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            b_test[i] += a[n*j+i] * x[j];
        }
    }

    for (int i = 0; i < n; ++i) {
        luinfo(dbg, "%d: %g = %g (%g)", i, b_test[i], b[i], x[i]);
    }

LU_CLEANUP
    free(b_test);
    LU_RETURN
}

int true(wheel *wheel) {

    LU_STATUS

    double *a = NULL, *b = NULL, *a_copy = NULL, *b_copy = NULL;  // we're going to solve ax = b
    int n = wheel->n_holes;
    lapack_int *pivot;
    double delta = 0.0;

    // we interleave x,y and the arrays are column major so for hole i
    // x[2i] = dx, x[2i+1] = dy
    // and the force on hole i in the x direction is the sum over forces from j
    // sum(j)(a[2n*2j+2i] * x[2j]) - b[2i]

    LU_ALLOC(dbg, a, 4 * n * n)
    LU_ALLOC(dbg, b, 2 * n)
    LU_ALLOC(dbg, pivot, 2 * n);

    // the length from (x1+dx1,y1+dy1) to (x2+dy2,y2+dy2) is approx (small dx,dy)
    // l + (x2-x1)*(dx2-dx1)/l + (y2-y1)*(dy2-dy1)/l
    // where l is the length from (x1,y1) to (x2,y2)
    // you can get this by expanding the usual pythag expression as a taylor
    // expansion in dx,dy or by simply looking at the trignometry involved.

    // the magnitude of the force on an element whose at rest length is L:
    // F = -(l + (x2-x1)*(dx2-dx1)/l + (y2-y1)*(dy2-dy1)/l - L)/L * E
    // (negative because if l > L then it's extended and wants to contract)
    //   = -((x2-x1)*(dx2-dx1) + (y2-y1)*(dy2-dy1)) * E/(Ll) - (l-L)/L * E
    // in the x direction, at x2, we need "cos theta":
    // Fx = (x2-x1)/l * F
    //    = -((x2-x1)*(dx2-dx1) + (y2-y1)*(dy2-dy1)) * (x2-x1)E/(Ll^2) - (x2-x1)(1/L-1/l)E
    //    = ax - b

    // for spokes x1,y1 is fixed so dx1=dy1=0
    for (int i_rim = 0; i_rim < n; ++i_rim) {
        int i_hub = wheel->rim_to_hub[i_rim];

        double x1 = wheel->hub[i_hub].x, y1 = wheel->hub[i_hub].y;
        double x2 = wheel->rim[i_rim].x, y2 = wheel->rim[i_rim].y;
        double delta_x = x2 - x1, delta_y = y2 - y1;
        double l = sqrt(delta_x * delta_x + delta_y * delta_y), lsq = l * l;
        double l0 = wheel->l_spoke[i_hub];
        double e = wheel->e_spoke;

        double fxx = delta_x * delta_x * e / (l0 * lsq);
        double fxy = delta_x * delta_y * e / (l0 * lsq);
        double fx =  delta_x * (1/l0 - 1/l) * e;
        double fyx = fxy;
        double fyy = delta_y * delta_y * e / (l0 * lsq);
        double fy =  delta_y * (1/l0 - 1/l) * e;

        inc_spoke(wheel, a, b, n, i_rim, fxx, fxy, fx, fyx, fyy, fy);
    }

    // hub chords
    for (int i_after = 0; i_after < n; ++i_after) {
        int i_before = (i_after-1+n) % n;

        double x1 = wheel->rim[i_before].x, y1 = wheel->rim[i_before].y;
        double x2 = wheel->rim[i_after].x, y2 = wheel->rim[i_after].y;
        double delta_x = x2 - x1, delta_y = y2 - y1;
        double l = sqrt(delta_x * delta_x + delta_y * delta_y), lsq = l * l;
        double l0 = wheel->l_chord;
        double e = wheel->e_rim;

        double fxx = delta_x * delta_x * e / (l0 * lsq);
        double fxy = delta_x * delta_y * e / (l0 * lsq);
        double fx =  delta_x * (1/l0 - 1/l) * e;
        double fyx = fxy;
        double fyy = delta_y * delta_y * e / (l0 * lsq);
        double fy =  delta_y * (1/l0 - 1/l) * e;

        inc_chord(wheel, a, b, n, i_before, i_after, fxx, fxy, fx, fyx, fyy, fy);
    }

    for (int i = 0; i < 2 * n; ++i) {
        printf("%d  ", i);
        for (int j = 0; j < 2 * n; ++j) {
            double value = a[i+2*n*j];
            if (value) {
                printf("%d %g  ", j, value);
            }
        }
        printf("\n");
    }

    LU_ALLOC(dbg, a_copy, 4 * n * n);
    memcpy(a_copy, a, 4 * n * n * sizeof(*a));
    LU_ALLOC(dbg, b_copy, 2 * n);
    memcpy(b_copy, b, 2 * n * sizeof(*b));

    LA_CHECK(dbg, LAPACKE_dgetrf(LAPACK_COL_MAJOR, 2*n, 2*n, a, 2*n, pivot))
    LA_CHECK(dbg, LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', 2*n, 1, a, 2*n, pivot, b, 2*n))

    double ferr, berr;
    LA_CHECK(dbg, LAPACKE_dgerfs(LAPACK_COL_MAJOR, 'N', 2*n, 1, a_copy, 2*n, a, 2*n, pivot, b_copy, 2*n, b, 2*n,
            &ferr, &berr));
    luinfo(dbg, "Errors %g %g", ferr, berr);

    for (int i = 0; i < n; ++i) {
        delta += sqrt(b[2*i]*b[2*i] + b[2*i+1]*b[2*i+1]);
        wheel->rim[i].x += b[2*i];
        wheel->rim[i].y += b[2*i+1];
    }
    ludebug(dbg, "Total shift %7.2gmm", delta);

    LU_CHECK(validate(2 * n, a_copy, b, b_copy))

LU_CLEANUP
    free(a);
    free(b);
    free(pivot);
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
