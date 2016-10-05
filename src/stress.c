
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "cairo/cairo.h"
#include "cblas.h"
#include "lapacke_utils.h"

#include "lu/status.h"
#include "lu/log.h"
#include "lu/files.h"
#include "lu/dynamic_memory.h"

#include "lib.h"

lulog *dbg = NULL;

#define LA_CHECK(log, stmt) if ((status = stmt)) {luerror(log, "LAPACKE error %d at %s:%d", status, __FILE__, __LINE__); goto exit;}

#define X 0
#define Y 1


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

int calculate_force(int n, double *a, double *x, double *b, double **result) {
    LU_STATUS
    LU_ALLOC(dbg, *result, 2 * n)
    memcpy(*result, b, 2 * n * sizeof(*b));
    cblas_dgemv(CblasColMajor, CblasNoTrans, 2*n, 2*n, 1, a, 2*n, x, 1, -1, *result, 1);
    LU_NO_CLEANUP
}

int draw_forces_initial(cairo_t *cr, wheel *wheel, double *a, double *b) {

    LU_STATUS

    int n = wheel->n_holes;

    double *dxy = NULL, *forces = NULL;
    LU_ALLOC(dbg, dxy, 2 * n)
    LU_CHECK(calculate_force(n, a, dxy, b, &forces))
    cairo_set_source_rgb(cr, 0, 1, 1);
    double mag = 0.1;
    for (int i = 0; i < n; ++i) {
        draw_line(cr, wheel->rim[i].x, wheel->rim[i].y, wheel->rim[i].x+mag*forces[2*i+X], wheel->rim[i].y+mag*forces[2*i+Y]);
    }

LU_CLEANUP
    free(dxy);
    free(forces);
    LU_RETURN
}

int plot_forces_initial(wheel *wheel, double *a, double *b, const char *path) {
    LU_STATUS
    cairo_surface_t *surface = NULL;
    cairo_t *cr = NULL;
    open_plot(wheel, 500, 500, wheel->r_rim / 100, &cr, &surface);
    draw_wheel(cr, wheel);
    LU_CHECK(draw_forces_initial(cr, wheel, a, b))
    close_plot(cr, surface, path);
    LU_NO_CLEANUP
}

int draw_forces_for_radial_movement(cairo_t *cr, wheel *wheel, double *a, double *b) {

    LU_STATUS

    int n = wheel->n_holes;

    double *dxy = NULL, *forces = NULL;
    LU_ALLOC(dbg, dxy, 2 * n)

    double extn = 1e-3;
    for (int i = 0; i < n; ++i) {
        dxy[2*i+X] = -extn * wheel->rim[i].x;
        dxy[2*i+Y] = -extn * wheel->rim[i].y;
    }
    LU_CHECK(calculate_force(n, a, dxy, b, &forces))
    cairo_set_source_rgb(cr, 1, 0, 0);
    double mag = 1e-3;
    for (int i = 0; i < n; ++i) {
        draw_line(cr, wheel->rim[i].x, wheel->rim[i].y, wheel->rim[i].x+mag*forces[2*i+X], wheel->rim[i].y+mag*forces[2*i+Y]);
    }

LU_CLEANUP
    free(dxy);
    free(forces);
    LU_RETURN
}

int plot_forces_for_radial_movement(wheel *wheel, double *a, double *b, const char *path) {
    LU_STATUS
    cairo_surface_t *surface = NULL;
    cairo_t *cr = NULL;
    open_plot(wheel, 500, 500, wheel->r_rim / 100, &cr, &surface);
    draw_wheel(cr, wheel);
    LU_CHECK(draw_forces_for_radial_movement(cr, wheel, a, b))
    close_plot(cr, surface, path);
    LU_NO_CLEANUP
}

void rotate(double dtheta, double x1, double y1, double *x2, double *y2) {
    double r = sqrt(x1*x1+ y1*y1), theta = atan2(y1, x1);
    theta += dtheta;
    *x2 = r * cos(theta);
    *y2 = r * sin(theta);
}

int draw_forces_for_tangential_movement(cairo_t *cr, wheel *wheel, double *a, double *b, double dtheta) {

    LU_STATUS

    int n = wheel->n_holes;

    double *dxy = NULL, *forces = NULL;
    LU_ALLOC(dbg, dxy, 2 * n)

    double x1, y1, x2, y2;

    x1 = wheel->rim[0].x;
    y1 = wheel->rim[0].y;
    rotate(1e5 * dtheta, x1, y1, &x2, &y2);
    draw_line(cr, x1, y1, x2, y2);

    double mag = 1e-1;
    for (int i = 0; i < n; ++i) {

        for (int j = 0; j < 2*n; ++j) dxy[j] = 0;

        x1 = wheel->rim[i].x;
        y1 = wheel->rim[i].y;
        rotate(dtheta, x1, y1, &x2, &y2);
        dxy[2*i+X] = x2 - x1;
        dxy[2*i+Y] = y2 - y1;

        LU_CHECK(calculate_force(n, a, dxy, b, &forces))
        draw_line(cr, wheel->rim[i].x, wheel->rim[i].y, wheel->rim[i].x+mag*forces[2*i+X], wheel->rim[i].y+mag*forces[2*i+Y]);
        free(forces); forces = NULL;
    }

LU_CLEANUP
    free(dxy);
    free(forces);
    LU_RETURN
}

int plot_forces_for_tangential_movement(wheel *wheel, double *a, double *b, const char *path) {
    LU_STATUS
    cairo_surface_t *surface = NULL;
    cairo_t *cr = NULL;
    open_plot(wheel, 500, 500, wheel->r_rim / 100, &cr, &surface);
    draw_wheel(cr, wheel);
    cairo_set_source_rgb(cr, 1, 0, 0);
    LU_CHECK(draw_forces_for_tangential_movement(cr, wheel, a, b, 1e-6))
    cairo_set_source_rgb(cr, 0, 1, 0);
    LU_CHECK(draw_forces_for_tangential_movement(cr, wheel, a, b, -1e-6))
    close_plot(cr, surface, path);
    LU_NO_CLEANUP
}

void draw_displacements(cairo_t *cr, wheel *wheel, double *b) {
    cairo_set_source_rgb(cr, 0, 1, 0);
    int n = wheel->n_holes;
//    double mag = 1e3;
    double mag = 1;
    for (int i = 0; i < n; ++i) {
        draw_line(cr, wheel->rim[i].x, wheel->rim[i].y, wheel->rim[i].x+mag*b[2*i+X], wheel->rim[i].y+mag*b[2*i+Y]);
    }
}

int plot_displacements(wheel *wheel, double *b, const char *path) {
    LU_STATUS
    cairo_surface_t *surface = NULL;
    cairo_t *cr = NULL;
    open_plot(wheel, 500, 500, wheel->r_rim / 100, &cr, &surface);
    draw_wheel(cr, wheel);
    draw_displacements(cr, wheel, b);
    close_plot(cr, surface, path);
    LU_NO_CLEANUP
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

void inc_a(double *a, int n, int xy_force, int i_force, int xy_posn, int i_posn, double value) {
    a[2*n*(2*i_posn+xy_posn) + 2*i_force+xy_force] += value;
}

void inc_b(double *b, int n, int xy_force, int i_force, double value) {
    b[2*i_force+xy_force] += value;
}

void inc_spoke(wheel *wheel, double *a, double *b, int n, int i_xy2,
        double t, double dtdx2, double dtdy2, double cos_theta, double sin_theta) {

    // x and y components of t.  force is ax-b so need to negate b, but
    // this is a tension, so is already "-ve".
    // we just need to take the components.
    inc_b(b, n, X, i_xy2,            t * cos_theta);
    inc_b(b, n, Y, i_xy2,            t * sin_theta);

    // the x and y components (left column) of the derivatives wrt x and y (right)
    // (in other words, left is where force is felt, right is who produces it)
    // this is tension so we need to negate for force.
    inc_a(a, n, X, i_xy2, X, i_xy2, -dtdx2 * cos_theta);
    inc_a(a, n, X, i_xy2, Y, i_xy2, -dtdy2 * cos_theta);
    inc_a(a, n, Y, i_xy2, X, i_xy2, -dtdx2 * sin_theta);
    inc_a(a, n, Y, i_xy2, Y, i_xy2, -dtdy2 * sin_theta);

    if (i_xy2 == 0) {
        int i_hub = wheel->rim_to_hub[i_xy2];
        xy rim = wheel->rim[i_xy2], hub = wheel->hub[i_hub];
        double l0 = wheel->l_spoke[i_hub], l = length(sub(rim, hub));
        ludebug(dbg, "Tension at hole %d (%.0f,%.0f) from spoke extended %4.2f%% with modulus %.0f",
                i_xy2, wheel->rim[i_xy2].x, wheel->rim[i_xy2].y, 100 * (l - l0) / l0, wheel->e_spoke);
        ludebug(dbg, "is %0.f with derivative (along spoke) of %.0f dx %+.0f dy (cos %5.3f, sin %5.3f)",
                t, dtdx2, dtdy2, cos_theta, sin_theta);
    }
}

void inc_chord(wheel *wheel, double *a, double *b, int n, int i_xy1, int i_xy2,
        double t, double dtdx2, double dtdy2, double cos_theta, double sin_theta) {

    // as for spokes, the t value is "already -ve" at xy2, but for xy1
    // we must reverse signs
    inc_b(b, n, X, i_xy2,            t * cos_theta);
    inc_b(b, n, Y, i_xy2,            t * sin_theta);
    inc_b(b, n, X, i_xy1,           -t * cos_theta);
    inc_b(b, n, Y, i_xy1,           -t * sin_theta);

    // again, for "xy2 xy2" we are the same as spokes
    inc_a(a, n, X, i_xy2, X, i_xy2, -dtdx2 * cos_theta);
    inc_a(a, n, X, i_xy2, Y, i_xy2, -dtdy2 * cos_theta);
    inc_a(a, n, Y, i_xy2, X, i_xy2, -dtdx2 * sin_theta);
    inc_a(a, n, Y, i_xy2, Y, i_xy2, -dtdy2 * sin_theta);

    // but for other combinations, signs swap accordingly
    inc_a(a, n, X, i_xy2, X, i_xy1,  dtdx2 * cos_theta);
    inc_a(a, n, X, i_xy2, Y, i_xy1,  dtdy2 * cos_theta);
    inc_a(a, n, Y, i_xy2, X, i_xy1,  dtdx2 * sin_theta);
    inc_a(a, n, Y, i_xy2, Y, i_xy1,  dtdy2 * sin_theta);

    inc_a(a, n, X, i_xy1, X, i_xy2,  dtdx2 * cos_theta);
    inc_a(a, n, X, i_xy1, Y, i_xy2,  dtdy2 * cos_theta);
    inc_a(a, n, Y, i_xy1, X, i_xy2,  dtdx2 * sin_theta);
    inc_a(a, n, Y, i_xy1, Y, i_xy2,  dtdy2 * sin_theta);

    inc_a(a, n, X, i_xy1, X, i_xy1, -dtdx2 * cos_theta);
    inc_a(a, n, X, i_xy1, Y, i_xy1, -dtdy2 * cos_theta);
    inc_a(a, n, Y, i_xy1, X, i_xy1, -dtdx2 * sin_theta);
    inc_a(a, n, Y, i_xy1, Y, i_xy1, -dtdy2 * sin_theta);

    if (i_xy2 == 0) {
        xy xy2 = wheel->rim[i_xy2], xy1 = wheel->rim[i_xy1];
        double l0 = wheel->l_chord, l = length(sub(xy2, xy1));
        ludebug(dbg, "Tension at hole %d (%.0f,%.0f) from rim chord at (%0.f,%0.f) extended %4.2f%% with modulus %.0f",
                i_xy2, wheel->rim[i_xy2].x, wheel->rim[i_xy2].y, wheel->rim[i_xy1].x, wheel->rim[i_xy1].y, 100 * (l - l0) / l0, wheel->e_rim);
        ludebug(dbg, "is %0.f with derivative (along chord) of %.0f dx %+.0f dy (cos %5.3f, sin %5.3f)",
                t, dtdx2, dtdy2, cos_theta, sin_theta);
    }
}

int validate(int n, double *a, double *x, double *b) {

    LU_STATUS

    double *forces = NULL;
    LU_CHECK(calculate_force(n, a, x, b, &forces))
    for (int i = 0; i < n; ++i) {
        luinfo(dbg, "%d: (%g, %g) (%g, %g)", i, forces[2*i+X], forces[2*i+Y], x[2*i+X], x[2*i+Y]);
    }

LU_CLEANUP
    free(forces);
    LU_RETURN
}

void calc_tension(xy xy1, xy xy2, double l0, double e,
        double *t, double *dtdx2, double *dtdy2, double *cos_theta, double *sin_theta) {

    // a line from x1,y1 to x2,y2 with length l makes an angle theta to
    // the horizontal where theta = atan2((y2-y1)/(x2-x1)),
    // cos(theta) = (x2-x1)/l and sin(theta) = (y2-y1)/l

    double delta_x = xy2.x - xy1.x, delta_y = xy2.y - xy1.y;
    double l = sqrt(delta_x * delta_x + delta_y * delta_y);
    *cos_theta = delta_x / l; *sin_theta = delta_y / l;

    // if xy2 moves by dx2, dy2 then the new length is approx (1st order)
    // l + dx2 * cos(theta) + dy2 * sin(theta)
    // this isn't obvious to me, but you can show it from the taylor expansion
    // of l^2 = (delta_x+dx)^2 + (delta_y+dy)^2
    // which gives
    // l' = l + delta_x * (dx2-dx1) / l + delta_y * (dy2-dy1) / l

    double dldx2 = delta_x / l, dldy2 = delta_y / l;

    // now tension is E(l-l0)/l0 where l0 is the original length
    // (tension is +ve when extended, ie when l > l0)

    *t = e * (l - l0) / l0;

    // the derivative of tension wrt dx and dy are then

    *dtdx2 = e * dldx2 / l0; *dtdy2 = e * dldy2 / l0;
}

int true(wheel *wheel, double damping) {

    LU_STATUS

    double *a = NULL, *b = NULL, *a_copy = NULL, *b_copy = NULL;  // we're going to solve ax = b
    int n = wheel->n_holes;
    lapack_int *pivot;
    double delta = 0.0;
    double t, dtdx2, dtdy2, cos_theta, sin_theta;

    ludebug(dbg, "e  spoke %g, rim %g", wheel->e_spoke, wheel->e_rim);
    ludebug(dbg, "chord length %g", wheel->l_chord);

    ludebug(dbg, "hub");
    for (int i = 0; i < n; ++i) {
        printf("%d  %g, %g\n", i, wheel->hub[i].x, wheel->hub[i].y);
    }

    ludebug(dbg, "rim");
    for (int i = 0; i < n; ++i) {
        printf("%d  %g, %g\n", i, wheel->rim[i].x, wheel->rim[i].y);
    }

    ludebug(dbg, "spokes");
    for (int i = 0; i < n; ++i) {
        double l = length(sub(wheel->hub[i], wheel->rim[wheel->hub_to_rim[i]]));
        printf("%d  %g, %g  (diff %g, force %g)\n", i, l, wheel->l_spoke[i], l - wheel->l_spoke[i],
                wheel->e_spoke * (l - wheel->l_spoke[i]) / wheel->l_spoke[i]);
    }

    ludebug(dbg, "chords");
    for (int i = 0; i < n; ++i) {
        int j = (i - 1 + wheel->n_holes) % wheel->n_holes;
        double l = length(sub(wheel->rim[i], wheel->rim[j]));
        printf("%d-%d  %g, %g  (diff %g, force %g)\n", j, i,
                l, wheel->l_chord, wheel->l_chord - l,
                wheel->e_rim * (wheel->l_chord - l) / wheel->l_chord);
    }

    // we interleave x,y and the arrays are column major so for hole i
    // x[2i+X] = dx, x[2i+Y] = dy
    // and the force on hole i in the x direction is the sum over forces from j
    // sum(j)(a[2n*j+2i] * x[j]) - b[2i]

    LU_ALLOC(dbg, a, 4 * n * n)
    LU_ALLOC(dbg, b, 2 * n)
    LU_ALLOC(dbg, pivot, 2 * n);

    for (int i_rim = 0; i_rim < n; ++i_rim) {
        // for spokes xy1 is fixed, so we only need the tension etc at xy2
        int i_hub = wheel->rim_to_hub[i_rim];
        calc_tension(wheel->hub[i_hub], wheel->rim[i_rim], wheel->l_spoke[i_hub], wheel->e_spoke,
                &t, &dtdx2, &dtdy2, &cos_theta, &sin_theta);
        inc_spoke(wheel, a, b, n, i_rim, t, dtdx2, dtdy2, cos_theta, sin_theta);
    }

    for (int i_xy2 = 0; i_xy2 < n; ++i_xy2) {
        // for rims the chord from xy1 to xy2 affects two holes
        int i_xy1 = (i_xy2-1+n) % n;
        calc_tension(wheel->rim[i_xy1], wheel->rim[i_xy2], wheel->l_chord, wheel->e_rim,
                &t, &dtdx2, &dtdy2, &cos_theta, &sin_theta);
        inc_chord(wheel, a, b, n, i_xy1, i_xy2, t, dtdx2, dtdy2, cos_theta, sin_theta);
    }

    ludebug(dbg, "A");
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

    ludebug(dbg, "b");
    for (int i = 0; i < n; ++i) {
        printf("%d  %g, %g\n", i, b[2*i+X], b[2*i+Y]);
    }

    double b_mag = 0;
    for (int i = 0; i < 2 * n; ++i) b_mag += fabs(b[i]);
    luinfo(dbg, "Total static force is %g", b_mag);

    LU_ALLOC(dbg, a_copy, 4 * n * n);
    memcpy(a_copy, a, 4 * n * n * sizeof(*a));
    LU_ALLOC(dbg, b_copy, 2 * n);
    memcpy(b_copy, b, 2 * n * sizeof(*b));

    LU_CHECK(plot_forces_initial(wheel, a_copy, b_copy, "forces-initial.png"))
    LU_CHECK(plot_forces_for_radial_movement(wheel, a_copy, b_copy, "forces-radial.png"))
    LU_CHECK(plot_forces_for_tangential_movement(wheel, a_copy, b_copy, "forces-tangential.png"))

    LA_CHECK(dbg, LAPACKE_dgetrf(LAPACK_COL_MAJOR, 2*n, 2*n, a, 2*n, pivot))
    LA_CHECK(dbg, LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', 2*n, 1, a, 2*n, pivot, b, 2*n))

//    double ferr, berr;
//    LA_CHECK(dbg, LAPACKE_dgerfs(LAPACK_COL_MAJOR, 'N', 2*n, 1, a_copy, 2*n, a, 2*n, pivot, b_copy, 2*n, b, 2*n,
//            &ferr, &berr));
//    luinfo(dbg, "Errors %g %g", ferr, berr);

    ludebug(dbg, "x");
    for (int i = 0; i < n; ++i) {
        printf("%d  %g, %g\n", i, b[2*i+X], b[2*i+Y]);
    }

    LU_CHECK(plot_displacements(wheel, b, "displacements.png"))

    for (int i = 0; i < n; ++i) {
        delta += damping * sqrt(b[2*i+X]*b[2*i+X] + b[2*i+Y]*b[2*i+Y]);
        wheel->rim[i].x += damping * b[2*i+X];
        wheel->rim[i].y += damping * b[2*i+Y];
    }
    ludebug(dbg, "Total shift %7.2gmm", delta);

    LU_CHECK(validate(n, a_copy, b, b_copy))


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
//    (*wheel)->l_chord *= 1.0001;
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
    holes = 4;
    LU_CHECK(make_wheel(offsets, length, holes, padding, type, &wheel));
    LU_CHECK(lace(wheel));
//    for (int i = 0; i < 100; ++i) LU_CHECK(true(wheel, 0.01));
    LU_CHECK(true(wheel, 1));
    LU_CHECK(true(wheel, 1));
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
