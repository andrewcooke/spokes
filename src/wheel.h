
#ifndef SPOKES_WHEEL_H
#define SPOKES_WHEEL_H

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

int make_wheel(lulog *dbg, int *offsets, int length, int holes, int padding, char type, wheel **wheel);
void free_wheel(wheel *wheel);
int copy_wheel(lulog *dbg, wheel *w, wheel **c);

xy add(xy a, xy b);
xy sub(xy a, xy b);
double dot(xy a, xy b);
double length(xy xy);
xy scalar_mult(double k, xy a);
xy norm(xy a);
xy xy_on_circle(double r, int hole, int n_holes);

void open_plot(wheel *wheel, int nx, int ny, double line_width, cairo_t **cr, cairo_surface_t **surface);
void close_plot(cairo_t *cr, cairo_surface_t *surface, const char *path);

void plot_wheel(wheel *wheel, const char *path);
void plot_deform(wheel *original, wheel *deformed, load *l, const char *path, double scale);
int plot_multi_deform(lulog *dbg, wheel *original, wheel *deformed, load *l, const char *pattern);

#endif
