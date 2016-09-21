
#ifndef SPOKES_LIB_H
#define SPOKES_LIB_H

#include "lu/log.h"

void draw_circle(cairo_t *cr, float r);
void draw_line(cairo_t *cr, float x0, float y0, float x1, float y1);
int unpack(lulog *dbg, const char *pattern, int **offsets, int *length, char *type, int *padding);
int dump_pattern(lulog *dbg, int *offsets, int length);
int rim_size(lulog *dbg, int length, int *holes);
int make_path(lulog *dbg, const char *pattern, char **path);

#endif
