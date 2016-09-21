
#ifndef SPOKES_LIB_H
#define SPOKES_LIB_H

#include "lu/log.h"

int unpack(lulog *dbg, const char *pattern, int **offsets, int *length, char *type, int *padding);
int dump_pattern(lulog *dbg, int *offsets, int length);
int rim_size(lulog *dbg, int length, int *holes);
int make_path(lulog *dbg, const char *pattern, char **path);

#endif
