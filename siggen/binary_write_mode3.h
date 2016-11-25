#ifndef _BINARY_WEITE_MODE3_H
#define _BINARY_WEITE_MODE3_H

#include <stdarg.h>

#define BASIS "basis_dat_copies/basis_0_45.dat"

typedef struct {
  int trclen;					 /* length of a single event */
  short iseg, r, p, z2;                          /* hit segment & integer cylindrical coordinates of grid point */
  float x, y, z;                                 /* cartesian coordinates of grid point */
  float t_drift;                                 /* PC drift time */
  float signal[NUM_SIGS][TIME_STEPS_C];          /* actual basis signals */
} Basis_Struct;

#endif
