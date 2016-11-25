#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "gdecomp.h"

Basis_Point   *basis;                 /* basis-signal data */
int           grid_pos_lu[SSEG][SRAD][SPHI][SZZZ];    /* basis-grid position look-up table */
int           maxir[SSEG], maxip[SSEG], maxiz[SSEG];  /* max. values of ir, ip, iz for each segment */
Adaptive_Grid ags1[SSEG];             /* Adaptive Grid Search coarse grid, 1 segment */
int           *bad_segs;              /* list of segment numbers that should be disabled */

int read_basis(char *basis_file_name)
{
/* routine to read decomposition basis signal files from 
   .dat unformatted binary file BASIS_FILE (defined in gdecomp.h)
   This file would normally have been created by
   the program convert_basis_sig2dat

   returns 0 on success, 1 on failure
   modifies:       Basis_Point basis[GRID_PTS];
                   int         grid_pos_lu[SX][SY][SZ];
      defined in gdecomp.h.

   Author:  D.C. Radford    Aug 2004
*/
  char  header[256], test[256];
  int   i, ii, j, t, s;
  FILE  *file;
  static int no_bad_segs[1] = {-1};

  bad_segs = no_bad_segs;
  // if (set_bad_segs()) bad_segs = no_bad_segs;
  
  /* malloc the space for ags_dc[seg1][seg2] */
  if (!(basis = (Basis_Point*) malloc(sizeof(Basis_Point) * GRID_PTS))) {
    printf("\nERROR  -  cannot malloc basis!\n\n");
    exit(-1);
  }

  printf("Reading basis signals from %s\n", BASIS_FILE);
  if (!(file=fopen(basis_file_name, "r"))) {
    printf("ERROR -- Cannot open file %s for reading.\n", basis_file_name);
    return 1;
  }

  fread(header, 256, 1, file);
  snprintf(test, 256,
	  "GRETA basis signals; Cylindrical; version 1.2\n"
	  "%d basis point, %d grid segments, %d time steps\n"
	  "grid_pos_lu_rpz SRAD SPHI SZ: %d %d %d\n",
	  GRID_PTS, GRID_SEGS, TIME_STEPS, SRAD, SPHI, SZZZ);

  /* test string should match the header in the .dat file */
  /* if it does, read the basis data and grid position look-up table */

  if (!strstr(header, test) ||
      fread(basis, sizeof(Basis_Point), GRID_PTS, file) != GRID_PTS ||
      fread(grid_pos_lu, SSEG*SRAD*SPHI*SZZZ*(sizeof(int)), 1, file) != 1) {
    /* something's wrong */
    printf("Something's wrong with the basis data file!\n");
    fclose(file);
    return 1;
  }

  fclose(file);

  /* set bad segment signals to zero */
  for (i=0; bad_segs[i]>=0; i++) {
    ii = bad_segs[i];
    for (j=0; j<GRID_PTS; j++) {
      for (t=0; t<TIME_STEPS; t++) basis[j].signal[ii][t] = 0.0f;
      basis[j].lo_time[ii] = TIME_STEPS;
      basis[j].hi_time[ii] = 0;
    }
  }

  /* find the maximum values of ir, ip and iz for each segment */
  for (i=0; i<SSEG; i++) {
    maxir[i] = maxip[i] = maxiz[i] = 0;
  }
  for (i=0; i<GRID_PTS; i++) {
    s = basis[i].iseg;
    if (maxir[s] < basis[i].ir) maxir[s] = basis[i].ir;
    if (maxip[s] < basis[i].ip) maxip[s] = basis[i].ip;
    if (maxiz[s] < basis[i].iz) maxiz[s] = basis[i].iz;
  }
  if (!quiet) {
    printf(" seg  max.ir  max.ip  max.iz\n");
    for (i=0; i<SSEG; i++) {
      printf(" %3d %7d %7d %7d\n", i, maxir[i], maxip[i], maxiz[i]);
    }
    printf(" All %7d %7d %7d\n\n", SRAD-1, SPHI-1, SZZZ-1);
  }

  return 0;
}
