#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "pdecomp.h"

int read_basis(char *basis_file_name, PDecomp *pdd)
{
/* routine to read decomposition basis signal files from
   .dat unformatted binary file BASIS_FILE (defined in pdecomp.h)
   This file would normally have been created by
   the program convert_basis_sig2dat

   returns 0 on success, 1 on failure
   modifies:       Basis_Point pdd.basis[GRID_PTS];
                   int         pdd.grid_pos_lu[SRAD][SPHI][SZZZ]];
		   int         pdd.seg_basis_num[GRID_PTS][NUM_SIGS];
      defined in pdecomp.h.

   Author:  D.C. Radford    Aug 2004
*/
  char    header[256], test[256];
  int     i, ii, j, k, t, n, npts[NUM_SIGS];
  FILE   *file;
  static int no_bad_segs[1] = {-1};


  pdd->bad_segs = no_bad_segs;

  /* malloc the space for the basis etc */
  if (!(pdd->basis = (Basis_Point*) malloc(sizeof(Basis_Point) * GRID_PTS))) {
    printf("\nERROR  -  cannot malloc basis!\n\n");
    exit(-1);
  }

  printf("Reading basis signals from %s\n", basis_file_name);
  if (!(file=fopen(basis_file_name, "r"))) {
    printf("ERROR -- Cannot open file %s for reading.\n", basis_file_name);
    return 1;
  }

  fread(header, sizeof(header), 1, file);
  snprintf(test, sizeof(test),
	  "Segmented Point Contact basis signals; Cylindrical; version 1.0\n"
	  "%d basis point, %d grid segments, %d time steps\n"
	  "SRAD SPHI SZZZ: %d %d %d\n",
	  GRID_PTS, NUM_SIGS, TIME_STEPS_C, SRAD, SPHI, SZZZ);

  /* test string should match the header in the .dat file */
  /* if it does, read the basis data and grid position look-up table */

  if (!strstr(header, test)) {
    printf("Something's wrong with the basis data file!\nHeader:\n%s", header);
    fclose(file);
    return 1;
  }

  if ((i=fread(pdd->basis, sizeof(Basis_Point), GRID_PTS, file)) != GRID_PTS) {
    /* something's wrong */
    printf("Something's wrong with the basis data file!\n");
    fclose(file);
    printf("%d %d\n", i, GRID_PTS);
    return 1;
  }

  fclose(file);

  /* populate grid_pos_lu[SRAD][SPHI][SZZZ] */
  for (i=0; i<SRAD; i++) {
    for (j=0; j<SPHI; j++) {
      for (k=0; k<SZZZ; k++) {
	pdd->grid_pos_lu[i][j][k] = -1;
      }
    }
  }
  for (i=0; i<NUM_SIGS; i++) npts[i]=1;
  k = 0;
  for (i=0; i<GRID_PTS; i++) {
    pdd->grid_pos_lu[pdd->basis[i].ir][pdd->basis[i].ip][pdd->basis[i].iz] = i;
    /* also find max value of t20-t80 (since t80-t20 can be < 0) */
    if (k < (pdd->basis[i].t20 - pdd->basis[i].t80)) k = pdd->basis[i].t20 - pdd->basis[i].t80;
  }

  /* populate t_order[GRID_PTS] and seg_basis_num[GRID_PTS][NUM_SIGS] */
  n = 0;
  for (t=0; n<GRID_PTS; t++) {
    for (i=0; i<GRID_PTS; i++) {
      if (pdd->basis[i].t80 - pdd->basis[i].t20 + k == t) {
	pdd->t_order[n++] = i;
	pdd->seg_basis_num[npts[pdd->basis[i].iseg]++][pdd->basis[i].iseg] = i;
      }
    }
  }
  for (i=0; i<NUM_SIGS; i++) pdd->seg_basis_num[0][i] = npts[i];
  if (!pdd->quiet) printf(" >> k, t, n = %d %d %d\n", k, t, n); fflush(stdout);

  /* set bad segment signals to zero */
  for (i=0; pdd->bad_segs[i]>=0; i++) {
    ii = pdd->bad_segs[i];
    for (j=0; j<GRID_PTS; j++) {
      for (t=0; t<TIME_STEPS_C; t++) pdd->basis[j].signal[ii][t] = 0.0f;
      pdd->basis[j].lo_time[ii] = TIME_STEPS_C;
      pdd->basis[j].hi_time[ii] = 0;
    }
  }

  return 0;
}
