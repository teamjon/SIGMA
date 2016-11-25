#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <ctype.h>

#include "pdecomp.h"
#include "binary_write_mode3.h"

int main(int argc, char **argv) {
  PDecomp pdd;
  int  i,j,k;
  char header[256];											// Global header init
  FILE *fp;

  Basis_Struct bp;

  fp = fopen("mode3.dat","w+");
  snprintf(header, sizeof(header),									// Global header declare
          "Segmented Point Contact basis signals; Cylindrical; version 1.0\n"
          "%d basis point, %d grid segments, %d time steps\n"
          "SRAD SPHI SZZZ: %d %d %d\n",
          GRID_PTS, NUM_SIGS, TIME_STEPS_C, SRAD, SPHI, SZZZ);
  fwrite(header, sizeof(header), 1, fp);								// Global header write

  /* The struct pdd.basis contains
	--> iseg,ir,ip,iz;
	--> x,y,z;
	--> t80,t20;
	--> t_hi,t_lo;
	--> signals;
	-->t_drift;
  */

  if (!read_basis(BASIS, &pdd)) printf("Success!\n");							// Read basis file
  for (i=1; i<GRID_PTS; i++) {
    bp.trclen  = sizeof(bp);
    bp.iseg    = pdd.basis[i].iseg;									// Store variables in bp struct
    bp.r       = pdd.basis[i].ir;
    bp.p       = pdd.basis[i].ip;
    bp.z2      = pdd.basis[i].iz;
    bp.x       = pdd.basis[i].x;
    bp.y       = pdd.basis[i].y;
    bp.z       = pdd.basis[i].z;
    bp.t_drift = pdd.basis[i].t_drift;
    for(j=0; j<NUM_SIGS; j++) {
      for(k=0; k<TIME_STEPS_C; k++) {
        bp.signal[j][k] = pdd.basis[i].signal[j][k];
      }
    }
    fwrite(&bp,sizeof(bp),1,fp);									// Write bp struct to file
  }
  fclose(fp);
  return 0;
}
