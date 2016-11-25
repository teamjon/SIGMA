#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <ctype.h>
#include <time.h>

#include "pdecomp.h"
#include "geant_siggen.h"

int main(int argc, char **argv) {
  Mode2_Struct m2;
  int  i,j,k;
  char header[256];                                                                                     // Global header init
  FILE *fp;
  srand(time(NULL));

  /* malloc the space for the basis etc */
/*  if (!(m2.basis = (Mode2_Struct*) malloc(sizeof(Mode2_Struct) * NEVNTS))) {
    printf("\nERROR  -  cannot malloc basis!\n\n");
    exit(-1);
  }
*/
  fp = fopen("mode2.dat","w+");
  snprintf(header, sizeof(header),                                                                      // Global header declare
          "Segmented Point Contact basis signals; Cylindrical; version 1.0\n"
          "%d basis point, %d grid segments, %d time steps\n"
          "SRAD SPHI SZZZ: %d %d %d\n",
          GRID_PTS, NUM_SIGS, TIME_STEPS_C, SRAD, SPHI, SZZZ);
  fwrite(header, sizeof(header), 1, fp);

  for(j=0; j<NEVNTS; j++) {
  m2.type		= 20;
  m2.crystal_id	= 1;
  m2.num		= 1;
  m2.tot_e		= 1;
  m2.timestamp	= 0;
  m2.t0		= 0;
  m2.chisq		= 0;
  m2.norm_chisq	= 0;
  for(i=0; i<MAX_INTPTS; i++) {
    if(i<m2.num) {
      m2.intpts[i].x        = (rand() % 700) / 10;
      m2.intpts[i].y        = (rand() % 700) / 10;
      m2.intpts[i].z        = (rand() % 800) / 10;
      m2.intpts[i].e        = 1;
      m2.intpts[i].seg      = 0;
      m2.intpts[i].seg_ener = 0;
    }
    if(i>=m2.num) {
      m2.intpts[i].x        = -1;
      m2.intpts[i].y        = -1;
      m2.intpts[i].z        = -1;
      m2.intpts[i].e        = -1;
      m2.intpts[i].seg      = -1;
      m2.intpts[i].seg_ener = -1;
    }
  }

  fwrite(&m2, sizeof(m2), 1, fp);
}
  fclose(fp);
  return 0;
}
