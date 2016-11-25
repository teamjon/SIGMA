#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <ctype.h>
#include <time.h>

//#include "pdecomp.h"
#include "geant_siggen.h"

int main(int argc, char **argv) {
  M2Decomp *m2;
  int  i,j,k;
  char header[256];                                                                                     // Global header init
  FILE *fp;
  srand(time(NULL));
printf("\n\n1\n\n");
  /* malloc the space for the basis etc */
  if (!(m2->basis = (Mode2_Struct*) malloc(sizeof(Mode2_Struct) * NEVNTS))) {
   printf("\nERROR  -  cannot malloc basis!\n\n");
    exit(-1);
  }
printf("\n\n2\n\n");
  fp = fopen("mode2.dat","w+");
  snprintf(header, sizeof(header),                                                                      // Global header declare
          "Segmented Point Contact basis signals; Cylindrical; version 1.0\n"
          "%d basis point, %d grid segments, %d time steps\n"
          "SRAD SPHI SZZZ: %d %d %d\n",
          GRID_PTS, NUM_SIGS, TIME_STEPS_C, SRAD, SPHI, SZZZ);
  fwrite(header, sizeof(header), 1, fp);

  for(i=0; i<NEVNTS; i++) {
	printf("\n\n3\n\n");
    m2->basis[i]->type = 20;
  printf("\n\n4\n\n");
    m2->basis[i]->crystal_id = 1;
    m2->basis[i]->num = 1;
    m2->basis[i]->tot_e = 1;
    m2->basis[i]->timestamp = 0;
    m2->basis[i]->t0 = 0;
    m2->basis[i]->chisq = 0;
    m2->basis[i]->norm_chisq = 0;
    for(j=0; j<MAX_INTPTS; j++) {
      if(j< m2->basis[i]->num) {
        m2->basis[i]->intpts[j].x = (rand() % 700) / 10;
        m2->basis[i]->intpts[j].y = (rand() % 700) / 10;
        m2->basis[i]->intpts[j].z = (rand() % 800) / 10;
        m2->basis[i]->intpts[j].e = 1;
        m2->basis[i]->intpts[j].seg = 0;
        m2->basis[i]->intpts[j].seg_ener = 0;
      }
      if(j>= m2->basis[i]->num) {
        m2->basis[i]->intpts[j].x = -1;
        m2->basis[i]->intpts[j].y = -1;
        m2->basis[i]->intpts[j].z = -1;
        m2->basis[i]->intpts[j].e = -1;
        m2->basis[i]->intpts[j].seg = -1;
        m2->basis[i]->intpts[j].seg_ener = -1;
      }
    }
  }
  fwrite(&m2, sizeof(m2), 1, fp);
  fclose(fp);
  return 0;
}
