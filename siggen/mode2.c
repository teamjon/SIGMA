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
  mode2_struct m2;
  g4Sim_emittedGamma g4;

  int  i,j,k,typ;
  char header[256];                                                                                     // Global header init
  FILE *fp;
  srand(time(NULL));

  fp = fopen("mode2.dat","w+");

  for(j=0; j<NEVNTS; j++) {
    typ = 0;//rand() % 2;

    if(typ == 0) {
      GHd.type 		= 1;
      GHd.length 	= sizeof(m2);
      m2.type		= 1;
      m2.crystal_id	= 1;
      m2.num		= (rand() % 3) + 1;
      m2.tot_e		= 1;
      m2.timestamp	= 0;
      m2.t0		= 0;
      m2.chisq		= 0;
      m2.norm_chisq	= 0;
      for(i=0; i<MAX_INTPTS; i++) {
        if(i < m2.num) {
          m2.intpts[i].x        = (rand() % 700) / 10 - 35;
          m2.intpts[i].y        = (rand() % 700) / 10 - 35;
          m2.intpts[i].z        = (rand() % 800) / 10;
          m2.intpts[i].e        = 1;
          m2.intpts[i].seg      = 10;
          m2.intpts[i].seg_ener = 100;
        }
        if(i >= m2.num) {
          m2.intpts[i].x        = -1;
          m2.intpts[i].y        = -1;
          m2.intpts[i].z        = -1;
          m2.intpts[i].e        = -1;
          m2.intpts[i].seg      = -1;
          m2.intpts[i].seg_ener = -1;
        }
      }
      fwrite(&GHd, sizeof(GHd), 1, fp);
      fwrite(&m2, sizeof(m2), 1, fp);
    }

    if(typ == 1) {
      GHd.type = 11;
      GHd.length = sizeof(g4);
      fwrite(&GHd, sizeof(GHd), 1, fp);
      fwrite(&g4, sizeof(g4), 1, fp);
    }
  }
  fclose(fp);
  return 0;
}
