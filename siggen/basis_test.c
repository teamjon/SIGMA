#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "pdecomp.h"

int main(int argc, char **argv) {
  PDecomp pdd;
  int     i, j, k, start, end, check, check2, seg = -1;
  char    str[256];
  FILE    *output;
  output = fopen("basis_135_180.txt","a+");

  if (!read_basis(BASIS_FILE, &pdd)) printf("Success!\n");
  for (i=1;i<GRID_PTS;i++) {
//  for (i=1000;i<10000;i++) {
    int x = pdd.basis[i].x;
    int y = pdd.basis[i].y;
    int z = pdd.basis[i].z;
    int ir = pdd.basis[i].ir;
    float ip = pdd.basis[i].ip * (M_PI/60.0 - 0.000001);
    int iz = pdd.basis[i].iz;
    int iseg = pdd.basis[i].iseg;
    int t80 = pdd.basis[i].t80;
    int t20 = pdd.basis[i].t20;
    int t_hi = *pdd.basis[i].hi_time;
    int t_lo = *pdd.basis[i].lo_time;
    int t_drift[GRID_PTS];
    check = check2 = start = end = 0;

    for (j=0;j<TIME_STEPS_C;j++) {
      float test = pdd.basis[i].signal[0][j];
      if (test > 0 && check == 0) {
        start = j;
        check++;
      }
      if (test > 0.95 && check2 == 0) {
        end = j;
        check2++;
      }
      if (start != 0 && end != 0) {t_drift[i] = end - start;}
    }

    float trace[NUM_SIGS][TIME_STEPS_C];
    memset(trace,0,TIME_STEPS_C*NUM_SIGS*sizeof(float));
if(iseg == 17) {
    fprintf(output,"\n%d, %d, %d, %d \n",x,y,z,iseg);
    fprintf(output,"%i, %f, %i\n",ir,ip,iz);
    fprintf(output,"%i \n",t_drift[i]);

    for (k=0;k<NUM_SIGS;k++) {
      for (j=0;j<TIME_STEPS_C;j++) {
        trace[k][j] = pdd.basis[i].signal[k][j];
        float xxx = trace[k][j];
        fprintf(output,"%f, ",xxx);
      }
      fprintf(output,"\n");
    }
}
  }
  fclose(output);
  return 0;
}
