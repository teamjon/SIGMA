#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "pdecomp.h"

int main(int argc, char **argv) {
  PDecomp pdd;
  int     i, j, k, seg = -1;
  char    str[256];
  FILE    *output;
  output = fopen("database.txt","a+");

  if (!read_basis(BASIS_FILE, &pdd)) printf("Success!\n");
  for (i=1;i<10/*GRID_PTS*/;i++) {
    int x = pdd.basis[i].x;
    int y = pdd.basis[i].y;
    int z = pdd.basis[i].z;
    int ir = pdd.basis[i].ir;
    int ip = pdd.basis[i].ip;
    int iz = pdd.basis[i].iz;
    int iseg = pdd.basis[i].iseg;
    int t80 = pdd.basis[i].t80;
    int t20 = pdd.basis[i].t20;
    int t_hi = *pdd.basis[i].hi_time;
    int t_lo = *pdd.basis[i].lo_time;

    float trace[NUM_SIGS][TIME_STEPS_C];
    memset(trace,0,TIME_STEPS_C*NUM_SIGS*sizeof(int));

    fprintf(output,"\n%d, %d, %d, %d \n",x,y,z,iseg);
    fprintf(output,"%d, %d, %d\n",ir,ip,iz);
    fprintf(output,"%d, %d, %d, %d \n",t20,t80,t_lo,t_hi);

/*    fprintf(output,"\nx = %dmm, y = %dmm, z = %dmm,seg = %d \n",x,y,z,iseg);
    fprintf(output,"r = %dmm, phi = %d rad, z = %dmm\n",ir,ip,iz);
    fprintf(output,"t20 = %dns, t80 = %dns, lo_time = %dns, hi_time = %dns \n",t20,t80,t_lo,t_hi);
*/    for (k=0;k<NUM_SIGS;k++) {
      for (j=0;j<TIME_STEPS_C;j++) {
        trace[k][j] = pdd.basis[i].signal[k][j];
        float xxx = trace[k][j];
        fprintf(output,"%f, ",xxx);
      }
      fprintf(output,"\n");
    }
  }
  fclose(output);
  return 0;
}
