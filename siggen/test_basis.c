#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "pdecomp.h"

int main(int argc, char **argv) {
  PDecomp pdd;
  int     i = 1, j, seg = -1;
  char    str[256];
  FILE    *file;

  if (!read_basis(BASIS_FILE, &pdd)) printf("Success!\n");

  while (i > 0) {
    printf("Basis point number = ? ");
    fgets(str, sizeof(str), stdin);
    if (strstr(str, "seg")) {
      seg = atoi(str+3);
      continue;
    }
    i = atoi(str);
    if (i > 0) {
      file = fopen("j.dat", "w+");
      fwrite (pdd.basis[i].signal[0], sizeof(pdd.basis[0].signal), 1, file);
      printf("x y z t80 t20 = %5.2f %5.2f %5.2f %4d %4d; seg = %d\n",
	     pdd.basis[i].x, pdd.basis[i].y, pdd.basis[i].z,
	     pdd.basis[i].t80, pdd.basis[i].t20, (int) pdd.basis[i].iseg);
      if (seg <= 0 || seg >= NUM_SIGS) {
	j = pdd.t_order[i];
      } else if (i > pdd.seg_basis_num[0][seg]) {
	printf("Only %d sigs for seg %d...\n", pdd.seg_basis_num[0][seg], seg);
	continue;
      } else {
	j = pdd.seg_basis_num[i][seg];
      }
      fwrite (pdd.basis[j].signal[0], sizeof(pdd.basis[0].signal), 1, file);
      printf("x y z t80 t20 = %5.2f %5.2f %5.2f %4d %4d; seg = %d  #%d\n",
	     pdd.basis[j].x, pdd.basis[j].y, pdd.basis[j].z,
	     pdd.basis[j].t80, pdd.basis[j].t20, (int) pdd.basis[j].iseg, j);
      fclose(file);
    }
  }

  return 0;
}
