#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include "pdecomp.h"

#include "geant_siggen.h"
#define MG4 "/local/SIGMA_Geant/GANIL/trunk/Marc_analysis/sigma/SigmaResponse.bin"

int G4_pars(void) {
  Marc_G4_Struct mg4s;

  if(fread(&mg4s, sizeof(mg4s), 1, MG4) != 1) return 0;

  int num = mg4s.SGnHits;
  printf("\n Number of hits = %i \n", num);
}
