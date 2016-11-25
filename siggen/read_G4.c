#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "geant_siggen.h"

int read_G4(FILE *file, SigGen_G4_Struct *m2)
{
  /* routine to read binary data from Marc's Geant
     simulations for use with SigGen

     Author: Jonathan Wright     Nov 2016
  */

  char header[256];
  char buffer [500];
  int i,j;

  if(!m2) return 0;

  if(fread(m2, sizeof(*m2), 1, file)!=1) return 0;
}
