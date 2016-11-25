#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "geant_siggen.h"

int read_mode2(FILE *file, mode2_struct *m2)
{
  /* routine to read mode 2 format binary data from Geant
     simulations for use with SigGen

     Author: Jonathan Wright     Aug 2016
  */

  char header[256];
  char buffer [500];
  int i,j;

  if(!m2) return 0;

  if(fread(&GHd, sizeof(GHd), 1, file)!=1) return 0;

  if(GHd.type == 1){
    if(sizeof(*m2)!=GHd.length*sizeof(char)) return 0;
    if(fread(m2, sizeof(*m2), 1, file)!=1) return 0;
    return GHd.type;
  }
  fseek(file, GHd.length, SEEK_CUR);
  return GHd.type;
}
