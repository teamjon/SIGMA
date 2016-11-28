#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include "pdecomp.h"

#include "geant_siggen.h"
#define G4DAT "/local/SIGMA_Geant/GANIL/trunk/Marc_analysis/sigma/SigmaResponse.bin"

/* 
	Simple parser script to convert SigmaResponse.bin binary structure to 
	SigG4.dat binary structure to work with SigGen
*/

int main(int type) {
  FILE *file, *fp;
  Marc_G4_Struct mg4;
  SigGen_G4_Struct sg4;
  int iii, e_tot, evn_test=0;

  printf("Reading basis signals from %s\n", G4DAT);				// Open binary file MODE2 for fread
  if (!(file=fopen(G4DAT, "r"))) {
    printf("ERROR -- Cannot open file %s for reading.\n", G4DAT);
  }
  fp = fopen("SigG4.dat", "w+");						// Open binary file SigG4.dat for output struct

  while(!feof(file)) {								// Loop over length of binary file
    fread(&mg4, sizeof(mg4), 1, file);						// Read 1 interaction from binary file - length described by struct Marc_G4_Struct

    int evn_id = mg4.SGevent;
    int nHits = mg4.SGnHits;

    sg4.num = nHits;

    if(evn_id != evn_test) {
      evn_test = evn_id;
      iii = 0;
      e_tot = 0;
    }

    if(evn_id == evn_test) {
      sg4.intpts[iii].x = mg4.SGx;
      sg4.intpts[iii].y = mg4.SGy;
      sg4.intpts[iii].z = mg4.SGz;
      sg4.intpts[iii].e = mg4.SGenergy;
      e_tot += mg4.SGenergy;
      iii++;
    }
    sg4.tot_e = e_tot;

    if(iii == nHits) fwrite(&sg4, sizeof(sg4), 1, fp);
  }
  fclose(file);
  fclose(fp);

  return 1;
}
