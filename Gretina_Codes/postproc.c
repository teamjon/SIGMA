/* program to post-process interactions from decomposed event signals for GRETINA / GRETA.

   Author:  D.C. Radford    Dec 2006
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "gdecomp.h"


int main(int argc, char **argv)
{
  double xfinal[MAX_PARS], yfinal[MAX_PARS], zfinal[MAX_PARS], efinal[MAX_PARS];
  double chisq, chisq2, total_energy, t0;
  int    i, number, event_number, found, photopk, nints, ppflag;
  char   line[256];
  char   in_fn[256] = "output/test_gdecomp.out", out_fn[256] = "output/postproc.out";
  FILE   *file_in = 0, *file_out = 0;

  Interaction   ints[2*MAX_SEGS];       /* interactions */
  int           nseg;                   /* number of segments that fired in the event */

  /* process possible command line options:
     -i <input_file_name>
     -o <output_file_name>
  */
  for (i=1; i<argc; i++) {
    if (*argv[i] == '-') {
      if (strstr(argv[i], "i")) {
	strncpy(in_fn, argv[++i], sizeof(in_fn));
      } else if (strstr(argv[i], "o")) {
	strncpy(out_fn, argv[++i], sizeof(out_fn));
      }
    }
  }

  /* open input interactions file for reading */
  if (!(file_in=fopen(in_fn, "r"))) {
    printf("ERROR -- Cannot open file %s for reading.\n", in_fn);
    return 1;
  }
  /* open output interactions file for writing */
  if (!(file_out=fopen(out_fn, "w+"))) {
    printf("ERROR -- Cannot open file %s for writing.\n", out_fn);
    return 1;
  }

  /* read in basis signal data */
  if (read_basis(BASIS_FILE)) return 1;

  /* loop over events */
  printf("processing...\n");

  for (number = 0; ; number++) {
    /* read in an event */
    if (!fgets(line, sizeof(line), file_in)) break;
    sscanf(line, " %d %d %lf %lf %lf %lf %d",
	   &nseg, &nints, &total_energy, &t0, &chisq, &chisq2, &event_number);
    if (!fgets(line, sizeof(line), file_in)) break;
    for (i=0; i<nints; i++) {
      fgets(line, sizeof(line), file_in);
      sscanf(line, " %d %lf %lf %lf %lf",
	     &ints[i].seg, &ints[i].r, &ints[i].p, &ints[i].z, &ints[i].e);
      ints[i].e /= total_energy;
    }

    /* post_process the results */
    // ppflag = 2; // combine interactions in the same seg based on energy or distance
    ppflag = 1; // do not combine interactions in the same seg based on energy or distance
    found = postprocess_events(ints, nints, total_energy, ppflag,
			       COAL_DIST_DEFAULT, xfinal, yfinal, zfinal, efinal);

    /* write results to file_out */
    photopk = 0; // no longer used; flag is available for reuse?
    fprintf(file_out,
	    " %2d %2d %9.2f  %6.2f %24.6f %9.4f %4d\n",
	    found, photopk, total_energy, t0, chisq, chisq2, event_number);
    for (i=0; i<found; i++) {
      fprintf(file_out,
	      " %7.2f %7.2f %7.2f %7.1f\n",
	      xfinal[i], yfinal[i], zfinal[i], efinal[i]);
    }
  }

  printf("  %d events processed.\n", number);
  free(basis);
  return 0;
} /* main */
