/* program to convert GEANT-simulated interactions and events
   into simulated event signals for testing the GRETINA / GRET
   A decomposition algorithm.

   Author:  D.C. Radford    Aug 2004
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gdecomp.h"

#define NEVENTS 5000   /* number of events to process*/
#define DET_FILE "signals/det_1MeV.dat"
#define SIG_FILE "signals/sig_1MeV.dat"

int main(int argc, char **argv)
{
  double x, y, z, e, x1, y1, x2, y2;
  int    i, j, event_number, column, hit, nhit;
  FILE   *det_file, *sig_file;

  Event_Signal event;

  /* seed random numbers */
  srand48(1954);

  /* open input (det) file for reading */
  if (!(det_file=fopen(DET_FILE, "r"))) {
    printf("ERROR -- Cannot open file %s for reading.\n", DET_FILE);
    return 1;
  }
  /* open output (sig) file for writing */
  if (!(sig_file=fopen(SIG_FILE, "w+"))) {
    printf("ERROR -- Cannot open file %s for writing.\n",SIG_FILE);
    return 1;
  }

  /* read in basis signal data */
  if (read_basis(BASIS_FILE)) return 1;

  printf("Processing events from %s...\n", DET_FILE);
  /* loop over events */
  for (event_number = 0; event_number < NEVENTS; event_number++) {

    /* read in the event */
    if (2 != fscanf(det_file, "%d %f",
		    &nhit, &event.total_energy)) break;
    // printf("%5d %d %f\n", event_number, nhit, event.total_energy);
    printf("%5d\r", event_number); fflush(stdout);
    if (nhit>20) {
      printf("\nACKK!! nhit > 20, event number %d\n", event_number);
      return 1;
    }
    for (i=0; i<TOT_SEGS; i++) {
      event.seg_energy[i] = 0.0;
      for (j=0; j<TIME_STEPS; j++) {
 	event.signal[i][j] = 7.0 * (drand48()-0.5) /event.total_energy ;
      }
    }

    /* loop over hits */
    for (hit=0; hit<nhit; hit++) {
      if (6 != fscanf(det_file, "%lf %lf %lf %lf %lf %lf",
		      &x1, &y1, &z, &e, &x, &y)) break;
      // printf("       %f %f %f %f\n", x1, y1, z, e);
      pars[0] = x;
      pars[1] = y;
      pars[2] = z;
      pars[3] = e/event.total_energy;

      /* figure which segment column was hit */

      for (column=1; column<7; column++) {
	if (x1 < -0.001 ||
	    y1 < -0.001 ||
	    y1 > x1 * 1.73205) {
	  /* rotate 60 deg about z to find set of 6 segments */
	  x2 = x1 * 0.5 + y1 * 0.8660254;
	  y2 = y1 * 0.5 - x1 * 0.8660254;
	  x1 = x2;
	  y1 = y2;
	} else {
	  break;
	}
      }
      if (column == 6) column = 0;
      if (column > 6) {
	printf("\nACKK!! column > 6, event number %d\n", event_number);
	return 1;
      }

      /* calculate the expected signals for this interaction */
      int_seg[0] = column * 6;
      npars = 4;
      eval_int_pos();
      for (i=0; i<TOT_SEGS; i++) {
	for (j=0; j<TIME_STEPS; j++) {
	  event.signal[i][j] += bsig.signal[i][j];
	}
	if (bsig.signal[i][TIME_STEPS-1] > 0.5*pars[3]) {
	  // printf("        >> segment %d\n", i);
	  event.seg_energy[i] += e;
	}
      }
    }

    /* write results to sig_file */
    fwrite(&event, sizeof(Event_Signal), 1, sig_file);

  }
  printf("\n Finished.\n");
  return 0;
} /* main */
