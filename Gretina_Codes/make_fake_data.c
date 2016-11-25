#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "gdecomp.h"

int quiet;               /* set to 1 to prevent output of diagnostic info */

int main(int argc, char **argv)
{
  double r, p, t0, x, y, z, te, e;
  int    seg, s, t, ie, is, ii;
  FILE   *datfile, *parfile;
  Interaction  ints;
  Event_Signal asig;
  Basis_Point  bsig1, bsig2;

  static int all[37] =
    {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

#define NEVENTS 10000   /* number of events to generate */
#define SEGM    2       /* segment multiplicity per crystal*/
#define INTM    2       /* interaction multiplicity per hit segment */


  /* seed random numbers */
  srand48(1954);

  /* read in basis signal data */
  if (read_basis(BASIS_FILE)) return 1;

  /* open output file */
  if (!(datfile=fopen("fake_data.dat", "w"))) {
    printf("ERROR -- Cannot open file fake_data.dat for writing.\n");
    return 1;
  }
  if (!(parfile=fopen("fake_data.txt", "w"))) {
    printf("ERROR -- Cannot open file fake_data.txt for writing.\n");
    return 1;
  }
  printf("\n");

  fprintf(parfile, "* enum, e_tot, tzero, N(seg, x, y, z, e)\n");
  for (ie=0; ie<NEVENTS; ie++) {
    if (ie%100 == 0) {printf("\r %d", ie); fflush(stdout);}
    asig.total_energy = 0.0f;
    for (s=0; s<GRID_SEGS; s++) {
      asig.seg_energy[s] = 0.0f;
      for (t=0; t<TIME_STEPS; t++) {
	asig.signal[s][t] = 8.0 * (drand48()-0.5); /* +- 4 keV noise */
      }
    }

    te = 1000.0;
    t0 = 2.0 * drand48();
    fprintf(parfile, "%5d %5.0f %5.3f\n", ie, te, t0);

    for (is=0; is<SEGM; is++) {
      seg = (int) ((double) TOT_SEGS) * drand48();

      for (ii=0; ii<INTM; ii++) {
	r = ((double) maxir[seg]-1) * drand48();
	p = ((double) maxip[seg]-1) * drand48();
	z = ((double) maxiz[seg]-1) * drand48();
	ints.seg = seg;
	ints.pos = -1;
	ints.r = r;
	ints.p = p;
	ints.z = z;
	cyl_to_cart(seg, &ints.r, &x, &y, &z, &e);
	/* find the nearest neighbor grid points, and interpolate to get
	   the signal and derivatives at (rp, pp, zp) */
	if (interpolate(seg, r, p, z, &bsig1, &bsig2, 0, all)) return 1;

	if (is == SEGM - 1 && ii == INTM - 1) {
	  /* last interaction; take all remaining energy */
	  e = te;
	} else {
	  /* take 1/3 to 2/3 of remaining total energy */
	  e = 0.33 * te * (1.0 + drand48());
	  te -= e;
	}
	fprintf(parfile, "           %2d %5.1f %5.1f %4.1f %4.0f\n", seg, x, y, z, e);

	asig.total_energy += e;
	for (s=0; s<GRID_SEGS; s++) {
	  asig.seg_energy[s] += e;
	  for (t=0; t<TIME_STEPS; t++) {
	    asig.signal[s][t] += bsig1.signal[s][t];
	  }
	}
      }
    }
    fwrite(&asig, sizeof(Event_Signal), 1, datfile);
  }

  printf("\n\n%d events written to fake_data.dat; seg mult = %d, int mult = %d\n",
	 ie, SEGM, INTM);
  fclose(datfile);
  fclose(parfile);

  return 0;
}
