/* Program to calculate potentials for SIGMA detector by relaxation
   Author:           D.C. Radford
   First Written:    Nov 2007
   Modified          Oct 2016, J Wright to match SIGMA detector

   This version simulates an inverted coax detector
   If the detector is n-type, we first treat it like a p-type, by
   changing the sign of the impurities and BV. This is so that the algorithm
   for finding undepleted regions still works. Then at the end of the process we
   just invert the potentials.
   The detector has a coax-style hole in the outer contact,
   and allows you to wrap the outer contact down around the passivated surface.
   It also allows you to add a taper to the outside of the detector,
   or to the inside of the core hole.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

float fminf(float x, float y);

#define MAX_L 300
#define MAX_R 150

#define RRH (RH + TH*(LH-L+z)/LH)
#define RTR (z <= L-TL ? R : R - TW*(TL+z-L)/TL)

int main(int argc, char **argv)
{
  /* --- default values, normally over-ridden by values in a *.conf file --- */
  int   R  = 70;        // radius of detector, in 0.5mm
  int   L  = 158;       // length of detector, in 0.5mm   - 158 instead of 160 due to charge drifting outside of detector in simulations
  int   TL = 120;       // length of radial tapered part of crystal, in 0.5mm
  int   TW = 21;        // width/amount of radial taper (decrease in radius) of crystal, in 0.5mm
  int   RC = 6;         // radius of central contact, in 0.5mm
  int   LC = 0;         // length of central contact, in 0.5mm
  int   RO = 30;        // radius of wrap-around outer (Li) contact, in 0.5mm
  int   LO = 0;         // length of ditch next to wrap-around outer (Li) contact, in 0.5mm
  int   WO = 0;         // width of ditch next to wrap-around outer (Li) contact, in 0.5mm
  int   RH = 10;        // radius of hole in outer (Li) contact, in 0.5mm
  int   TH = 0;         // radial taper (increase in radius from inside to outside) of hole, in 0.5mm
  int   LH = 108;       // length of hole in wrap-around outer (Li) contact, in 0.5mm
  float BV = -3500.0;   // bias voltage
  float N  = 1.31;      // charge density at z=0 in units of e+10/cm3
  float M  = -0.1;      // charge density gradient, in units of e+10/cm4

  int   WV = 1;         // 0: do not write the V and E values to ppc_ev.dat
                        // 1: write the V and E values to ppc_ev.dat
                        // 2: write the V and E values for both +r, -r (for gnuplot)
  int   WP = 1;         // 0: do not calculate the weighting potential
                        // 1: calculate the WP and write the values to ppc_wp.dat
  /* ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  --- */

  float v[2][MAX_L+1][MAX_R+1], dif, sum_dif, max_dif, a, b, mean, min;
  float eps[MAX_L+1][MAX_R+1], eps_dr[MAX_L+1][MAX_R+1], eps_dz[MAX_L+1][MAX_R+1];
  float s1[MAX_R], s2[MAX_R], E_r, E_z, factor, eps_sum, v_sum, bubble_volts;
  float e_over_E = 0.7072;   // actually (delta^2/4)*1e10*e/epsilon, where delta=0.5mm
  int   i, r, z, iter, old, new, rr;
  int   ntype = 0;
  char  undepleted[MAX_R+1][MAX_L+1], line[256];
  FILE  *file;

  /* initialize */

  if (argc%2 != 1) {
    printf("Possible options:\n"
	   "      -c config_file_name\n"
	   "      -b bias_volts\n"
	   "      -w {0,1,2}    (for WV options)\n"
	   "      -p {0,1}      (for WP options)\n");
    return 1;
  }

  for (i=1; i<argc-1; i+=2) {
    if (strstr(argv[i], "-c")) {
      if (!(file = fopen(argv[i+1], "r"))) {
	printf("\nERROR: config file %s does not exist?\n", argv[i+1]);
	return 1;
      }
      /* read config file */
      printf("\nReading values from config file %s\n", argv[i+1]);
      while (fgets(line, sizeof(line), file)) {
	// printf("%s", line);
	if ((!strncmp(line, "R ", 2) && (1 != sscanf(line+2," %d", &R ))) ||
	    (!strncmp(line, "L ", 2) && (1 != sscanf(line+2," %d", &L ))) ||
	    (!strncmp(line, "TL", 2) && (1 != sscanf(line+2," %d", &TL))) ||
	    (!strncmp(line, "TW", 2) && (1 != sscanf(line+2," %d", &TW))) ||
	    (!strncmp(line, "RC", 2) && (1 != sscanf(line+2," %d", &RC))) ||
	    (!strncmp(line, "LC", 2) && (1 != sscanf(line+2," %d", &LC))) ||
	    (!strncmp(line, "RO", 2) && (1 != sscanf(line+2," %d", &RO))) ||
	    (!strncmp(line, "LO", 2) && (1 != sscanf(line+2," %d", &LO))) ||
	    (!strncmp(line, "WO", 2) && (1 != sscanf(line+2," %d", &WO))) ||
	    (!strncmp(line, "RH", 2) && (1 != sscanf(line+2," %d", &RH))) ||
	    (!strncmp(line, "LH", 2) && (1 != sscanf(line+2," %d", &LH))) ||
	    (!strncmp(line, "TH", 2) && (1 != sscanf(line+2," %d", &TH))) ||
	    (!strncmp(line, "BV", 2) && (1 != sscanf(line+2," %f", &BV))) ||
	    (!strncmp(line, "N ", 2) && (1 != sscanf(line+2," %f", &N ))) ||
	    (!strncmp(line, "M ", 2) && (1 != sscanf(line+2," %f", &M ))) ||
	    (!strncmp(line, "WV", 2) && (1 != sscanf(line+2," %d", &WV))) ||
	    (!strncmp(line, "WP", 2) && (1 != sscanf(line+2," %d", &WP)))) {
	  printf("ERROR in config file %s\n"
		 "   ...line is: %s", argv[i+1], line);
	  return 1;
	}
      }
      fclose(file);
    } else if (strstr(argv[i], "-b")) {
      BV = atof(argv[i+1]);   // bias volts
    } else if (strstr(argv[i], "-w")) {
      WV = atoi(argv[i+1]);   // write-out options
    } else if (strstr(argv[i], "-p")) {
      WP = atoi(argv[i+1]);   // weighting-potential options
    } else {
      printf("Possible options:\n"
	     "      -c config_file_name\n"
	     "      -b bias_volts\n"
	     "      -w {0,1,2}    (for WV options)\n"
	     "      -p {0,1}      (for WP options)\n");
      return 1;
    }
  }
  if (WV < 0 || WV > 2) WV = 0;

  /* make sure bias voltage has correct polarity */
  if (BV < 0) BV *= -1.0f;

  if (RH <= 0.0 && TH <= 0.0 && LH <= 0.0) {
    LH = -1.0;
  }

  if (RO <= 0.0 || RO >= R) {
    RO = R + 2;
    printf("\n\n"
	   "    Crystal: Radius x Length: %.1f x %.1f mm\n"
	   "Outer Taper: Length %.1f mm at %.1f mm/cm\n"
	   "       Hole: Radius x length: %.1f x %.1f mm, taper %.1f mm/cm\n"
	   "No wrap-around contact or ditch...\n"
	   "       Bias: %.0f V\n"
	   " Impurities: (%.3f + %.3fz) e10/cm3\n\n",
	   0.5 * (float) R, 0.5 * (float) L,
	   0.5 * (float) TL, 10.0 * (float) TW / (float) TL,
	   0.5 * (float) RH, 0.5 * (float) LH, 0.5 * (float) TH / (float) LH,
	   BV, N, M);
  } else {
    printf("\n\n"
	   "    Crystal: Radius x length: %.1f x %.1f mm\n"
	   "Outer Taper: Length %.1f mm at %.1f mm/cm\n"
	   "       Hole: Radius x length: %.1f x %.1f mm, taper %.1f mm/cm\n"
	   "Wrap-around: Radius x ditch x gap:  %.1f x %.1f x %.1f mm\n"
	   "       Bias: %.0f V\n"
	   " Impurities: (%.3f + %.3fz) e10/cm3\n\n",
	   0.5 * (float) R, 0.5 * (float) L,
	   0.5 * (float) TL, 10.0 * (float) TW / (float) TL,
	   0.5 * (float) RH, 0.5 * (float) LH, 0.1 * (float) TH / (float) LH,
	   0.5 * (float) RO, 0.5 * (float) LO, 0.5 * (float) WO, BV, N, M);
  }

  /* if detector is n-type, make it p-type for now */
/*  if (N + M * L / 40.0f > 0) { // middle of detector is n-type
    printf(" Detector is n-type.\n\n");
    ntype = 1;
    N *= -1.0f;
    M *= -1.0f;
  }*/

  memset(undepleted, ' ', sizeof(undepleted));
  old = 1;
  new = 0;

  // permittivity
  /* boundary condition at Ge-vacuum interface:
     epsilon0 * E_vac = espilon_Ge * E_Ge
  */
  for (z=0; z<L+1; z++) {
    for (r=0; r<R+1; r++) {
      eps[z][r] =  16;   // inside Ge
      if (z < LO  && r < RO && r > RO-WO-1) eps[z][r] =  1;  // inside vacuum
      if (r > 0) eps_dr[z][r-1] = (eps[z][r-1]+eps[z][r])/2.0f;
      if (z > 0) eps_dz[z-1][r] = (eps[z-1][r]+eps[z][r])/2.0f;
    }
  }

  // boundary conditions and starting guesses
  // initial wild guess at final potential:
  for (z=0; z<L; z++) {
    a = BV * (float) (z) / (float) L;
    for (r=0; r<R; r++) {
      v[0][z][r] =  a + (BV - a) * (float) (r) / (float) R;
    }
  }
  // outside contact:
  for (r=0;  r<R+1; r++) v[0][L][r] = v[1][L][r] = 0.0f;//BV;				//Bottom segments
  for (r=RO; r<R+1; r++) v[0][0][r] = v[1][0][r] = 0.0f;//BV;				//Top segments
  for (z=0;  z<L+1; z++) {
    for (r=RTR; r<R+1; r++) {
      v[0][z][r] = v[1][z][r] = 0.0f;//BV;						//Side segments
    }
  }
  for (z=L-LH; z<L+1; z++) {
    for (r=0; r<RRH+1; r++) {
      v[0][z][r] = v[1][z][r] = 0.0f;//BV;						//Core
    }
  }
  // inside contact:
  for (z=0; z<LC+1; z++) {
    for (r=0; r<RC+1; r++) {
      v[0][z][r] = v[1][z][r] = BV;//0.0f;						//PC
    }
  }

  // weighting values for the relaxation alg. as a function of r
  s1[0] = 2.0;
  s2[0] = 0.0;
  for (r=1; r<R; r++) {
    s1[r] = 1.0 + 0.5 / (float) r;   //  for r+1
    s2[r] = 1.0 - 0.5 / (float) r;   //  for r-1
  }

  // now do the actual relaxation
  for (iter=0; iter<50000; iter++) {
    if (old == 0) {
      old = 1;
      new = 0;
    } else {
      old = 0;
      new = 1;
    }
    sum_dif = 0.0f;
    max_dif = 0.0f;
    bubble_volts = 0.0f;

    for (z=0; z<L; z++) {
      for (r=0; r<RTR; r++) {
	if (z == 0 && r >= RO) continue;      // outside contact
	if (z >= L-LH && r <= RRH) continue;  // outside contact hole
	if (z < LC+1 && r < RC+1) continue;   // inside contact

	/* boundary condition at Ge-vacuum interface:
	   epsilon0 * E_vac = espilon_Ge * E_Ge
	*/
	factor = 1.0;
	undepleted[r][z] = '.';
	if (z < LO && r < RO && r > RO-WO-1) {
	  factor = 0.0;
	  undepleted[r][z] = ' ';
	}
	v_sum = v[old][z+1][r]*eps_dz[z][r] + v[old][z][r+1]*eps_dr[z][r]*s1[r];
	eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r];
	min = fminf(v[old][z+1][r], v[old][z][r+1]);
	if (z > 0) {
	  v_sum += v[old][z-1][r]*eps_dz[z-1][r];
	  eps_sum += eps_dz[z-1][r];
	  min = fminf(min, v[old][z-1][r]);
	} else {
	  v_sum += v[old][z+1][r]*eps_dz[z][r];  // reflection symm around z=0
	  eps_sum += eps_dz[z][r];
	}
	if (r > 0) {
	  v_sum += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];
	  eps_sum += eps_dr[z][r-1]*s2[r];
	  min = fminf(min, v[old][z][r-1]);
	} else {
	  v_sum += v[old][z][r+1]*eps_dr[z][r]*s1[r];  // reflection symm around r=0
	  eps_sum += eps_dr[z][r]*s1[r];
	}

	mean = v_sum / eps_sum;
	v[new][z][r] = mean + factor * (N + M * ((float) z)/20.0) * e_over_E;
	if (v[new][z][r] <= 0.0f) {
	  v[new][z][r] = 0.0f;
	  undepleted[r][z] = '*';
/* 	} else if (v[new][z][r] <= min + 0.05) {   // FIXME */
/* 	  v[new][z][r] = min + 0.05f; */
/* 	  undepleted[r][z] = '*'; */
	} else if (v[new][z][r] < min) {
	  if (bubble_volts == 0.0f) bubble_volts = min + 0.5f;
	  v[new][z][r] = bubble_volts;
	  undepleted[r][z] = '*';
	}
	dif = v[old][z][r] - v[new][z][r];
	if (dif < 0.0f) dif = -dif;
	sum_dif += dif;
	if (max_dif < dif) max_dif = dif;
      }
    }
    // report results for some iterations
    if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0)
      printf("%5d %d %d %.10f %.10f\n", iter, old, new, max_dif, sum_dif/(float) (L*R));
    if (max_dif == 0.0f) break;
  }
  printf(">> %d %.16f\n\n", iter, sum_dif);

  ntype = 1;
  /* for n-type detectors, invert the polarity of the potential */
  if (ntype) {
    for (z=0;  z<L+1; z++) {
      for (r=0; r<R+1; r++) {
	v[new][z][r] *= -1.0f;
      }
    }
  }

  // give V and E along the axes r=0 and z=0
  printf("  z(mm)(r=0)      V   E(V/cm) |  r(mm)(z=0)      V   E(V/cm)\n");
  a = b = v[new][0][0];
  for (z=0; z<L+1; z++) {
    printf("%10.1f %8.1f %8.1f  |",
	      ((float) z)/2.0f, v[new][z][0], 20.0*(v[new][z][0] - a));
    a = v[new][z][0];
    if (z > R) {
      printf("\n");
    } else {
      r = z;
      printf("%10.1f %8.1f %8.1f\n",
	     ((float) r)/2.0f, v[new][0][r], 20.0*(v[new][0][r] - b));
      b = v[new][0][r];
    }
  }
  if (bubble_volts > 0.0f) printf("\nPinch-off at %.0f volts\n\n", bubble_volts);

  // write a little file that shows any undepleted voxels in the crystal
  file = fopen("undepleted.txt", "w");
  for (r=R; r>=0; r--) {
    undepleted[r][L] = '\0';
    fprintf(file, "%s\n", undepleted[r]);
  }
  fclose(file);

  if (WV) {
    // write potential and field to output file
    // use 1.0-mm or 0.5-mm grid, depending on value of WV
    file = fopen("ppc_ev.dat", "w");
    fprintf(file, "## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)\n");

    for (rr = (WV==2?-R:0); rr<R+1; rr+=WV) {
      r = abs(rr);
      for (z=0; z<L+1; z+=WV) {
	// calc E in r-direction
	if (r==0) {
	  E_r = (v[new][z][r] - v[new][z][r+1])/0.05;
	} else if (r==R) {
	  E_r = (v[new][z][r-1] - v[new][z][r])/0.05;
	} else {
	  E_r = (v[new][z][r-1] - v[new][z][r+1])/0.1;
	}
	if (z < LO  && r < RO && r > RO-WO-1) {   // inside vacuum
	  if (r+2 > RO-WO/2) {
	    E_r = 999.0;
	  } else {
	    E_r = -999.0;
	  }
	}
	// calc E in z-direction
	if (z==0) {
	  E_z = (v[new][z][r] - v[new][z+1][r])/0.05;
	} else if (z==L) {
	  E_z = (v[new][z-1][r] - v[new][z][r])/0.05;
	} else {
	  E_z = (v[new][z-1][r] - v[new][z+1][r])/0.1;
	}
	fprintf(file, "%6.1f %6.1f %7.1f %7.1f %7.1f %7.1f\n",
		((float) rr)/2.0f,  ((float) z)/2.0f, v[new][z][r],
		sqrt(E_r*E_r + E_z*E_z), E_r, E_z);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }

  // now calculate the weighting potential for the central contact
  // if so desired

  if (WP) {
    printf("\nCalculating weighting potential...\n\n");

    // boundary conditions and starting guesses
    // initial wild guess at final potential, and outside contact:
    for (z=0; z<L+1; z++) {
      for (r=0; r<R+1; r++) {
	v[0][z][r] = v[1][z][r] = 0.0f;
      }
    }
    // inside contact:
    for (z=0; z<LC+1; z++) {
      for (r=0; r<RC+1; r++) {
	v[0][z][r] = v[1][z][r] = 1.0f;
      }
    }

    // now do the actual relaxation
    for (iter=0; iter<50000; iter++) {
      if (old == 0) {
	old = 1;
	new = 0;
      } else {
	old = 0;
	new = 1;
      }
      sum_dif = 0.0f;
      max_dif = 0.0f;

      for (z=0; z<L; z++) {
	for (r=0; r<RTR; r++) {
	  if (z == 0 && r >= RO) continue;      // outside contact
	  if (z >= L-LH && r <= RRH) continue;  // outside contact hole
	  if (z < LC+1 && r < RC+1) continue;   // inside contact

	  v_sum = v[old][z+1][r]*eps_dz[z][r] + v[old][z][r+1]*eps_dr[z][r]*s1[r];
	  eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r];
	  if (z > 0) {
	    v_sum += v[old][z-1][r]*eps_dz[z-1][r];
	    eps_sum += eps_dz[z-1][r];
	  } else {
	    v_sum += v[old][z+1][r]*eps_dz[z][r];  // reflection symm around z=0
	    eps_sum += eps_dz[z][r];
	  }
	  if (r > 0) {
	    v_sum += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];
	    eps_sum += eps_dr[z][r-1]*s2[r];
	  } else {
	    v_sum += v[old][z][r+1]*eps_dr[z][r]*s1[r];  // reflection symm around r=0
	    eps_sum += eps_dr[z][r]*s1[r];
	  }

	  mean = v_sum / eps_sum;
	  v[new][z][r] = mean;
	  dif = v[old][z][r] - v[new][z][r];
	  if (dif < 0.0f) dif = -dif;
	  sum_dif += dif;
	  if (max_dif < dif) max_dif = dif;
	}
      }
      // report results for some iterations
      if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0)
	printf("%5d %d %d %.10f %.10f\n", iter, old, new, max_dif, sum_dif/(float) (L*R));
      if (max_dif == 0.0f) break;
    }
    printf(">> %d %.16f\n\n", iter, sum_dif);

    // write WP values to ppc_wp.dat, on 0.5-mm grid
    file = fopen("ppc_wp.dat", "w");
    fprintf(file, "## r (mm), z (mm), WP\n");
    for (r=0; r<R+1; r++) {
      for (z=0; z<L+1; z++) {
	fprintf(file, "%6.1f %6.1f %10.6f\n",
		((float) r)/2.0f,  ((float) z)/2.0f, v[new][z][r]);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }

  return 0;
}
