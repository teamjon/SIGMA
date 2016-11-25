/* program to calculate weighting potentials of segments
   in inverted coax pc Ge detectors by relaxation
   author:           D.C. Radford
   first written:    April 2013
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define MAX_L 162
#define MAX_R 72

#define RRH (RH + TH*(LH-L+z)/LH)
#define RTR (z <= L-TL ? R : R - TW*(TL+z-L)/TL)  // FIXME - rounding??

#define D2R    0.01745329
#define SIN(x) sin(D2R*(float)x)
#define COS(x) cos(D2R*(float)x)

int main(int argc, char **argv)
{
  /* --- default values, normally over-ridden by values in a *.conf file --- */
  /* all dimensions are in unts of 0.5mm */
  int   R  = 70;        // radius of detector
  int   L  = 160;       // length of detector
  int   TL = 120;       // length of radial tapered part of crystal
  int   TW = 21;        // width/amount of radial taper (decrease in radius) of crystal
  int   RC = 6;         // radius of central contact
  int   LC = 0;         // length of central contact
  int   RO = 30;        // radius of wrap-around outer (Li) contact
  int   LO = 0;         // length of ditch next to wrap-around outer (Li) contact
  int   WO = 0;         // width of ditch next to wrap-around outer (Li) contact
  int   RH = 10;        // radius of hole in outer (Li) contact
  int   TH = 0;         // radial taper (increase in r from inside to outside) of hole
  int   LH = 110;       // length of hole in wrap-around outer (Li) contact

  /* all dimensions are in unts of 0.5mm */
  int   SG   = 1;       // segment gap
  int   SR1  = 65;      // segs 1-8 outer; inner = RO
  int   SL9  = 20;      // seg  9
  int   SL10 = 19;      // seg 10
  int   SL11 = 19;      // seg 11
  int   SL12 = 19;      // seg 12
  int   SL13 = 19;      // seg 13
  int   SL14 = 19;      // seg 14
  int   SL15 = 19;      // seg 15
  int   SL16 = 19;      // seg 16
  int   SRA  = 45;      // seg 17 outer / 16 inner
  int   SRB  = 30;      // seg 18 outer / 17 inner
  int   SRC  = 15;      // seg 19 outer / 18 inner

  int   seg, sz_start[20], sz_end[20], sz_len[20];
  /* ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  --- */

  float *v[2][360][MAX_L+1], dif, sum_dif, max_dif, mean;
  float eps[MAX_L+1][MAX_R+1], eps_dr[MAX_L+1][MAX_R+1], eps_dz[MAX_L+1][MAX_R+1];
  float s1[MAX_R], s2[MAX_R], sa[MAX_R], eps_sum, v_sum;
  int   i, a, r, z, iter, old, new, aa, rr;
  char  line[256];
  FILE  *file;
  double s3, s4;
  float dx, dy, d, xx[360][MAX_R+1], yy[360][MAX_R+1];

  /* initialize */
  for (i=0; i<2; i++) {
    for (a=0; a<360; a++) {
      for (z=0; z<MAX_L+1; z++) {
	v[i][a][z] = malloc((MAX_R+1) * sizeof(float));
      }
    }
  }

  sz_len[ 9] = SL9;
  sz_len[10] = SL10;
  sz_len[11] = SL11;
  sz_len[12] = SL12;
  sz_len[13] = SL13;
  sz_len[14] = SL14;
  sz_len[15] = SL15;
  sz_len[16] = SL16;

  if (argc%2 != 1) {
    printf("Possible options:\n"
	   "      -c config_file_name\n");
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

	    (!strncmp(line, "SR1", 3) && (1 != sscanf(line+3," %d", &SR1))) ||
	    (!strncmp(line, "SRA", 3) && (1 != sscanf(line+3," %d", &SRA))) ||
	    (!strncmp(line, "SRB", 3) && (1 != sscanf(line+3," %d", &SRB))) ||
	    (!strncmp(line, "SRC", 3) && (1 != sscanf(line+3," %d", &SRC))) ||

	    (!strncmp(line, "SL", 2) && ((1 != sscanf(line+2,"%d", &i) ||
					  i < 9 || i > 16 ||
					  1 != sscanf(line+4," %d", &sz_len[i])))) ||

	    (!strncmp(line, "SG", 2) && (1 != sscanf(line+2," %d", &SG)))) {
	  printf("ERROR in config file %s\n"
		 "   ...line is: %s", argv[i+1], line);
	  return 1;
	}
      }
      fclose(file);
    } else {
      printf("Possible options:\n"
	     "      -c config_file_name\n");
      return 1;
    }
  }

  sz_start[9] = 0;
  sz_end[9] = sz_start[9] + sz_len[9];
  for (i=10; i<17; i++) {
    sz_start[i] = sz_end[i-1] + SG;
    sz_end[i] = sz_start[i] + sz_len[i];
  }

  if (RH <= 0.0 && TH <= 0.0 && LH <= 0.0) {
    LH = -1.0;
  }

  if (RO <= 0.0 || RO >= R) {
    RO = R + 2;
    printf("\n\n"
	   "    Crystal: Radius x Length: %.1f x %.1f mm\n"
	   "Outer Taper: Length %.1f mm at %.1f mm/cm\n"
	   "       Hole: Radius x length: %.1f x %.1f mm, taper %.1f mm/cm\n"
	   "No wrap-around contact or ditch...\n\n",
	   0.5 * (float) R, 0.5 * (float) L,
	   0.5 * (float) TL, 10.0 * (float) TW / (float) TL,
	   0.5 * (float) RH, 0.5 * (float) LH, 0.5 * (float) TH / (float) LH);
  } else {
    printf("\n\n"
	   "    Crystal: Radius x length: %.1f x %.1f mm\n"
	   "Outer Taper: Length %.1f mm at %.1f mm/cm\n"
	   "       Hole: Radius x length: %.1f x %.1f mm, taper %.1f mm/cm\n"
	   "Wrap-around: Radius x ditch x gap:  %.1f x %.1f x %.1f mm\n\n",
	   0.5 * (float) R, 0.5 * (float) L,
	   0.5 * (float) TL, 10.0 * (float) TW / (float) TL,
	   0.5 * (float) RH, 0.5 * (float) LH, 0.1 * (float) TH / (float) LH,
	   0.5 * (float) RO, 0.5 * (float) LO, 0.5 * (float) WO);
  }

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

  /* weighting values for the relaxation alg. as a function of r
     These effectively correct for the different volumes of the 
     neighboring voxels, and the different distance separations of the voxel centers,
     both of which vary as a function of r
     They rely on delta-r = delta-z = 1 unit (0.5mm in this case)
     Added sa[r] for azimuthal segmentation, delta-phi = 1 degree <<== NOTE
     They are relative to the delta-z side, for which the volumes are
     always equal and the separation is always 1 unit.

     field * area, summed over all 6 surfaces, is zero
     field = delta-V / distance (dx, dr, or r*dphi)

     s1[r] is for r -> r+1, s2[r] is for r -> r-1, sa[r] is for (r,p) -> (r, p+/-1)
  */
  s1[0] = 4.0;  // r=0; R=dr/2 at outer edge; 2*pi*R / pi*R*R = 2*dz/(dr/2) = 4
  s2[0] = 0.0;  // not used
  sa[0] = 0.0;  // not used
  for (r=1; r<R; r++) {
    s1[r] = 1.0 + 0.5 / (float) r;   //  dz*2*pi*(r+0.5) / dr*2*pi*r ; dr = dz
    s2[r] = 1.0 - 0.5 / (float) r;   //  dz*2*pi*(r-0.5) / dr*2*pi*r ; dr = dz
    sa[r] = 180.0 / (3.141593 * (float) r); // dr*dz / dr*2*pi*r/360 ; dz = 1 unit
    sa[r] *= sa[r]; // squared
  }

  /* --------------- calculate the weighting potentials --------------- */

  for (seg = 8; seg < 20; seg++) {
  //for (seg = 8; seg < 9; seg++) {
    printf("\nCalculating weighting potential for seg %d...\n\n", seg);

    // boundary conditions and starting guesses
    // initial wild guess at final potential, and outside contact:
    for (z=0; z<L+1; z++) {
      for (r=0; r<R+1; r++) {
	v[0][0][z][r] = v[1][0][z][r] = 0.0f;
      }
    }
    // point contact:
    for (z=0; z<LC+1; z++) {
      for (r=0; r<RC+1; r++) {
	v[0][0][z][r] = v[1][0][z][r] = 0.0f;
      }
    }
    // segment contact:
    if (seg <= 8) {   // back pie-slices
      for (r = RO; r < SR1; r++) {
	v[0][0][0][r] = v[1][0][0][r] = 1.0f;
	printf("seg %d; z, r: %.1f %.1f\n", seg, 0., r/2.0f);
      }
      v[0][0][0][SR1] = v[1][0][0][SR1] = 0.5f;

    } else if (seg <= 16) {
      if (seg == 9) {  // seg 9 wraps around the back just a little
	v[0][0][0][SR1] = v[1][0][0][SR1] = 0.5f;
	z = 0;
	for (r=SR1+1; r<RTR; r++) {
	  v[0][0][0][r] = v[1][0][0][r] = 1.0f;
	  printf("seg %d; z, r: %.1f %.1f\n", seg, 0., r/2.0f);
	}
      } else {
	z = sz_start[seg] - 1;
	if (z >= 0) v[0][0][z][RTR] = v[1][0][z][RTR] = 0.5f;
      }
      for (z = sz_start[seg]; z < sz_end[seg]; z++) {  // cylindrical outside
	v[0][0][z][RTR] = v[1][0][z][RTR] = 1.0f;
	printf("seg %d; z, r: %.1f %.1f\n", seg, z/2.0f, RTR/2.0f);
      }
      if (seg == 16) {  // seg 16 wraps around the front just a little
	z = L;
	for (r=RTR; r>SRA; r--) {
	  v[0][0][z][r] = v[1][0][z][r] = 1.0f;
	  printf("seg %d; z, r: %.1f %.1f\n", seg, z/2.0f, r/2.0f);
	}
	v[0][0][L][SRA] = v[1][0][L][SRA] = 0.5f;
      } else {
	z = sz_end[seg];
	if (z >= 0) v[0][0][z][RTR] = v[1][0][z][RTR] = 0.5f;
      }

    } else if (seg == 17) { // outer front ring
      v[0][0][L][SRA] = v[1][0][L][SRA] = 0.5f;
      v[0][0][L][SRB] = v[1][0][L][SRB] = 0.5f;
      for (r=SRA-1; r>SRB; r--) {
	v[0][0][L][r] = v[1][0][L][r] = 1.0f;
	printf("seg %d; z, r: %.1f %.1f\n", seg, L/2.0f, r/2.0f);
      }
    } else if (seg == 18) { // inner front ring
      v[0][0][L][SRB] = v[1][0][L][SRB] = 0.5f;
      v[0][0][L][SRC] = v[1][0][L][SRC] = 0.5f;
      for (r=SRB-1; r>SRC; r--) {
	v[0][0][L][r] = v[1][0][L][r] = 1.0f;
	printf("seg %d; z, r: %.1f %.1f\n", seg, L/2.0f, r/2.0f);
      }
    } else if (seg == 19) { // core
      v[0][0][L][SRC] = v[1][0][L][SRC] = 0.5f;
      for (r=SRC-1; r>RH; r--) {  // core wraps around front just a little
	v[0][0][L][r] = v[1][0][L][r] = 1.0f;
	printf("seg %d; z, r: %.1f %.1f\n", seg, L/2.0f, r/2.0f);
      }
      for (z=L; z>L-LH; z--) {
	v[0][0][z][r] = v[1][0][z][r] = 1.0f;
	printf("seg %d; z, r: %.1f %.1f\n", seg, z/2.0f, r/2.0f);
      }
      for (r=RH; r>=0; r--) {
	v[0][0][z][r] = v[1][0][z][r] = 1.0f;
	printf("seg %d; z, r: %.1f %.1f\n", seg, z/2.0f, r/2.0f);
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

	  v_sum = v[old][0][z+1][r]*eps_dz[z][r] + v[old][0][z][r+1]*eps_dr[z][r]*s1[r];
	  eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r];
	  if (z > 0) {
	    v_sum += v[old][0][z-1][r]*eps_dz[z-1][r];
	    eps_sum += eps_dz[z-1][r];
	  } else {
	    v_sum += v[old][0][z+1][r]*eps_dz[z][r];  // reflection symm around z=0
	    eps_sum += eps_dz[z][r];  // CHECKME / FIXME; use 1 instead?
	  }
	  if (r > 0) {
	    v_sum += v[old][0][z][r-1]*eps_dr[z][r-1]*s2[r];
	    eps_sum += eps_dr[z][r-1]*s2[r];
	  }

	  mean = v_sum / eps_sum;
	  v[new][0][z][r] = mean;
	  dif = v[old][0][z][r] - v[new][0][z][r];
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

    // write WP values to file seg_wp_<nn>.dat, on 0.5-mm grid
    sprintf(line, "seg_wp_%2.2d.dat", seg);
    file = fopen(line, "w");
    fprintf(file, "## r (mm), z (mm), WP\n");
    for (r=0; r<R+1; r++) {
      for (z=0; z<L+1; z++) {
	fprintf(file, "%6.1f %6.1f %10.6f\n",
		((float) r)/2.0f,  ((float) z)/2.0f, v[new][0][z][r]);
      }
      fprintf(file, "\n");
    }
    fclose(file);

    if (seg == 8) {  /* do azimuthal segmentation */

      for (a=0; a<360; a++) {
	for (r=0; r<R+1; r++) {
	  xx[a][r] = r*SIN(a);
	  yy[a][r] = r*COS(a);
	}
      }

      /* boundary conditions and initial estimates */
      for (z=0; z<L+1; z++) {
	for (r=0; r<RTR; r++) {
	  if (r == 0) {  // along axis
	    s3 = v[new][0][z][r];
	    for (a=0; a<360; a++) {
	      v[0][a][z][r] = v[1][a][z][r] = s3/8.0f;
	    }
	  } else if ((z == 0 && r >= RO) ||      // outside contact
		     (z >= L) ||                 // outside contact
		     (z >= L-LH && r <= RRH) ||  // outside contact hole
		     (z < LC+1 && r < RC+1)) {   // inside contact
	    /* boundary conditions */
	    s3 = v[new][0][z][r];
	    for (a=1; a<45; a++) {
	      v[0][a][z][r] = v[1][a][z][r] = s3;
	    }
	    for (a=46; a<360; a++) {
	      v[0][a][z][r] = v[1][a][z][r] = 0.0f;
	    }
	    v[0][ 0][z][r] = v[1][ 0][z][r] = s3/2.0;
	    v[0][45][z][r] = v[1][45][z][r] = s3/2.0;
	  }
	}
      }

      for (z=0; z<L; z++) {
	printf("z=%d\n", z);
	for (r=0; r<RTR; r++) {
	  if (r == 0 ||
	      (z == 0 && r >= RO) ||      // outside contact
	      (z >= L-LH && r <= RRH) ||  // outside contact hole
	      (z < LC+1 && r < RC+1)) {   // inside contact
	    continue;
	  } else {
	    /* initial estimate for bulk */
	    for (a=23; a<203; a++) {
	      s3 = s4 = 0;
	      for (rr=RO; rr<SR1; rr++) {
		for (aa=1; aa<45; aa++) {
		  dx = xx[a][r] - xx[aa][rr];
		  dy = yy[a][r] - yy[aa][rr];
		  // d = sqrt(z*z + dx*dx + dy*dy);
		  d = z*z + dx*dx + dy*dy;
		  s3 += ((double) rr)/d;
		}
		for (aa=46; aa<360; aa++) {
		  dx = xx[a][r] - xx[aa][rr];
		  dy = yy[a][r] - yy[aa][rr];
		  // d = sqrt(z*z + dx*dx + dy*dy);
		  d = z*z + dx*dx + dy*dy;
		  s4 += ((double) rr)/d;
		}
		dx = xx[a][r] - xx[0][rr];
		dy = yy[a][r] - yy[0][rr];
		// d = sqrt(z*z + dx*dx + dy*dy);
		d = z*z + dx*dx + dy*dy;
		s3 += 0.5*((double) rr)/d;
		s4 += 0.5*((double) rr)/d;
		dx = xx[a][r] - xx[45][rr];
		dy = yy[a][r] - yy[45][rr];
		// d = sqrt(z*z + dx*dx + dy*dy);
		d = z*z + dx*dx + dy*dy;
		s3 += 0.5*((double) rr)/d;
		s4 += 0.5*((double) rr)/d;
	      }
	      v[old][a][z][r] = v[new][0][z][r] * s3 / (s3+s4);
	      if (a > 22 && a <= 45)
		v[old][45-a][z][r] = v[old][a][z][r];
	      if (a > 45 && a < 203)
		v[old][405-a][z][r] = v[old][a][z][r];
	    }
	    for (a=0; a<360; a++) {
	      v[new][a][z][r] = v[old][a][z][r];
	    }
	  }
	}
      }

      /* can use this code to save initial estimate for comparison with final result */
      /* */
      // write WP values to file seg_wp_<nn>.dat, on 0.5-mm grid
      sprintf(line, "seg_wp_07.dat");
      file = fopen(line, "w");
      fprintf(file, "## r (mm), theta (deg),  z (mm), WP\n");
      for (r=0; r<R+1; r++) {
	for (a=0; a<360; a++) {
	  for (z=0; z<L+1; z++) {
	    fprintf(file, "%6.1f %3d %6.1f %10.6f\n",
		    ((float) r)/2.0f, a, ((float) z)/2.0f, v[new][a][z][r]);
	  }
	  fprintf(file, "\n");
	}
      }
      fclose(file);
      /* */

      /* this part is just for making 2d projections to look at with gnuplot */
      /*
      for (j=2; j<6; j++) {
	a = 22 + (j-2)*90;
	sprintf(line, "seg_wp_%2.2d.dat", j);
	file = fopen(line, "w");
	fprintf(file, "## r (mm), theta (deg),  z (mm), WP\n");
	for (r=0; r<R+1; r++) {
	  for (z=0; z<L+1; z++) {
	    fprintf(file, "%6.1f %3d %6.1f %10.6f\n",
		    ((float) r)/2.0f, a, ((float) z)/2.0f, v[new][a][z][r]);
	  }
	  fprintf(file, "\n");
	}
	fclose(file);
      }
      */

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
	  for (r=1; r<RTR; r++) {  // exclude r=0 (axis)
	    if (z == 0 && r >= RO) continue;      // outside contact
	    if (z >= L-LH && r <= RRH) continue;  // outside contact hole
	    if (z < LC+1 && r < RC+1) continue;   // inside contact

	    for (a=23; a<203; a++) {

	      v_sum = v[old][a][z+1][r]*eps_dz[z][r] +
	              v[old][a][z][r+1]*eps_dr[z][r]*s1[r] +
		      v[old][a][z][r-1]*eps_dr[z][r-1]*s2[r] +
		     (v[old][a-1][z][r] + v[old][a+1][z][r])*eps[z][r]*sa[r];
	      eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r] +
		        eps_dr[z][r-1]*s2[r] + 2.0f*eps[z][r]*sa[r];
	      if (z > 0) {
		v_sum += v[old][a][z-1][r]*eps_dz[z-1][r];
		eps_sum += eps_dz[z-1][r];
	      } else {
		v_sum += v[old][a][z+1][r]*eps_dz[z][r];  // reflection symm around z=0
		eps_sum += eps_dz[z][r];
	      }

	      mean = v_sum / eps_sum;
	      v[new][a][z][r] = mean;
	      dif = v[old][a][z][r] - v[new][a][z][r];
	      if (dif < 0.0f) dif = -dif;
	      sum_dif += dif;
	      if (max_dif < dif) max_dif = dif;
	    }
	    v[new][22][z][r] = v[new][23][z][r];
	    v[new][203][z][r] = v[new][202][z][r];
	  }
	}
	// report results for some iterations
	if (iter < 10 || (iter < 600 && iter%10 == 0) || iter%50 == 0)
	  printf("%5d %d %d %.10f %.10f\n", iter, old, new, max_dif, sum_dif/(float) (L*R*180));
	if (max_dif == 0.0f) break;
      }
      printf(">> %d %.16f\n\n", iter, sum_dif);
      for (z=0; z<L+1; z++) {
	for (r=0; r<RTR; r++) {
	  for (a=23; a<=45; a++) {
	    v[new][45-a][z][r] = v[new][a][z][r];
	  }
	  for (a=46; a<203; a++) {
	    v[new][405-a][z][r] = v[new][a][z][r];
	  }
	}
      }

      // write WP values to file seg_wp_01.dat, on 0.5-mm grid
      sprintf(line, "seg_wp_01.dat");
      file = fopen(line, "w");
      fprintf(file, "## r (mm), theta (deg),  z (mm), WP\n");
      for (r=0; r<R+1; r++) {
	for (a=0; a<360; a++) {
	  for (z=0; z<L+1; z++) {
	    fprintf(file, "%6.1f %3d %6.1f %10.6f\n",
		    ((float) r)/2.0f, a, ((float) z)/2.0f, v[new][a][z][r]);
	  }
	  fprintf(file, "\n");
	}
      }
    }

  }

  return 0;
}
