/* grid_basis.c
 *  David Radford
 *  Initial attempt at a grid for the basis signals for the PSIC detector
 *
 * to compile:
 *  gcc -o grid_basis grid_basis.c point.c cyl_point.c calc_signal.c\
 *    fields_ppc.c detector_geometry.c signal_calc_util.c -lm -lreadline
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <ctype.h>

#include "calc_signal.h"
#include "cyl_point.h"
#include "signal_calc_util.h"
#include "detector_geometry.h"
#include "fields.h"
#include "pdecomp.h"

#define GEOMETRY_FILE "geometry_setup_ppc.dat"
#define FIELDS_DIR_FILE "field_setup_ppc.dat"
#define SIGNAL_PARS_FILE "calc_signal_setup_ppc.dat"

#define RZ_STEP 1.0
#define P_STEP (M_PI/60.0 - 0.000001)

#define CYL 0
#define CART 1

#define CHATTY

static int rc_integrate(float **s_in, float **s_out, float tau);

static int time_steps;
static int nsegments;
static float tau, tau7;

int main(void) {
  float r, p, z, zz, f1, f2, t_drift;
  int   i, j, seg, seg2, t, t80, t20, old_seg, n=0, longest_drift=0;
  int   c1, c2, t1, t2;
  struct point cart;
  struct cyl_pt cyl;
  static float **s;
  FILE   *fp;
  char   header[256];

  Basis_Point bp;

  tau = 3.0;
  tau7 = 4.0;

  //set_signal_calc_output(ANNOYINGLY_CHATTY, vprintf);

  if (signal_calc_init(GEOMETRY_FILE, FIELDS_DIR_FILE,
		       SIGNAL_PARS_FILE, &time_steps, &nsegments) != 0){
    return 1;
  }

  if ((s = malloc(nsegments*sizeof(*s))) == NULL){
    printf("Malloc failed\n");
    return 1;
  }
  for (i = 0; i < nsegments; i++){
    if ((s[i] = malloc(time_steps*sizeof(*s[i]))) == NULL){
      printf("malloc failed\n");
      return 1;
    }
  }

  fp = fopen(BASIS_FILE, "w+");
  snprintf(header, sizeof(header),
	  "Segmented Point Contact basis signals; Cylindrical; version 1.0\n"
	  "%d basis point, %d grid segments, %d time steps\n"
	  "SRAD SPHI SZZZ: %d %d %d\n",
	  GRID_PTS, NUM_SIGS, TIME_STEPS_C, SRAD, SPHI, SZZZ);
  fwrite(header, sizeof(header), 1, fp);

  for (i=0; i<NUM_SIGS; i++) {
    bp.lo_time[i] = 0;
    bp.hi_time[i] = TIME_STEPS_C-1;
  }
  bp.ip = -1;

  printf("  r   phi   z    seg drift ; t80  t20 ;    x     y   number\n");
  for (p = 0*M_PI/4.0 + 0.001/* 0.0f*/; p < 1*M_PI/4.0 + 0.001; p += P_STEP) {
    bp.ip++;
    cyl.phi = p;
    bp.ir = -1;
    for (r = 0.0f; r < 35.1; r += RZ_STEP) {
      bp.ir++;
      cyl.r = r;
      old_seg = -1;
      bp.iz = -1;
      for (z = 0.0f; z < 80.1; z += RZ_STEP) {
	bp.iz++;
	cyl.z = z;

	cart = cyl_to_cart(cyl);
	t = get_signal(cart, s);
        t_drift = t;
	seg = -1;

	/* s[seg][time_steps-1] > 0.9f changed to -s[seg][time_steps-1] > 0.9f
           to account for different polarity signals in SIGMA to LBNL PPC
	*/

	if (t >= 0) {
	  for (seg = 1; seg < nsegments; seg++) {
	    if (-s[seg][time_steps-1] > 0.9f) break;
	  }
	  if (t > longest_drift) longest_drift = t;
	}

	if (seg >= 0 && seg < nsegments) {
	  rc_integrate(s, s, tau);
	  /* find 80% time of PC (seg 0) and 20% time of hit segment */
	  for (t80 = 0; t80 < time_steps && s[0][t80]/s[0][time_steps-1] < 0.8f; t80++) ;
	  f1 = s[0][t80-1]/s[0][time_steps-1];
	  f2 = s[0][t80  ]/s[0][time_steps-1];
	  t80 = 10*(t80-1) + lrintf(10.0f * (0.8f-f1)/(f2-f1));
	  for (t20 = 0; t20 < time_steps && s[seg][t20]/s[seg][time_steps-1] < 0.2f; t20++) ;
	  f1 = s[seg][t20-1]/s[seg][time_steps-1];
	  f2 = s[seg][t20  ]/s[seg][time_steps-1];
	  t20 = 10*(t20-1) + lrintf(10.0f * (0.2f-f1)/(f2-f1));
	  printf("%4.1f %4.2f %5.2f %3d %4d ; %4d %4d ; %4.2f %4.2f  %8d\n",
		 cyl.r, cyl.phi, cyl.z, seg, t, t80, t20, cart.x, cart.y, n++);

	  c1 = c2 = t1 = t2 = 0;
	  for (j=0; j<TIME_STEPS_C; j++) {
	    float test_sig = s[0][j];
	    if (test_sig > 0 && c1 == 0) {
	      t1 = j;
	      c1 ++;
	    }
	    if (test_sig > 0.95 && c2 == 0) {
	      t2 = j;
	      c2 ++;
	    }
	    if (t1 != 0 && t2 != 0) {
	      t_drift = t2 - t1;
	      break;
	    }
	  }

	  bp.iseg = seg;
	  bp.x = cart.x;
	  bp.y = cart.y;
	  bp.z = cart.z;
	  bp.t80 = t80;
	  bp.t20 = t20;
	  bp.t_drift = t;
	  for (t=0; t<TIME_STEPS_C; t++) bp.signal[0][t] = s[0][t];
	  for (i=1; i<NUM_SIGS; i++) {
	    for (t=0; t<TIME_STEPS_C; t++) bp.signal[i][t] = -s[i][t];   //invert neg segment signals
	  }
	  fwrite(&bp, sizeof(bp), 1, fp);

	} else if (old_seg > 0) {
	  // printf("%4.1f %4.2f %4.1f\n", cyl.r, cyl.phi, cyl.z);
	}

	if (old_seg > 0 && seg != old_seg) {
	  for (zz = z - RZ_STEP + 0.01f; zz <= z+0.001f; zz += 0.01f) {
	    cyl.z = zz;
	    cart = cyl_to_cart(cyl);
	    t = get_signal(cart, s);
	    seg2 = -1;
	    if (t >= 0) {
	      for (seg2 = 1; seg2 < nsegments; seg2++) {
		if (s[seg2][time_steps-1] > 0.9f) break;
	      }
	    }
	    if (/* seg2 < nsegments && */ seg2 != old_seg) {
	      rc_integrate(s, s, tau);
	      /* find 80% time of PC (seg 0) and 20% time of hit segment */
	      for (t80 = 0; t80 < time_steps && s[0][t80]/s[0][time_steps-1] < 0.8f; t80++) ;
	      f1 = s[0][t80-1]/s[0][time_steps-1];
	      f2 = s[0][t80  ]/s[0][time_steps-1];
	      t80 = 10*(t80-1) + lrintf(10.0f * (0.8f-f1)/(f2-f1));
	      if (seg > 0) {
		for (t20 = 0; t20 < time_steps && s[seg][t20]/s[seg][time_steps-1] < 0.2f; t20++) ;
		f1 = s[seg][t20-1]/s[seg][time_steps-1];
		f2 = s[seg][t20  ]/s[seg][time_steps-1];
		t20 = 10*(t20-1) + lrintf(10.0f * (0.8f-f1)/(f2-f1));
	      } else {
		t20 = -1;
	      }

	      printf(">> %4.1f %4.2f %5.2f %3d %4d ; %4d %4d ; %4.2f %4.2f\n",
		     cyl.r, cyl.phi, cyl.z, seg, t, t80, t20, cart.x, cart.y);
	      old_seg = seg2;
	    }
	  }
	}

	old_seg = seg;

      }
    }
  }

  printf(">>> Longest drift time: %d\n"
	 ">>> Total of %d points inside detector.\n", longest_drift, n);
  if (n != GRID_PTS) {
    printf("  !!!!!  WARNING !!!!!! GRID_PTS = %d.\n"
	   " Be sure to change value of GRID_PTS in pdecomp.h before reading this basis.\n",
	   GRID_PTS);
    rewind(fp);
    snprintf(header, sizeof(header),
	     "Segmented Point Contact basis signals; Cylindrical; version 1.0\n"
	     "%d basis point, %d grid segments, %d time steps\n"
	     "SRAD SPHI SZZZ: %d %d %d\n",
	     n, NUM_SIGS, TIME_STEPS_C, SRAD, SPHI, SZZZ);
    fwrite(header, sizeof(header), 1, fp);
  }

  bp.ip++; bp.ir++; bp.iz++;
  if (bp.ip != SPHI || bp.ir != SRAD || bp.iz != SZZZ)
    printf("  !!!!!  ERROR !!!!!!\n"
	   " ip = %d, SPHI = %d ;  ir = %d, SRAD = %d ;  iz = %d, SZZZ = %d\n",
	   bp.ip, SPHI, bp.ir, SRAD, bp.iz, SZZZ);

  fclose(fp);
  return 0;
}

static int rc_integrate(float **s_in, float **s_out, float tau){
  int i, j;
  float local_tau;
  float s_in_old, s;  /* DCR: added so that it's okay to
			 call this function with s_out == s_in */

  if (tau < 0.2f || tau7 < 0.2f) return 0;
  for (i = 0; i < nsegments; i++){
    local_tau = tau;
    if (i == 7) local_tau = tau7;
    s_in_old = s_in[i][0];
    s_out[i][0] = 0.0;
    for (j = 1; j < time_steps; j++){
      s = s_out[i][j-1] + (s_in_old - s_out[i][j-1])/local_tau;//Apr 26, KL -- BUG
      s_in_old = s_in[i][j];
      s_out[i][j] = s;
    }
  }
  return 0;
}
