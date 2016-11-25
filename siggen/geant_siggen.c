/* grid_basis.c
 *  Jonathan Wright
 *  A simple code to combine the output of Geant simulations with the input
 *  of SigGen to simulated "real" data. Reads binary file using "read_mode2.c"
 *  and inputs positions into SigGen func "get_signal". Outputs information from
 *  Geant in addition to pulse shapes and drift times from SigGen.
 *  -> Works for multiple interactions, with pulses scaled using energy
 *  -> Output pulse shapes are for total energy on each segment (multiples)
 *
 * to compile:
 *  gcc -o g_s geant_siggen.c point.c cyl_point.c calc_signal.c fields_ppc.c\
 *     detector_geometry_ppc.c signal_calc_util.c read_mode2.c -lm -lreadline
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
#include "geant_siggen.h"

#define GEOMETRY_FILE "geometry_setup_ppc.dat"
#define FIELDS_DIR_FILE "field_setup_ppc.dat"
#define SIGNAL_PARS_FILE "calc_signal_setup_ppc.dat"

//#define CHATTY

static int rc_integrate(float **s_in, float **s_out, float tau);

static int time_steps;
static int nsegments;
static float tau, tau7;

int main(void) {

  float f1, f2, t_drift, x, y, z, e_corr;
  int   i, j, k, seg, seg2, t, t80, t20, old_seg, n=0, longest_drift=0;
  int   c1, c2, t1, t2;
  struct point cart;
  struct cyl_pt cyl;
  static float **s;
  FILE   *fp, *file;
  char   header[256];

  Basis_Point bp;
  mode2_struct m2;
  SigGen_Output sgo;

  tau = 3.0;
  tau7 = 4.0;

  printf("Reading basis signals from %s\n", MODE2);				// Open binary file MODE2 for fread
  if (!(file=fopen(MODE2, "r"))) {
    printf("ERROR -- Cannot open file %s for reading.\n", MODE2);
    return 1;
  }

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

  fp = fopen(BASIS, "w+");							// Opens output file BASIS for fwrite
  int test,evt_cnt=0;
  old_seg = -1;
  int status = 0;
  while((status = read_mode2(file,&m2))) {
//    if(evt_cnt>10000) break;							// evt_cnt used to test small portion of data
    if(status != 1) continue;							// If event is not type 1, continue to next event
    memset(sgo.signal_mult, 0, sizeof(sgo.signal_mult));
    sgo.type = 25;
    sgo.num = m2.num;
    sgo.tot_e = m2.tot_e;
    sgo.nseg = NUM_SIGS;
    sgo.time_step = TIME_STEPS_C;
    for(j=0; j<MAX_INTPTS; j++) {

      int nint = m2.num;
      cart.x = m2.intpts[j].x;							// Set x,y,z to struct cart
      cart.y = m2.intpts[j].y;
      cart.z = 80 - (m2.intpts[j].z - 40);					// -40 to correct from Heather sim --> 80  - () to change beam from bottom to top

      /* Adds a 2.5 degree rotation to the x and y coords */
/*      int rad = 24;
      float phi = 0;
      cart.x += rad*sin(phi);
      cart.y += rad*cos(phi);*/
      /* ------------------------------------------------ */
      t = get_signal(cart, s);							// Get signal/drift time s from coords cart
      seg = -1;
      if (t >= 0) {
        for (seg = 1; seg < nsegments; seg++) {
          if (s[seg][time_steps-1] > 0.9f) break;
        }
        if (t > longest_drift) longest_drift = t;
      }

      sgo.intpts[j].x = cart.x;
      sgo.intpts[j].y = cart.y;
      sgo.intpts[j].z = cart.z;
      sgo.intpts[j].e = m2.intpts[j].e;
      sgo.intpts[j].seg = m2.intpts[j].seg;
      sgo.intpts[j].seg_ener = m2.intpts[j].seg_ener;
      sgo.intpts[j].t_drift = get_signal(cart,s);
      e_corr = sgo.intpts[j].e;// / sgo.tot_e;

      for (t=0; t<TIME_STEPS_C; t++) {
        bp.signal[0][t] = -s[0][t] * e_corr;                             // invert neg PC signal
        sgo.signal_mult[0][t] += bp.signal[0][t];
      }
      for (k=1; k<NUM_SIGS; k++) {
        for (t=0; t<TIME_STEPS_C; t++) {
          bp.signal[k][t] = s[k][t] * e_corr;
          sgo.signal_mult[k][t] += bp.signal[k][t];
        }
      }
    }	                                                                // end of for j = 1->MAX_INTPTS loop
    fwrite(&sgo, sizeof(sgo), 1, fp);
    evt_cnt++;
  }									// end of while file open loop
  fclose(fp);
  fclose(file);
  return 0;
}

static int rc_integrate(float **s_in, float **s_out, float tau) {
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
