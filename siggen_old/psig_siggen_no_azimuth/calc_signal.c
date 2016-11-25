/* calc_signal.c -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module contains the main interface to the signal calculation
 * code. 
 *
 * To use: 
 * -- call signal_calc_init. This will initialize geometry, fields,
 *       drift velocities etc.
 * -- call hit_segment/get_signal
 *
 */
/* TODO: see FIXME's below
   charge_trapping is just a placeholder ATM. Should it be defined
   in the fields module?
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "calc_signal.h"
#include "point.h"
#include "detector_geometry.h"
#include "fields.h"
#include "signal_calc_util.h"
#ifndef SYMMETRIC_CRYSTAL
#include "assym_detector.h"
#endif

#define HOLE_CHARGE 1.0
#define ELECTRON_CHARGE -1.0

/* module global variables*/
static int nsegments; /*number of segments, including central contact*/
static float tstep_sz; /* length of time step, ns*/
static int ntime_steps; /*number of time steps used in calculations*/
static int ntsteps_out; /*number of time steps in output signal*/
static float zmax; /*max possible z value = detector length*/

static point *dpath_e, *dpath_h;

/* prototypes for module-private functions*/
static int segment_max_wp(float *wp, float thresh);
static int make_signal(point pt, float **signal, float q);
static int sum_signal(float **s);//modifies s
static float charge_trapping(point pt, float distance, float q);//trapping
static int setup_signal_calc(char *fname);
static int zero_signal(float **signal);
static int rc_differentiate(float **s_in, float **s_out, float tau);
static int rc_integrate(float **s_in, float **s_out, float tau);

int drift_time;


/* signal_calc_init
   given names of files that define the geometry (as required by the
   geometry module), the electric field and weighting potential,
   and the signals (time step sizes, for instance), initialize
   the calculation, and set the variables pointed to by time_steps 
   and nsegments to number of time steps in output signal and total 
   number of segments, respectively
   returns 0 for success
*/
int signal_calc_init(char *geometry_fname, char *fields_fname, 
		     char *signal_pars_fname, int *time_steps,
		     int *nsegs){
  tell(TERSE, "Reading detector geometry from file: %s\n", geometry_fname);
  if ((nsegments = geometry_init(geometry_fname)) <= 0){
    error("setup of detector geometry failed\n");
   return -1;
  }
  zmax = zmax_detector();
  tell(TERSE, "Reading signal calculation parameters from file: %s\n", 
       signal_pars_fname);
  if (setup_signal_calc(signal_pars_fname) != 0){
    error("signal calculation setup failed\n");
    return -1;
  }
  tell(TERSE, "Reading field configuration data from file: %s\n", 
       fields_fname);
  if (field_setup(fields_fname, nsegments) != 0){
    error("Field setup failed\n");
    return -1;
  }
  *nsegs = nsegments;
  *time_steps = ntsteps_out;
  
  if ((dpath_e = malloc(ntime_steps*sizeof(*dpath_e))) == NULL
      ||(dpath_h = malloc(ntime_steps*sizeof(*dpath_h))) == NULL){
    error("Malloc failed\n");
    return -1;
  }

  tell(TERSE, "Setup of signal calculation done\n");
  return 0;
}

/* setup_signal_calc
   read information about time steps etc from file -- config. information
   specific for this module
   returns 0 for success
*/
static int setup_signal_calc(char *fname){
  FILE *fp;
  char line[MAX_LINE];

  if ((fp = fopen(fname, "r")) == NULL){
    error("Failed to open signal config file: %s\n", fname);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%d", &ntime_steps) != 1){
    error("Failed to read maximum number of time steps to be used in calculations from file: %s\n", 
	  fname);
    fclose(fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
	|| sscanf(line, "%f", &tstep_sz) != 1){
    error("Failed to read size of time steps from file: %s\n", fname);
    fclose(fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
	|| sscanf(line, "%d", &ntsteps_out) != 1){
    error("Failed to read number of time steps in output signal from file:: %s\n", 
	  fname);
    fclose(fp);
    return -1;
  }
  tell(NORMAL, "Will use %d time steps in calculations, each %.2f ns long\n", 
       ntime_steps, tstep_sz);
  tell(NORMAL, "The output signals will have %d time steps\n", 
       ntsteps_out);

  fclose(fp);
  return 0;
}

/* hit_segment
   return the segment number for interaction at point pt
   returns -1 if point is outside crystal
*/
int hit_segment(point pt){
  return get_signal(pt, NULL);
}

/* get_signal
   calculate signal for point pt. Result is placed in signal_out array
   which is guaranteed to have the appropriate size (nsegments * ntsteps_out)
   returns segment number or -1 if outside crystal
   if signal_out == NULL => no signal is stored
*/
int get_signal(point pt, float **signal_out){
  static float **signal;
  static int tsteps = 0, nsegs = 0;
  int initial_segment;
  float **tmp = NULL;
  char tmpstr[MAX_LINE];
  int segment;
  int i, j;
  int comp_f;

  /*first time -- allocate signal array*/
  if (nsegs != nsegments || tsteps != ntime_steps){
    if (nsegs != 0) tmp = signal;
    if ((signal = malloc(nsegments*sizeof(*signal))) == NULL){
      error("malloc failed in hit_segment\n");
      return -1;
    }
    for (j = 0; j < nsegments; j++){
      if (tsteps != 0) free(tmp[j]);
      if ((signal[j] = malloc(ntime_steps*sizeof(*signal[j]))) == NULL){
	return -1;
      }
    }
    if (tmp) free(tmp);
    nsegs = nsegments;
    tsteps = ntime_steps;
  }

  zero_signal(signal);
  initial_segment = segment_number(pt);
  if (initial_segment < 0) return initial_segment;
  pt_to_str(tmpstr, MAX_LINE, pt);
  tell(ANNOYINGLY_CHATTY, 
       " Signal for %s, segment %d\n", tmpstr, initial_segment);

  segment = -1;

  memset(dpath_e, 0, ntime_steps*sizeof(*dpath_e));
  memset(dpath_h, 0, ntime_steps*sizeof(*dpath_h));

  pt_to_str(tmpstr, MAX_LINE, pt);
  tell(ANNOYINGLY_CHATTY, " @@@@@ %s\n", tmpstr);
  drift_time = 0;
  if (make_signal(pt, signal, HOLE_CHARGE) 
      || make_signal(pt, signal, ELECTRON_CHARGE)) return -1;
  /* make_signal returns 0 for success */
  sum_signal(signal);

  /*figure out which segment has net charge*/
#ifndef PPC
  for (i = 0; i < nsegments-1; i++){
#ifdef WEIRD_COAXIAL_GEOMETRY
    if (-signal[i][ntime_steps-1] > NET_SIGNAL_THRESH){
#else
    if (signal[i][ntime_steps-1] > NET_SIGNAL_THRESH){
#endif
      if (segment >= 0){
	error("found more than one segment with net charge\n");
	return -1;
      }
      segment = i;
    }
  }
#else
  if (signal[0][ntime_steps-1] > NET_SIGNAL_THRESH ||
      /* DCR allow for n-type PC detector, where signal is negative */
      signal[0][ntime_steps-1] < -NET_SIGNAL_THRESH) {
    segment = 0;
  }
#endif
  if (signal_out != NULL){
    /*now, compress the signal and place it in the signal_out array*/
    comp_f = ntime_steps/ntsteps_out;
    for (i = 0; i < nsegments; i++){
      for (j = 0; j < ntsteps_out; j++){
	signal_out[i][j] = 0;
      }
      /*I truncate the signal if ntime_steps % ntsteps_out != 0*/
      for (j = 0; j < ntsteps_out*comp_f; j++){
	signal_out[i][j/comp_f] += signal[i][j]/comp_f;
      }
    }
  }
  return drift_time;
}


static int zero_signal(float **signal){
  int i, j;
  
  for (i = 0; i < nsegments; i++){
    for (j = 0; j < ntime_steps; j++){
      signal[i][j] = 0.0;
    }
  }
  return 0;
}


/* make_signal
   Generates the signal originating at point pt, for charge q
   returns 0 for success
*/
static int make_signal(point pt, float **signal, float q){
  int initial_segment, segment;
  char tmpstr[MAX_LINE];
  static float *wpot, *wpot_old, *dwpot;
  static int wp_size;
  point new_pt, prev_pt;
  int t, n;
  vector v, dx;
  float dist;
  int largest_wp_seg;
  int i,j;
  int dv_res;

  if (wp_size != nsegments){/*first time called*/
    if (wpot != NULL) free(wpot);
    if (wpot_old != NULL) free(wpot_old);
    if (dwpot != NULL) free(dwpot);
    if ((wpot = malloc(nsegments*sizeof(*wpot))) == NULL
	|| (wpot_old = malloc(nsegments*sizeof(*wpot_old))) == NULL
	|| (dwpot = malloc(nsegments*sizeof(*dwpot))) == NULL){
      error("malloc failed in make_signal\n");
      exit(1); //OK?
    }
    wp_size = nsegments;
  }

  initial_segment = segment_number(pt);
  prev_pt = new_pt = pt;
  for (t = 0; (dv_res = drift_velocity(new_pt, q, &v)) >= 0; t++){ 
    if (q > 0) dpath_h[t] = new_pt;
    else dpath_e[t] = new_pt;
    tell(ANNOYINGLY_CHATTY, "pt: (%.2f %.2f %.2f), ", new_pt.x,new_pt.y, new_pt.z);
    tell(ANNOYINGLY_CHATTY, "v: (%e %e %e)\n", v.x, v.y, v.z);
    if (t >= ntime_steps - 2){
      tell(ANNOYINGLY_CHATTY,
	   "Exceeded maximum number of time steps (%d).\n", ntime_steps);
      return -1; // FIXME DCR: calculate signal for ntime_steps?
    }
    if (wpotentials(new_pt, wpot) != 0) {
      pt_to_str(tmpstr, MAX_LINE, new_pt);
      tell(NORMAL,
	   "Can calculate velocity but not weighting potentials at %s!\n",
	   tmpstr);
      return -1;
    }
    for (i = 0; i < nsegments; i++){
      if (t > 0) signal[i][t] += q*(wpot[i] - wpot_old[i]);
      wpot_old[i] = wpot[i];
    }
    dx = vector_scale(v, tstep_sz);
    prev_pt = new_pt;
    new_pt = vector_add(new_pt, dx);
    dist = vector_length(dx);
    q = charge_trapping(new_pt, dist, q);//FIXME
  }
  if (t == 0) {
    pt_to_str(tmpstr, MAX_LINE, pt);
    tell(ANNOYINGLY_CHATTY, "The starting point %s is outside the field.\n", tmpstr);
    return -1;
  }
#ifndef PPC
#ifndef WEIRD_COAXIAL_GEOMETRY
  /*check if we have drifted out of back of detector (not to contact)*/
  if (new_pt.z >= zmax){
      tell(ANNOYINGLY_CHATTY,
	   "Drifted out of back end of detector.\n");
      return -1; 
  }
#endif
#endif
  pt_to_str(tmpstr, MAX_LINE, new_pt);
  segment = segment_number(new_pt);
  tell(ANNOYINGLY_CHATTY, 
       "Drifted to edge of field grid, point: %s segment: %d q: %.2f\n", 
       tmpstr, segment, q);

  /* now we are outside the electric grid. figure out how much we must
    drift to get to the crystal boundary*/
  for (n = 0; segment >= 0 && n+t < ntime_steps; n++){
    new_pt = vector_add(new_pt, dx);
    segment = segment_number(new_pt);
    if (q > 0) dpath_h[t+n] = new_pt;
    else dpath_e[t+n] = new_pt;
  }
  if (n == 0) n = 1;/*always drift at least one more step*/
  // tell(ANNOYINGLY_CHATTY, 
  tell(NORMAL, 
       "q: %.1f t: %d n: %d ((%.2f %.2f %.2f)=>(%.2f %.2f %.2f))\n", 
       q,t,n, pt.x, pt.y, pt.z, new_pt.x, new_pt.y, new_pt.z);
  if (drift_time < t) drift_time = t;
  if (n + t >= ntime_steps){
    tell(ANNOYINGLY_CHATTY,
	 "Exceeded maximum number of time steps (%d)\n", ntime_steps);
    return -1;  /* FIXME DCR: does this happen? could this be improved? */
  }
  /*weighting pot. is 1 at edge for hit segment, 0 for other segments.
    Make it so, gradually if applicable*/
  largest_wp_seg = segment_max_wp(wpot, WP_THRESH);
  for (i = 0; i < nsegments; i++){
    dwpot[i] = ((i == largest_wp_seg) - wpot[i])/n;
  }

  /*now drift the final n steps*/
  dx = vector_scale(v, tstep_sz);
  dist = vector_length(dx);
  for (i = 1; i <= n; i++){
    for (j = 0; j < nsegments; j++){
      signal[j][i+t-1] += q*dwpot[j];
    }
    q = charge_trapping(prev_pt,dist, q);//FIXME    
  }
  pt_to_str(tmpstr, MAX_LINE, pt);
  tell(ANNOYINGLY_CHATTY, "q:%.2f pt: %s segment: %d\n", q, tmpstr, segment);

  return 0;
}

/*modifies s. each time step will contain the summed signals of 
  all previous time steps*/
static int sum_signal(float **s){
  int i, j;

  for (i = 0; i < nsegments; i++){
    for (j = 1; j < ntime_steps; j++){
      s[i][j] += s[i][j-1];
    }
  }
  return 0;
}

//FIXME -- placeholder function. Even parameter list is dodgy
static float charge_trapping(point pt, float distance, float q){
  return q;
}

/* segment_max_wp 
 *  Return the segment number corresponding to the largest w.p.
 * make sure no more than one segment has wp higher than thresh
 */
static int segment_max_wp(float *wp, float thresh){
  int n, i;
  int segno;
  float wpmax;

  n = 0;
  for (i = 0; i < nsegments; i++){
    if (wp[i] > thresh) {
      segno = i;
      n++;
      tell(ANNOYINGLY_CHATTY, " segment %d over threshold\n", i);
    }
  }
  if (n == 1) return segno;
  if (n > 1){
    error(" %d segments over threshold. Your weigthing potential is broken\n",
	  n);
    return -1;
  }
  n = 0;
  wpmax = thresh/10; //OK? -- FIXME
  for (i = 0; i < nsegments; i++){
    if (wp[i] > wpmax){
      segno = i;
      wpmax = wp[i];
      n++;
    }
  }
  if (n){
    tell(ANNOYINGLY_CHATTY, "largest wp for segment %d\n", segno);
    return segno;
  }
  tell(ANNOYINGLY_CHATTY, "segment_max_wp: no charge collected!\n");
  return -1;
}

//OK? -- not used ATM, so not tested.
static int rc_integrate(float **s_in, float **s_out, float tau){
  int i, j;
  float s_in_old, s;  /* DCR: added so that it's okay to
			 call this function with s_out == s_in */

  if (tau < 0.2f) return 0;
  for (i = 0; i < nsegments; i++){
    s_in_old = s_in[i][0];
    s_out[i][0] = 0.0;
    for (j = 1; j < ntime_steps; j++){
      s = s_out[i][j-1] + (s_in_old - s_out[i][j-1])/tau;//Apr 26, KL -- BUG
      s_in_old = s_in[i][j];
      s_out[i][j] = s;
    }
  }
  return 0;
}

//OK? -- not used ATM, so not tested.
static int rc_differentiate(float **s_in, float **s_out, float tau){
  int i, j;
  float s_in_old, s;  /* DCR: added so that it's okay to
			 call this function with s_out == s_in */

  if (tau < 0.2f) return 0;
  for (i = 0; i < nsegments; i++){
    s_in_old = s_in[i][0];
    s_out[i][0] = 0.0;
    for (j = 1; j < ntime_steps; j++){
      s = s_out[i][j-1] + s_in[i][j] - s_in_old - s_out[i][j-1]/tau;
      s_in_old = s_in[i][j];
      s_out[i][j] = s;
    }
  }
  return 0;
}


/* signal_calc_finalize
 * Clean up (free arrays, close open files...)
 */
int signal_calc_finalize(void){
  fields_finalize();
  geometry_finalize();
  free(dpath_h);
  free(dpath_e);
  return 0;
}


#ifndef SYMMETRIC_CRYSTAL
#ifndef WEIRD_COAXIAL_GEOMETRY
/*
 * set_detector_type
 * crystal type A or B?
 */
int set_detector_type(int type){
  if (type >= 0 && type < N_CRYSTAL_TYPES){
    set_crystal_geometry(type);
    read_fields(type);
  }
  return get_crystal_geometry();
}
#endif
#endif


int drift_path_e(point **pp){
  *pp = dpath_e;

  return ntime_steps;
}
int drift_path_h(point **pp){
  *pp = dpath_h;

  return ntime_steps;
}
