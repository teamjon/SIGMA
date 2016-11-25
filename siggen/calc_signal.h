/* calc_signal.h -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module contains the main interface to the signal calculation
 * code. 
 *
 * To use: 
 * -- call signal_calc_init. This will initialize geometry, fields,
 *       drift velocities etc.
 * -- call hit_segment/get_signal
 */
#ifndef _CALC_SIGNAL_H
#define _CALC_SIGNAL_H

#include <stdarg.h>
#include "point.h"

#define NET_SIGNAL_THRESH 0.55
#define WP_THRESH 0.55
#define WP_THRESH_ELECTRONS 1e-4 /*electrons are considered collected if
				   they stop drifting where the wp is < this*/

typedef struct{
  float **s;
  int *t_lo;
  int *t_hi;
}Signal;

/* signal_calc_init
   given names of files that define the geometry (as required by the
   geometry module), the electric field and weighting potential,
   and the signals (time step sizes, for instance), initialize
   the calculation, and set the variables pointed to by time_steps 
   and nsegments to number of time steps in output signal and total 
   number of segments, respectively
   returns 0 for success
*/
int signal_calc_init(char *geometry_fname, char *fields_dir_fname, 
		     char *signal_pars_fname, int *time_steps,
		     int *nsegments);
/* hit_segment
 * return the segment number for an interaction at point pt
 * -1 if outside crystal
 */
int hit_segment(point pt);

/* get_signal calculate signal for point pt. Result is placed in
 * signal array which is assumed to have at least (number of segments)
 * x ( number of time steps) elements
 * returns segment number or -1 if outside crystal
 */
int get_signal(point pt, float **signal);

/* signal_calc_finalize
 * Clean up
 */
int signal_calc_finalize(void);

#ifndef SYMMETRIC_CRYSTAL
/*
 * set_detector_type
 * crystal type A or B?
 */
int set_detector_type(int type);
#endif

#ifdef TRAPPING   // RJC 04/12/08
float field_trapping(point new_pt, vector v_d, float q);
#endif


/*drift paths for last calculated signal.
  after the call, "path" will point at a 1D array containing the points
  (one per time step) of the drift path. 
  freeing that pointer will break the code.
*/
int drift_path_e(point **path);
int drift_path_h(point **path);

#endif /*#ifndef _CALC_SIGNAL_H*/
