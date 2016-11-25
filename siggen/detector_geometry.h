/* detector_geometry.h -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module keeps track of the detector geometry
 * See also: assym_detector.h
 */
#ifndef _DETECTOR_GEOMETRY_H
#define _DETECTOR_GEOMETRY_H

#include "point.h"

/* geometry_init
   reads information about detector geometry from file given by geometry_fname
   returns total number of segments, -1 for failure
*/
int geometry_init(char *geometry_fname);

/* geometry_finalize
   Clean up at end of program
*/
int geometry_finalize(void);


/* segment_number
   returns the (geometrical) segment number at point pt, or -1 
   if outside crystal
*/
int segment_number(point pt);

/* zmax_detector
   returns the maximum z value for detector
*/
float zmax_detector(void);

/* rmax_detector
   returns the maximum radius for detector
*/
float rmax_detector(void);




#endif /*#ifndef _DETECTOR_GEOMETRY_H*/
