/* detector_geometry_ppc.c -- for "ppc" geometry
 * Karin Lagergren
 *
 * This module keeps track of the detector geometry
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "detector_geometry.h"
#include "point.h"
#include "calc_signal.h"
#include "signal_calc_util.h"

/*coaxial detector*/
static float zlen; /*z = 0 => front of detector, zlen=>length*/
static float max_r; /*detector radius, constant*/
static float bullet_radius; /*bulletization radius @ front of det (z=zmax)*/
static float contact_l, contact_r; /*length, radius of central contact*/
static float central_hole_l, central_hole_r; /*length and radius of central hole*/
static float taper_z;       /* z where taper begins */
static float taper_dr;      /* reduction in r at z=zlen due to tapering */

static int in_crystal(point pt);
static int nsegs;


/* geometry_init
   reads information about detector geometry from file given by geometry_fname
   returns total number of segments, -1 for failure

   ppc type detectors are assumed.
   Format for init file for crystal:
   float     -- size in z direction, mm
   float     -- radius, mm
   float     -- bulletization radius, mm (@ front of det, z=zmax)

   returns 1 (one segment)
*/
int geometry_init(char *geometry_fname){
  FILE *fp;
  char line[MAX_LINE];

  if ((fp = fopen(geometry_fname, "r")) == NULL){
    error("failed to open geometry configuration file: %s\n", geometry_fname);
    return -1;
  }
  if (read_setup_line(fp,line) != 0
      || sscanf(line, "%f", &zlen) != 1){
    error("failed to read z length from %s\n", geometry_fname);
    fclose(fp);
    return -1;
  }
  if (read_setup_line(fp,line) != 0
      || sscanf(line, "%f", &max_r) != 1){
    error("failed to read detector radius from %s\n", geometry_fname);
    fclose(fp);
    return -1;
  }
  if (read_setup_line(fp,line) != 0
      || sscanf(line, "%f", &bullet_radius) != 1){
    error("failed to read bulletization radius from %s\n", geometry_fname);
    fclose(fp);
    return -1;
  }
  if (read_setup_line(fp,line) != 0
      || sscanf(line, "%f %f", &contact_l, &contact_r) != 2){
    error("failed to read contact dimensions from %s\n", geometry_fname);
    fclose(fp);
    return -1;
  }
  if (read_setup_line(fp,line) != 0
      || sscanf(line, "%f %f", &central_hole_l, &central_hole_r) != 2){
    error("failed to read central hole size from %s\n", geometry_fname);
    fclose(fp);
    return -1;
  }

  if (read_setup_line(fp, line) != 0) {
    taper_z = zlen + 1.0f;
  } else if (sscanf(line, "%f", &taper_z) != 1) {
    error("Failed to read taper start z from file: %s\n",geometry_fname);
    fclose (fp);
    return -1;
  }
  if (read_setup_line(fp, line) != 0) {
    taper_dr = 0;
  } else if (sscanf(line, "%f", &taper_dr) != 1){
    error("Failed to read outer surface taper r from file: %s\n",geometry_fname);
    fclose (fp);
    return -1;
  }

  if (read_setup_line(fp, line) != 0) {
    nsegs = 1;
  } else if (sscanf(line, "%d", &nsegs) != 1){
    error("Failed to read number of signals from file: %s\n", geometry_fname);
    fclose (fp);
    return -1;
  }

  /* ------------------------- additional segment contact info ----------
  if (read_setup_line(fp, line) != 0) {
    nseg_phi = 0;
  } else if (sscanf(line, "%d", &nseg_phi) != 1){
    error("Failed to read # azimuthal segments from file: %s\n",geometry_fname);
    fclose (fp);
    return -1;
  } else {
    if ((seg_no_phi = malloc(360*sizeof(*seg_no_phi))) == NULL){
      error("malloc failed\n");
      exit(1);
    }
    for (i = 0; i < 360; i++){
      seg_no_phi[i] = (int)(nseg_phi*i/360.0);
    }
  }
  if (read_setup_line(fp, line) != 0) {
    seg_start_r = max_r;
  } else if (sscanf(line, "%f", &seg_start_r) != 1){
    error("Failed to read azimuthal segmentation inner r from file: %s\n",geometry_fname);
    fclose (fp);
    return -1;
  }

  if (read_setup_line(fp, line) != 0) {
    nsegs_tip = zlen;
  } else if (sscanf(line, "%d", &nsegs_tip) != 1){
    error("Failed to read number of segments @ narrow end from file: %s\n",geometry_fname);
    fclose (fp);
    return -1;
  }

  if ((zmax_segment = malloc(zlen*sizeof(*zmax_segment))) == NULL){
    error("Malloc failed\n", geometry_fname);
    exit(1);
  }
  if (read_setup_line(fp, line) != 0) {
    nseg_z = 0;
  } else {
    cp = line;
    for (nseg_z = 0; nseg_z < zlen; nseg_z++){
      zmax_segment[nseg_z] = strtof(cp, &ep);
      if (ep == cp) break;
      if (zmax_segment[nseg_z] <= 0 || zmax_segment[nseg_z] > zlen){
	error("invalid segment width: %f\n",zmax_segment[nseg_z]);
	fclose(fp);
	return -1;
      }
      cp = ep;
    }
  }
  if (nseg_z == zlen){
    error("too many segments in z direction\n");
    fclose(fp);
    return -1;
  }
  for (i = 1; i < nseg_z; i++) zmax_segment[i] += zmax_segment[i-1];
  if (zmax_segment[nseg_z-1] > zlen){
    error("sum of segment thicknesses is larger than detector size\n");
    fclose(fp);
    return -1;
  }
  if ((seg_no_z = malloc((int)(zlen+1)*sizeof(*seg_no_z))) == NULL){
    error("Realloc failed\n");
    exit(1);
  }
  j = 0;
  for (i = 0; i < nseg_z; i++) {
    for (; j < zmax_segment[i] && j < zlen+1; j++) {
      seg_no_z[j] = i;
    }
  }
  for (; j < zlen+1; j++) {
    seg_no_z[j] = nseg_z-1;
  }

  ncontacts = nseg_phi + nseg_z + nsegs_tip + 2;
  printf("ncontacts: %d\n", ncontacts);
  free(zmax_segment);
 ----------------------------------------------------------- */

  tell(ANNOYINGLY_CHATTY, "r: %.2f  z: %.2f b_r %.2f ch_l %.2f ch_r %.2f\n",
       max_r, zlen, bullet_radius, central_hole_l, central_hole_r);
  if (taper_z < zlen)
    tell(ANNOYINGLY_CHATTY, "taper: %.2f by %.2f\n",
	 zlen - taper_z, taper_dr);
  if (nsegs > 1)
    tell(ANNOYINGLY_CHATTY, "  %d segments (including PC)\n", nsegs);

  fclose (fp);
  return nsegs;
}

/* geometry_finalize
   Clean up at end of program
*/
int geometry_finalize(void){
  zlen = max_r = 0;
  return 0;
}


/* segment_number
   returns the (geometrical) segment number at point pt, or -1
   if outside crystal
*/
int segment_number(point pt){
  if (!in_crystal(pt)) return -1;
  return 0;
}

#define SQ(x) ((x)*(x))

/*returns 0 (false) or 1 (true) depending on whether pt is inside the crystal*/
static int in_crystal(point pt){
  float r, z, br;

  z = pt.z;
  if (z >= zlen || z < 0) return 0;

  r = sqrt(SQ(pt.x)+SQ(pt.y));
  if (r > max_r) return 0;
  if (z > taper_z &&
      r > max_r - (taper_dr * (z - taper_z) / (zlen - taper_z))) return 0;

  br = bullet_radius;
  if (z > zlen - br){
    if (r > (max_r - br) + sqrt(SQ(br)- SQ(z-(zlen - br))))
      return 0;
  }
  if (contact_r > 0){
    if (z <= contact_l && r <= contact_r)
      return 0;
  }
  if (central_hole_r > 0){
    if (z >= zlen - central_hole_l
	&& r <= central_hole_r){
      return 0;
    }
  }
  return 1;
}
#undef SQ

/* zmax_detector
   returns the maximum z value for detector
*/
float zmax_detector(void){
  return zlen;
}

float rmax_detector(void){
  return max_r;
}



