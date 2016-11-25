/* fields_ppc.c -- based on m3d2s.f by I-Yang Lee 
 * Karin Lagergren
 *
 * This module handles the electric field and weighting potential and 
 * calculates drift velocities
 *
 */

//TODO: Add setup_done flag & check it before trying to access data ?
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "point.h"
#include "cyl_point.h"
#include "fields.h"
#include "detector_geometry.h"
#include "calc_signal.h"
#include "signal_calc_util.h"

#define MAX_FNAME_LEN 512

static float rmin, rmax, rstep;
static float zmin, zmax, zstep;
static float phimin, phimax, phistep;
static float ***wpot;
static float ***az_wpot;  // only one segment actually stored
static int   nsegs;
static float temperature;

static int setup_wp(char fnames[20][MAX_FNAME_LEN]);

int inside(float r, float z);

/* calculate segement capacitances, and segment-segment capacitances */
int seg_caps(void) {
  int    i, seg1, seg2, ir, iz, ir2, iz2;
  int    rlen, zlen, alen;
  double sum[50], er, ez, sum2;
  double pi=3.14159, Epsilon=(8.85*16.0/1000.0);  // permittivity of Ge in pF/mm

  rlen = 1.5 + (rmax - rmin)/rstep;  // round up and add 1
  zlen = 1.5 + (zmax - zmin)/zstep;
  alen = 2.5 + (phimax - phimin)/phistep;  // round up and add 2
  printf("rlen, zlen, alen, nsegs = %d %d %d %d\n", rlen, zlen, alen, nsegs);

  for (i=0; i<nsegs; i++) sum[i] = 0.0;
  for (ir=1; ir<rlen-1; ir++) {
    for (iz=1; iz<zlen-1; iz++) {
      if (inside(ir*rstep, iz*zstep) &&
	  (inside((ir+1)*rstep, iz*zstep) && inside((ir-1)*rstep, iz*zstep)) &&
	  (inside(ir*rstep, (iz+1)*zstep) && inside(ir*rstep, (iz-1)*zstep))) {
	for (seg1=9; seg1<nsegs; seg1++) {
	  er = 0.5*(wpot[ir+1][iz][seg1] - wpot[ir-1][iz][seg1]) / rstep;
	  ez = 0.5*(wpot[ir][iz+1][seg1] - wpot[ir][iz-1][seg1]) / zstep;
	  sum[seg1] += (er*er + ez*ez) * 2.0*pi * rstep*(double)ir;
	}
      }
    }
  }
  for (seg2=9; seg2<nsegs; seg2++)
    printf(" %7.3f", sum[seg2] * 0.01 * Epsilon * rstep*rstep*zstep * 4.0*pi);
  // 0.01 converts (V/cm)2 to (V/mm)2,  rstep*rstep*zstep converts grid^3 to  mm3
  // FIXME!!  Why the extra factor of 4*pi? It seems to be needed to get agreement with other method below...
  printf("\n");

  for (i=0; i<nsegs; i++) sum[i] = 0.0;
  for (ir=1; ir<rlen-1; ir++) {
    for (iz=1; iz<zlen-1; iz++) {
      if (inside(ir*rstep, iz*zstep) &&
	  inside((ir+1)*rstep, iz*zstep) &&
	  inside(ir*rstep, (iz+1)*zstep)) {
	for (seg1=9; seg1<nsegs; seg1++) {
	  er = (wpot[ir+1][iz][seg1] - wpot[ir][iz][seg1]) / rstep;
	  ez = (wpot[ir][iz+1][seg1] - wpot[ir][iz][seg1]) / zstep;
	  sum[seg1] += (er*er + ez*ez) * 2.0*pi * rstep*(double)ir;
	}
      }
    }
  }
  for (seg2=9; seg2<nsegs; seg2++)
    printf(" %7.3f", sum[seg2] * 0.01 * Epsilon * rstep*rstep*zstep * 4.0*pi);
  // 0.01 converts (V/cm)2 to (V/mm)2,  rstep*rstep*zstep converts grid^3 to  mm3
  // FIXME!!  Why the extra factor of 4*pi? It seems to be needed to get agreement with other method below...
  printf("\n");
  printf("\n");

  printf("\n"
	 "# seg  capacitances (nF)\n");
  for (seg1=9; seg1<nsegs; seg1++) {
    for (i=0; i<nsegs; i++) sum[i] = 0.0;
    for (ir=0; ir<rlen; ir++) {
      for (iz=0; iz<zlen; iz++) {
	//if (inside(ir*rstep, iz*zstep) && wpot[ir][iz][seg1] == 1.0f) {
	if (wpot[ir][iz][seg1] == 1.0f) {
	  ir2 = ir;
	  iz2 = iz;
	  if (ir < rlen-1 && inside((ir+1)*rstep, iz*zstep) && wpot[ir+1][iz][seg1] < 1.0f) ir2 = ir+1;
	  if (ir > 0      && inside((ir-1)*rstep, iz*zstep) && wpot[ir-1][iz][seg1] < 1.0f) ir2 = ir-1;
	  if (iz < zlen-1 && inside(ir*rstep, (iz+1)*zstep) && wpot[ir][iz+1][seg1] < 1.0f) iz2 = iz+1;
	  if (iz > 0      && inside(ir*rstep, (iz-1)*zstep) && wpot[ir][iz-1][seg1] < 1.0f) iz2 = iz-1;
	  //if (ir < rlen-1 && wpot[ir+1][iz][seg1] < 1.0f) ir2 = ir+1;
	  //if (ir > 0      && wpot[ir-1][iz][seg1] < 1.0f) ir2 = ir-1;
	  //if (iz < zlen-1 && wpot[ir][iz+1][seg1] < 1.0f) iz2 = iz+1;
	  //if (iz > 0      && wpot[ir][iz-1][seg1] < 1.0f) iz2 = iz-1;
	  if (ir2 == ir && iz2 == iz) continue;
	  //printf("%d %.1f %.1f    %.1f %.1f\n", seg1, ir*rstep, iz*zstep, ir2*rstep, iz2*zstep);

	  for (seg2=9; seg2<nsegs; seg2++) {
	    er = (wpot[ir][iz][seg2] - wpot[ir2][iz][seg2]) / rstep;
	    ez = (wpot[ir][iz][seg2] - wpot[ir][iz2][seg2]) / zstep;
	    sum[seg2] += sqrt(er*er + ez*ez) * pi * rstep*(double)(ir+ir2);
	  }
	}
      }
    }
    printf("%2d", seg1);
    sum2 = 0.0;
    for (seg2=9; seg2<nsegs; seg2++) {
      sum2 += sum[seg2];
      printf(" %7.3f", sum[seg2] * 0.1 * Epsilon * rstep*zstep);
      // 0.1 converts (V/cm) to (V/mm),  rstep*zstep converts grid^2 to  mm2
    }
    printf(" %9.3f", sum[seg1] * 0.1 * Epsilon * rstep*zstep);
    printf(" %9.3f", (sum2 - sum[seg1]) * 0.1 * Epsilon * rstep*zstep);
    printf(" %8.3f\n", (2.0*sum[seg1] - sum2) * 0.1 * Epsilon * rstep*zstep);
  }
  printf("\n");
  return 0;

}

/* field_setup
   given a field directory file, read electic field and weighting
   potential tables from files listed in directory
   check that number of segments =1
   returns 0 for success
*/
int field_setup(char *fields_fname, int nsegs_in){
  FILE *fp;
  char line[MAX_LINE];
  char wp_fnames[20][MAX_FNAME_LEN];
  char efield_fname[MAX_FNAME_LEN];
  char velo_fname[MAX_FNAME_LEN];

  // fudge_radial_field_1 = fudge_radial_field_2 = 0.0f;

  tell(NORMAL, "reading field data from file: %s\n", fields_fname);
  if ((fp = fopen(fields_fname, "r")) == NULL){
    error("failed to open field init file: %s\n", fields_fname);
    return -1;
  }
  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f %f %f", &rmin, &rmax, &rstep) != 3){
    error("failed to read r limits, step from file: %s\n", fields_fname);
    fclose(fp);
    return -1;
  }
  tell(NORMAL, "rmin: %.1f rmax: %.1f, rstep: %.1f\n", rmin, rmax, rstep);

  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f %f %f", &zmin, &zmax, &zstep) != 3){
    error("failed to read z limits, step from file: %s\n", fields_fname);
    fclose(fp);
    return -1;
  }
  tell(NORMAL, "zmin: %.1f zmax: %.1f, zstep: %.1f\n", zmin, zmax, zstep);

  if (read_setup_line(fp, line) != 0
      || sscanf(line, "%f %f %f", &phimin, &phimax, &phistep) != 3){
    error("failed to read angle limits, step from file: %s\n", fields_fname);
    fclose(fp);
    return -1;
  }
  if (((int) (phimax/phistep + 1.5f)) % 8 != 0){
    error("angle step (%f) not a divisor of 45 degrees in %s\n", phistep, fields_fname);
    fclose(fp);
    return -1;
  }
  phimin  /= 180.0/3.14159265;  // degrees to radians
  phimax  /= 180.0/3.14159265;
  phistep /= 180.0/3.14159265;

  tell(NORMAL, "phimin: %.1f phimax: %.1f, phistep: %.1f\n",
                phimin, phimax, phistep);

  if (read_setup_line(fp, line) != 0 ||
      sscanf(line, "%f", &temperature) != 1){
    error("Failed to read detector temperature from file: %s\n", fields_fname);
    return -1;
  }

  if (read_setup_line(fp, velo_fname) != 0){
    error("failed to read drift velocity data file name: %s\n", 
	  velo_fname);
    fclose(fp);
    return -1;
  }
  if (read_setup_line(fp, efield_fname) != 0){
    error("failed to read electric field data file name: %s\n", 
	  efield_fname);
    fclose(fp);
    return -1;
  }
  for (nsegs=0; nsegs<20; nsegs++) {
    if (read_setup_line(fp, wp_fnames[nsegs]) != 0) break;
    tell(ANNOYINGLY_CHATTY, "Got file name # %d = %s\n", nsegs, wp_fnames[nsegs]);
    /* use only one file for the 8 azimuthal segments (45-degree symmetry) */
    if (nsegs == 1) nsegs += 7;
  }
  if (nsegs != nsegs_in){
    error("Failed to read correct number of file names (%d - 7) for wp from file %s\n",
	  nsegs_in, fields_fname);
    fclose(fp);
    return -1;
  }
  if (setup_wp(wp_fnames) != 0){
    error("failed to read WPs\n");
    fclose(fp);
    return -1;
  }

  fclose(fp);
  return 0;
}


/*setup_wp
  read weighting potential values from files. returns 0 on success*/
static int setup_wp(char fnames[20][MAX_FNAME_LEN]){
  int fno, segno;
  FILE *fp;
  char line[MAX_LINE];
  int i, j, k, npts;
  int rlen, zlen, alen;
  cyl_pt cyl;
  point cart;
  int lineno;
  float wp;
  char *cp;

  rlen = 1.5 + (rmax - rmin)/rstep;  // round up and add 1
  zlen = 1.5 + (zmax - zmin)/zstep;
  alen = 2.5 + (phimax - phimin)/phistep;  // round up and add 2
  printf("rlen, zlen, alen, nsegs = %d %d %d %d\n", rlen, zlen, alen, nsegs);

  if (wpot == NULL){//assuming rlen, zlen, never change as for setup_efld
    if ((wpot = malloc(rlen*sizeof(*wpot))) == NULL){
      error("Malloc failed in setup_wp\n");
      return 1;
    }
    for (i = 0; i < rlen; i++){
      if ((wpot[i] = malloc(zlen*sizeof(*wpot[i]))) == NULL){
	error("Malloc failed in setup_wp\n");
	//NB: memory leak here.
	return 1;
      }
      for (j = 0; j < zlen; j++){
	if ((wpot[i][j] = malloc(nsegs*sizeof(*wpot[i][j]))) == NULL){
	  error("Malloc failed in setup_wp\n");
	  //memory leak again
	  return 1;
	}
	for (segno = 0; segno < nsegs; segno++){
	  wpot[i][j][segno] = 0.0;
	}
      }
    }
  }

  /* only one azimuthal segment is calculated and stored; use 45-deg symmetry */
  if (az_wpot == NULL){//assuming rlen, zlen, alen never change as for setup_efld
    if ((az_wpot = malloc(rlen*sizeof(*az_wpot))) == NULL){
      error("Malloc failed in setup_wp\n");
      return 1;
    }
    for (i = 0; i < rlen; i++){
      if ((az_wpot[i] = malloc(zlen*sizeof(*az_wpot[i]))) == NULL){
	error("Malloc failed in setup_wp\n");
	//NB: memory leak here.
	return 1;
      }
      for (j = 0; j < zlen; j++){
	if ((az_wpot[i][j] = malloc(alen*sizeof(*az_wpot[i][j]))) == NULL){
	  error("Malloc failed in setup_wp\n");
	  //memory leak again
	  return 1;
	}
	for (k = 0; k < alen; k++){
	  az_wpot[i][j][k] = 0.0;
	}
      }
    }
  }

  /*now read the tables*/
  for (fno = segno = 0; segno < nsegs; segno++, fno++){
    /* data in files is for segments 0,1,2, etc*/
    if ((fp = fopen(fnames[fno], "r")) == NULL){
      error("failed to open file: %s\n", fnames[fno]);
      return -1;
    }
    tell(NORMAL, "Reading weighting potential from file: %s\n",
	 fnames[fno]);
    lineno = npts = 0;

    if (segno == 1) {  // azimuthal segment
      while (fgets(line, MAX_LINE, fp) != NULL){
	lineno++;
	for (cp = line; isspace(*cp) && *cp != '\0'; cp++);
	if (*cp == '#' || !strlen(cp)) continue;
	if (sscanf(line, "%f %f %f %f\n",&cyl.r, &cyl.phi, &cyl.z, &wp) != 4){ 
	  error("failed to read azimuthal weighting potential from line %d\n"
		"line: %s", lineno, line);
	  fclose(fp);
	  return 1;
	}
	cyl.phi /= 180.0/3.14159265;  // degrees to radians
	i = 0.5 + (cyl.r - rmin)/rstep; // round up
	j = 0.5 + (cyl.z - zmin)/zstep;
	k = 0.5 + (cyl.phi - phimin)/phistep;
	if (i < 0 || i >= rlen || j < 0 || j >= zlen || k < 0 || k >= alen) continue;
	cart = cyl_to_cart(cyl);
	//if (segment_number(cart) < 0) continue;
	//printf("%d %d %.2f %.2f %.4f \n", i, j, cyl.r, cyl.z, wp); fflush(stdout);
	az_wpot[i][j][k] = wp;
	npts++;
      }

      fno += 7;
      segno += 7;
      if (npts != rlen*zlen*(alen-1)) {
	error("number of points read (%d) != expected (%d*%d*%d)\n",
	      npts, rlen, zlen, (alen-1));
	fclose(fp);
	return 1;
      }
      tell(ANNOYINGLY_CHATTY,
	   "i, j, k; rlen, zlen, alen; npts: %d %d %d; %d %d %d; %d\n",
	   i, j, k, rlen, zlen, alen, npts);
      /* copy phimin values to phimax+phistep */
      for (i=0; i<rlen; i++) {
	for (j=0; j<zlen; j++) {
	  az_wpot[i][j][alen-1] = az_wpot[i][j][0];
	}
      }

    } else {  // rotationally symmetric segments
      while (fgets(line, MAX_LINE, fp) != NULL){
	lineno++;
	for (cp = line; isspace(*cp) && *cp != '\0'; cp++);
	if (*cp == '#' || !strlen(cp)) continue;
	if (sscanf(line, "%f %f %f\n",&cyl.r, &cyl.z, &wp) != 3){ 
	  error("failed to read weighting potential from line %d\n"
		"line: %s", lineno, line);
	  fclose(fp);
	  return 1;
	}
	i = 0.5 + (cyl.r - rmin)/rstep; // round up
	j = 0.5 + (cyl.z - zmin)/zstep;
	if (i < 0 || i >= rlen || j < 0 || j >= zlen) continue;
	cyl.phi = 0;
	cart = cyl_to_cart(cyl);
	//if (segment_number(cart) < 0) continue;
	//printf("%d %d %.2f %.2f %.4f \n", i, j, cyl.r, cyl.z, wp); fflush(stdout);
	wpot[i][j][segno] = wp;
	npts++;
      }
      if (npts != rlen*zlen) {
	error("number of points read (%d) != expected (%d*%d)\n",
	      npts, rlen, zlen);
	fclose(fp);
	return 1;
      }
      tell(ANNOYINGLY_CHATTY,
	   "i, j, rlen, zlen, npts: %d %d %d %d %d\n", i, j, rlen, zlen, npts);
    }
    tell(NORMAL,"Read %d lines from file %s\n", lineno, fnames[fno]);
    fclose(fp);
  }
  tell(NORMAL, "Done reading weighting potential\n");
  return 0;
}
