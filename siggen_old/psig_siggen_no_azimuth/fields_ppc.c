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

struct velocity_lookup{
  float e;
  float e100;
  float e110;
  float e111;
  float h100;
  float h110;
  float h111;
#ifdef DRIFT_VEL_ANISOTROPY
  float ea; //coefficients for anisotropic drift 
  float eb;
  float ec;
  float ebp;
  float ecp;
  float ha;
  float hb;
  float hc;
  float hbp;
  float hcp;
#endif
  float hcorr;
  float ecorr;
};

static float rmin, rmax, rstep;
static float zmin, zmax, zstep;
static cyl_pt **efld;
static float ***wpot;
static int nsegs;
static struct velocity_lookup *v_lookup;
static int v_lookup_len;

static int nearest_field_grid_index(cyl_pt pt, cyl_int_pt *ipt);
static int grid_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2]);
static cyl_pt efield(cyl_pt pt, cyl_int_pt ipt);
static int setup_efield(char *fname);
static int setup_wp(char fnames[20][MAX_FNAME_LEN]);
static int setup_velo(char *fname);
static cyl_int_pt field_grid_index(cyl_pt pt);
static float tcorr(float q, float e);
static float temperature = REF_TEMP; //77K is the default/
static float fudge_radial_field_1;
static float fudge_radial_field_2;

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

  fudge_radial_field_1 = fudge_radial_field_2 = 0.0f;

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
      || sscanf(line, "%f", &temperature) != 1){
    error("Failed to read detector temperature from file: %s\n", fields_fname);
    return -1;
  }
  tell(NORMAL, "detector temperature is set to %.1f K\n", temperature);

  if (read_setup_line(fp, velo_fname) != 0
	|| setup_velo(velo_fname) != 0){
    error("failed to read drift velocity data from file: %s\n", 
	  velo_fname);
    fclose(fp);
    return -1;
  }
  if (read_setup_line(fp, efield_fname) != 0
	|| setup_efield(efield_fname) != 0){
    error("failed to read electric field data from file: %s\n", 
	  efield_fname);
    fclose(fp);
    return -1;
  }
  for (nsegs=0; nsegs<20; nsegs++) {
    if (read_setup_line(fp, wp_fnames[nsegs]) != 0) break;
    tell(ANNOYINGLY_CHATTY, "Got file name # %d = %s\n", nsegs, wp_fnames[nsegs]);
  }
  if (nsegs != nsegs_in){
    error("Failed to readcorrect number of file names (%d) for wp from file %s\n",
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

static int efield_exists(cyl_pt pt){
  int segment;
  char ptstr[MAX_LINE];
  int i, j;
  int ir, iz;
  cyl_int_pt ipt;
  int rlen, zlen;
  point cart;

  sprintf(ptstr, "(%.1f,%.1f %.1f)", pt.r, pt.phi, pt.z);
  cart = cyl_to_cart(pt);
  segment = segment_number(cart);
  if (segment < 0){
    tell(ANNOYINGLY_CHATTY, "point %s is outside crystal\n", ptstr);
    return 0;
  }else{
    tell(ANNOYINGLY_CHATTY, "point %s is in crystal\n", ptstr);
  }
  ipt = field_grid_index(pt);
  rlen = (int) rint( (rmax - rmin)/rstep );
  zlen = (int) rint( (zmax - zmin)/zstep );
  
  if (ipt.r < 0 || ipt.r + 1 >= rlen
      || ipt.z < 0 || ipt.z + 1 >= zlen){
    tell(ANNOYINGLY_CHATTY, "point %s is outside wp table\n", ptstr);
    return 0;
  }
  for (i = 0; i < 2 ; i++){
    ir = ipt.r + i;
    for (j = 0; j < 2; j++){
      iz = ipt.z + j;
      if (efld[ir][iz].r == 0.0 && efld[ir][iz].z == 0.0) {
	tell(ANNOYINGLY_CHATTY, "point %s has no efield\n", ptstr);
	return 0;
      }
    }
  }
  return 1;
}


/* wpotentials
   gives (interpolated) weighting potential for each segment
   at point pt. These values are stored in wp, which is assumed to 
   have one element per segment. returns 0 for success, 1 on failure
*/
int wpotentials(point pt, float *wp){
  float w[2][2];
  int i, j, k;
  cyl_int_pt ipt;
  cyl_pt cyl;
  int res;

  cyl = cart_to_cyl(pt);
  res = nearest_field_grid_index(cyl, &ipt);
  if (res < 0) return 1;
  grid_weights(cyl, ipt, w);
  for (k=0; k<nsegs; k++) {
    wp[k] = 0.0;
    for (i = 0; i < 2; i++){
      for (j = 0; j < 2; j++){
	wp[k] += w[i][j]*wpot[ipt.r+i][ipt.z+j][k];
      }
    }
  }

  return 0;
}

/* drift_velocity
   calculates drift velocity for charge q at point pt
   returns 0 on success, 1 on success but extrapolation was necessary,
   and -1 for failure
   anisotropic drift: crystal axes are assumed to be (x,y,z)
*/
int drift_velocity(point pt, float q, vector *velo){
  cyl_pt e, en;
  point cart_en;
  cyl_int_pt ipt;
  int i;
  float abse;
  float f;
  float a, b, c;
  float absv;
  int sign;
  cyl_pt cyl;
  int nfgi_ans;
#ifdef DRIFT_VEL_ANISOTROPY
  float bp, cp;
  float en4, en6;
#else
  float  v100, v110, v111;
#endif
  float tempf;
  
  cyl = cart_to_cyl(pt);
  if ((nfgi_ans = nearest_field_grid_index(cyl, &ipt)) < 0) return -1;
  e = efield(cyl, ipt);
  abse = vector_norm_cyl(e, &en);
  en.phi = cyl.phi;
  cart_en = cyl_to_cart(en);
#ifdef DRIFT_VEL_ANISOTROPY
  /* find location in table to interpolate / extrapolate from*/
  for (i = 0; i < v_lookup_len - 2 && abse > v_lookup[i+1].e; i++);
  if (abse == v_lookup[i].e){// FIXME: why is this a special case? 
                             // slightly faster in rare cases, yes, but why?
    if (q > 0){ /*hole*/
      a = v_lookup[i].ha;
      b = v_lookup[i].hb;
      c = v_lookup[i].hc;
      bp = v_lookup[i].hbp;
      cp = v_lookup[i].hcp;
    }else{
      a = v_lookup[i].ea;
      b = v_lookup[i].eb;
      c = v_lookup[i].ec;
      bp = v_lookup[i].ebp;
      cp = v_lookup[i].ecp;
    }
  }else{/*interpolate/extrapolate*/
    f = (abse - v_lookup[i].e)/(v_lookup[i+1].e - v_lookup[i].e); 
    if (q > 0){
      a = (v_lookup[i+1].ha - v_lookup[i].ha)*f+v_lookup[i].ha;
      b = (v_lookup[i+1].hb- v_lookup[i].hb)*f+v_lookup[i].hb;
      c = (v_lookup[i+1].hc - v_lookup[i].hc)*f+v_lookup[i].hc;
      bp = (v_lookup[i+1].hbp- v_lookup[i].hbp)*f+v_lookup[i].hbp;
      cp = (v_lookup[i+1].hcp - v_lookup[i].hcp)*f+v_lookup[i].hcp;
    }else{
      a = (v_lookup[i+1].ea - v_lookup[i].ea)*f+v_lookup[i].ea;
      b = (v_lookup[i+1].eb- v_lookup[i].eb)*f+v_lookup[i].eb;
      c = (v_lookup[i+1].ec - v_lookup[i].ec)*f+v_lookup[i].ec;
      bp = (v_lookup[i+1].ebp- v_lookup[i].ebp)*f+v_lookup[i].ebp;
      cp = (v_lookup[i+1].ecp - v_lookup[i].ecp)*f+v_lookup[i].ecp;
    }
  }

#define POW4(x) ((x)*(x)*(x)*(x))
#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
  en4 = POW4(cart_en.x) + POW4(cart_en.y) + POW4(cart_en.z);
  en6 = POW6(cart_en.x) + POW6(cart_en.y) + POW6(cart_en.z);
  absv = a + b*en4 + c*en6;
  sign = (q < 0 ? -1 : 1);
  tempf = tcorr(q, abse);
  velo->x = sign*tempf*cart_en.x*(absv + 
				  bp*4*(cart_en.x*cart_en.x - en4) +
				  cp*6*(POW4(cart_en.x) - en6));
  velo->y = sign*tempf*cart_en.y*(absv +
				  bp*4*(cart_en.y*cart_en.y - en4) +
				  cp*6*(POW4(cart_en.y) - en6));
  velo->z = sign*tempf*cart_en.z*(absv +
				  bp*4*(cart_en.z*cart_en.z - en4) +
				  cp*6*(POW4(cart_en.z) - en6));
#undef POW4
#undef POW6
#else //#ifdef DRIFT_VEL_ANISOTROPY
  /* find location in table to interpolate / extrapolate from*/
  for (i = 0; i < v_lookup_len - 2 && abse > v_lookup[i+1].e; i++);
  if (abse == v_lookup[i].e){// FIXME: why is this a special case? 
                             // slightly faster in rare cases, yes, but why?
    if (q > 0){ /*hole*/
      v100 = v_lookup[i].h100;
      v110 = v_lookup[i].h110;
      v111 = v_lookup[i].h111;      
    }else{
      v100 = v_lookup[i].e100;
      v110 = v_lookup[i].e110;
      v111 = v_lookup[i].e111;      
    }
  }else{/*interpolate/extrapolate*/
    f = (abse - v_lookup[i].e)/(v_lookup[i+1].e - v_lookup[i].e); 
    if (q > 0){
      v100 = (v_lookup[i+1].h100 - v_lookup[i].h100) * f + v_lookup[i].h100;
      v110 = (v_lookup[i+1].h110 - v_lookup[i].h110) * f + v_lookup[i].h110;
      v111 = (v_lookup[i+1].h111 - v_lookup[i].h111) * f + v_lookup[i].h111;
    }else{
      v100 = (v_lookup[i+1].e100 - v_lookup[i].e100) * f + v_lookup[i].e100;
      v110 = (v_lookup[i+1].e110 - v_lookup[i].e110) * f + v_lookup[i].e110;
      v111 = (v_lookup[i+1].e111 - v_lookup[i].e111) * f + v_lookup[i].e111;
    }
  }
  /*some matrix algebra (inversion)*/
  a =  0.5 * v100 -  4 * v110 +  4.5 * v111;
  b = -2.5 * v100 + 16 * v110 - 13.5 * v111;
  c =  3.0 * v100 - 12 * v110 +  9.0 * v111;
  /*Approximation:velocity is assumed to be in the direction of the el. field*/
#define POW4(x) ((x)*(x)*(x)*(x))
#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
  absv = a + b * (POW4(cart_en.x) + POW4(cart_en.y) + POW4(cart_en.z))
           + c * (POW6(cart_en.x) + POW6(cart_en.y) + POW6(cart_en.z)); 
#undef POW4
#undef POW6
  tempf = tcorr(q, abse);
  sign = (q < 0 ? -1 : 1);
  *velo = vector_scale(cart_en, sign*absv*tempf);

#endif
  return 0;
}

/* Find (interpolated or extrapolated) electric field for this point */
static cyl_pt efield(cyl_pt pt, cyl_int_pt ipt){
  cyl_pt e, zero = {0,0,0}, ef;
  float w[2][2], ee;
  int i, j;
  static int k=0;

  grid_weights(pt, ipt, w);
  e = zero;
  for (i = 0; i < 2; i++){
    for (j = 0; j < 2; j++){
      ef = efld[ipt.r + i][ipt.z +j];
      e.r += ef.r*w[i][j];
      e.z += ef.z*w[i][j];
    }
  }
  
  // FIXME - Check, is this a bug in the original???
  // e.phi = cyl.phi;  // cyl was declared but never initialized
  e.phi = pt.phi;
  // FIXME - fudge to add quadrople component of field
  if (fabs(fudge_radial_field_1) > 0.01) {
    ee = e.r;
    e.r += fudge_radial_field_1 * (pt.r * (24.0 - pt.r) / 12.0f) *
          (fudge_radial_field_2 + cos(4.0f*pt.phi));
    if (k++ % 50 == 0) {
      tell(ANNOYINGLY_CHATTY, "r phi cos(phi) E, dE: %f %f %f %f %f\n",
	   pt.r, pt.phi, cos(4.0f*pt.phi), e.r, e.r-ee);
    }
  }
  return e;
}


/* Find weights for 8 voxel corner points around pt for e/wp field*/
/* DCR: modified to work for both interpolation and extrapolation */
static int grid_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2]){
  float r,z;

  r = (pt.r - rmin)/rstep - ipt.r;
  z = (pt.z - zmin)/zstep - ipt.z;

  out[0][0] = (1.0 - r) * (1.0 - z);
  out[0][1] = (1.0 - r) *        z;
  out[1][0] =        r  * (1.0 - z);
  out[1][1] =        r  *        z;
  return 0;
}


/*find existing integer field grid index closest to pt*/
/* added DCR */
static int nearest_field_grid_index(cyl_pt pt, cyl_int_pt *ipt){
  /* returns <0 if outside crystal or too far from a valid grid point
              0 if interpolation is okay
              1 if we can find a point but extrapolation is needed
  */
  static cyl_pt  last_pt;
  static cyl_int_pt last_ipt;
  static int    last_ret = -99;
  cyl_pt new;
  int r, z;
  float d[3] = {0.0, -1.0, 1.0};
  point cart;
  float seg;

  if (last_ret != -99 &&
      pt.r == last_pt.r && pt.z == last_pt.z) {
    *ipt = last_ipt;
    return last_ret;
  }
  last_pt = pt;
  last_ret = -2;

  cart = cyl_to_cart(pt);
  if ((seg = segment_number(cart)) < 0) {
    last_ret = -1;
  } else{
    new.phi = 0.0;
    for (z=0; z<3; z++) {
      new.z = pt.z + d[z]*zstep;
      for (r=0; r<3; r++) {
	new.r = pt.r + d[r]*rstep;
	if (efield_exists(new)) {
	  *ipt = last_ipt = field_grid_index(new);
	    if (r == 0 && z == 0) {
	      last_ret = 0;
	    } else {
	      last_ret = 1;
	    }
	    return last_ret;
	}
      }
    }
  }

  return last_ret;
}

/*find integer field grid index corresponding to pt*/
static cyl_int_pt field_grid_index(cyl_pt pt){
  cyl_int_pt ipt;

  ipt.r = (pt.r - rmin)/rstep;
  ipt.phi = 0.0;
  ipt.z = (pt.z - zmin)/zstep;
  return ipt;
}

/* setup_velo
   set up drift velocity calculations (read in table)
*/
static int setup_velo(char *fname){
  char line[MAX_LINE];
  FILE *fp;
  int vlook_sz = 8;
  struct velocity_lookup *tmp;
#ifdef DRIFT_VEL_ANISOTROPY
  struct velocity_lookup v, v0;
  int i;
  float sumb_e, sumc_e, sumb_h, sumc_h;
#endif
  if ((v_lookup = malloc(vlook_sz*sizeof(*v_lookup))) == NULL){
    error("malloc failed in setup_velo\n");
    return -1;
  }
  if ((fp = fopen(fname, "r")) == NULL){
    error("failed to open velocity lookup table file: '%s'\n", fname);
    return -1;
  }
  if (read_setup_line(fp,line) != 0){
    error("Failed to read velocity lookup table from file: %s\n", fname);
    fclose(fp);
    return -1;
  }
  tell(ANNOYINGLY_CHATTY, "Drift velocity table:\n");
  tell(ANNOYINGLY_CHATTY, 
       "  e          e100    e110    e111    h100    h110    h111    ecorr    hcorr\n");   
  for (v_lookup_len = 0; ;v_lookup_len++){
    if (v_lookup_len == vlook_sz -1){
      vlook_sz *= 2;
      if ((tmp = realloc(v_lookup, vlook_sz*sizeof(*v_lookup))) == NULL){
	error("realloc failed in setup_velo\n");
	fclose(fp);
	return -1;
      }
      v_lookup = tmp;
    }
    if (sscanf(line, "%f %f %f %f %f %f %f %f %f", 
	       &v_lookup[v_lookup_len].e,
	       &v_lookup[v_lookup_len].e100,
	       &v_lookup[v_lookup_len].e110,
	       &v_lookup[v_lookup_len].e111,
	       &v_lookup[v_lookup_len].h100,
	       &v_lookup[v_lookup_len].h110,
	       &v_lookup[v_lookup_len].h111,
	       &v_lookup[v_lookup_len].ecorr,
	       &v_lookup[v_lookup_len].hcorr) != 9){
      if (sscanf(line, "%f %f %f %f %f %f %f", 
		 &v_lookup[v_lookup_len].e,
		 &v_lookup[v_lookup_len].e100,
		 &v_lookup[v_lookup_len].e110,
		 &v_lookup[v_lookup_len].e111,
		 &v_lookup[v_lookup_len].h100,
		 &v_lookup[v_lookup_len].h110,
		 &v_lookup[v_lookup_len].h111) != 7){
	break;//asume EOF
      }else{
	//no correction;
	printf("setting tcorr terms to 0\n");
	v_lookup[v_lookup_len].ecorr = v_lookup[v_lookup_len].hcorr = 0;
      }
    }
    //v_lookup[v_lookup_len].e *= 100; /*V/m*/
    tmp = &v_lookup[v_lookup_len];
    /* DCR FUDGE of e mobilities */
    if (0) {
      v_lookup[v_lookup_len].e100 *= 1.135f;
      v_lookup[v_lookup_len].e110 *= 1.045f;
      v_lookup[v_lookup_len].e111 /= 1.045f;
    } else if (1) { // scale to v110
      if (v_lookup[v_lookup_len].e < 999) {
	v_lookup[v_lookup_len].e100 = v_lookup[v_lookup_len].e110 * 10.1f/9.5f;
	v_lookup[v_lookup_len].e111 = v_lookup[v_lookup_len].e110 * 7.8f/9.5f;
      }
    } else if (0) { // scale to v111
      if (v_lookup[v_lookup_len].e < 999) {
	v_lookup[v_lookup_len].e100 = v_lookup[v_lookup_len].e111 * 10.1f/7.8f;
	v_lookup[v_lookup_len].e110 = v_lookup[v_lookup_len].e111 * 9.5f/7.8f;
      }
    }
    /*  ------------------------ */
    tell(ANNOYINGLY_CHATTY, "%10.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
	 tmp->e, tmp->e100, tmp->e110, tmp->e111, tmp->h100, tmp->h110,
	 tmp->h111, tmp->ecorr, tmp->hcorr);
    if (read_setup_line(fp, line) != 0){
      break;
    }
  }
  if (v_lookup_len == 0){
    error("Failed to read velocity lookup table from file: %s\n", fname);
    return -1;
  }  
  v_lookup_len++;
  if (vlook_sz != v_lookup_len){
    if ((tmp = realloc(v_lookup, v_lookup_len*sizeof(*v_lookup))) == NULL){
      error("realloc failed in setup_velo. This should not happen\n");
      fclose(fp);
      return -1;
    }
    v_lookup = tmp;
    vlook_sz = v_lookup_len;
  }
  tell(NORMAL,"Drift velocity table has %d rows of data\n", v_lookup_len);
  fclose(fp);
#ifdef DRIFT_VEL_ANISOTROPY
  for (i = 0; i < vlook_sz; i++){
    v = v_lookup[i];
    v_lookup[i].ea =  0.5 * v.e100 -  4 * v.e110 +  4.5 * v.e111;
    v_lookup[i].eb = -2.5 * v.e100 + 16 * v.e110 - 13.5 * v.e111;
    v_lookup[i].ec =  3.0 * v.e100 - 12 * v.e110 +  9.0 * v.e111;
    v_lookup[i].ha =  0.5 * v.h100 -  4 * v.h110 +  4.5 * v.h111;
    v_lookup[i].hb = -2.5 * v.h100 + 16 * v.h110 - 13.5 * v.h111;
    v_lookup[i].hc =  3.0 * v.h100 - 12 * v.h110 +  9.0 * v.h111;
  }
  v_lookup[0].ebp = v_lookup[0].ecp = v_lookup[0].hbp = v_lookup[0].hcp = 0.0;
  // sumb_e = sumc_e = sumb_h = sumc_h = 0.0;
  sumb_e = -10;  // FIXME Fudge by DCR
  sumc_e =  5;   // FIXME Fudge by DCR
  sumb_h = sumc_h = 0.0;
  for (i = 1; i < vlook_sz; i++){
    v0 = v_lookup[i-1];
    v = v_lookup[i];
    sumb_e += (v.e - v0.e)*(v0.eb+v.eb)/2;
    sumc_e += (v.e - v0.e)*(v0.ec+v.ec)/2;
    sumb_h += (v.e - v0.e)*(v0.hb+v.hb)/2;
    sumc_h += (v.e - v0.e)*(v0.hc+v.hc)/2;
    v_lookup[i].ebp = sumb_e/v.e;
    v_lookup[i].ecp = sumc_e/v.e;
    v_lookup[i].hbp = sumb_h/v.e;
    v_lookup[i].hcp = sumc_h/v.e;
  }
#endif
  return 0;
}


/* This may or may not break if we switch to a non-integer grid*/
/*setup_efield
  read electric field data from file, apply sanity checks
  returns 0 for success
*/
static int setup_efield(char *fname){
  FILE *fp;
  char line[MAX_LINE];
  int i, j;
  int rlen, zlen;
  int lineno;
#ifdef OLD_FORMAT
  float **v, vtmp;
  char **vexists;
#else
  float v, eabs, er, ez;
  char *cp;
#endif
  cyl_pt cyl;
  point cart;

  if ((fp = fopen(fname, "r")) == NULL){
    error("failed to open electric field table: %s\n", fname);
    return 1;
  }
  
  rlen = (int)rint((rmax- rmin)/rstep);
  zlen = (int)rint((zmax - zmin)/zstep);

  if (efld == NULL){// here I assume that r, zlen never
                    // change from their initial values, which is reasonable
    if ((efld = malloc(rlen*sizeof(*efld))) == NULL){
      error("Malloc failed in setup_efield\n");
      fclose(fp);
      return 1;
    }
    for (i = 0; i < rlen; i++){
      if ((efld[i] = malloc(zlen*sizeof(*efld[i]))) == NULL){
	
	error("Malloc failed in setup_efield\n");
	//NB: potential memory leak here.
	fclose(fp);
	return 1;
      }
      memset(efld[i], 0, zlen*sizeof(*efld[i]));
    }
  }
#ifdef OLD_FORMAT
  if ((v = malloc(rlen*sizeof(*v))) == NULL
      || (vexists = malloc(rlen*sizeof(*vexists))) == NULL){
    error("Malloc failed in setup_efield\n");
    fclose(fp);
    return 1;
  }
  for (i = 0; i < rlen; i++){
    if ((v[i] = malloc(zlen*sizeof(*v[i]))) == NULL
	||(vexists[i] = malloc(zlen*sizeof(*vexists[i]))) == NULL){
      
      error("Malloc failed in setup_efield\n");
      //NB: potential memory leak here.
      fclose(fp);
      return 1;
    }
    memset(vexists[i], 0, zlen*sizeof(*vexists[i]));
  }
#endif
  tell(NORMAL, "Reading electric field data from file: %s\n",
       fname);
  lineno = 0;

  /*now read the table*/
  while(fgets(line, MAX_LINE, fp) != NULL){
    lineno++;
#ifdef OLD_FORMAT
    if (sscanf(line, "%f %f %f", &cyl.r, &cyl.z, &vtmp) != 3){
      error("failed to read electric field data from line no %d\n"
	    "of file %s\n", lineno, fname);
      fclose(fp);
      return 1;
    }
#else
    for (cp = line; isspace(*cp) && *cp != '\0'; cp++);
    if (*cp == '#' || !strlen(cp)) continue;
    if (sscanf(line, "%f %f %f %f %f %f", 
	       &cyl.r, &cyl.z, &v, &eabs, &er, &ez) != 6){
      error("failed to read electric field data from line no %d\n"
	    "of file %s\n", lineno, fname);
      fclose(fp);
      return 1;
    }
#endif
    i = (int)rint((cyl.r - rmin)/rstep);
    j = (int)rint((cyl.z - zmin)/zstep);
    if (i < 0 || i >= rlen || j < 0 || j >= zlen){
      continue;
    }
    cyl.phi = 0;
    cart = cyl_to_cart(cyl);
    if (segment_number(cart) < 0) continue;
#ifdef OLD_FORMAT
    v[i][j] = vtmp;
    vexists[i][j] = 1;
#else
    efld[i][j].r = er;
    efld[i][j].z = ez;
    efld[i][j].phi = 0;
#endif
  }      

#ifdef OLD_FORMAT
  for (i = 0; i < rlen; i++){
    for (j = 0; j < zlen; j++){
      /*differentiate potential to get electric field*/
      if (!vexists[i][j]) continue;
      if (i - 1 >= 0 && i + 1 < rlen
	  && vexists[i-1][j] && vexists[i+1][j]){
	efld[i][j].r = -(v[i+1][j] - v[i-1][j])/2/rstep*10;//V/cm 
      }else if ((i - 1 < 0 || !vexists[i-1][j])
		&& i + 1 < rlen && vexists[i+1][j]){
	efld[i][j].r = -(v[i+1][j] - v[i][j])/rstep*10;//V/cm
      }else if ((i + 1 >= rlen || !vexists[i+1][j])
		&& i - 1 >= 0 && vexists[i-1][j]){
	efld[i][j].r = -(v[i][j] - v[i-1][j])/rstep*10;//V/cm
      }
      if (j - 1 >= 0 && j + 1 < zlen
	  && vexists[i][j-1] && vexists[i][j+1]){
	efld[i][j].z = -(v[i][j+1] - v[i][j-1])/2/zstep*10;//V/cm 
      }else if ((j - 1 < 0 || !vexists[i][j-1])
		&& j + 1 < zlen && vexists[i][j+1]){
	efld[i][j].z = -(v[i][j+1] - v[i][j])/zstep*10;//V/cm
      }else if ((j + 1 >= zlen || !vexists[i][j+1])
		&& j - 1 >= 0 && vexists[i][j-1]){
	efld[i][j].z = -(v[i][j] - v[i][j-1])/zstep*10;//V/cm
      }
      /*tell(ANNOYINGLY_CHATTY, "e[%d][%d] = %f %f\n", 
	i,j,efld[i][j].r,efld[i][j].z);*/
    }
  }
#endif

#ifdef OLD_FORMAT
  for (i = 0; i < rlen; i++){
    free(v[i]);
    free(vexists[i]);
  }
  free(v);
  free(vexists);
#endif

  tell(NORMAL, "Done reading %d lines of electric field data\n", lineno);
  fclose(fp);
  return 0;
}


/*setup_wp
  read weighting potential values from files. returns 0 on success*/
static int setup_wp(char fnames[20][MAX_FNAME_LEN]){
  int fno, segno;
  FILE *fp;
  char line[MAX_LINE];
  int i, j;
  int rlen, zlen;
  cyl_pt cyl;
  point cart;
  int lineno;
  float wp;
  char *cp;

  rlen = (int)rint((rmax - rmin)/rstep);
  zlen = (int)rint((zmax - zmin)/zstep);
  printf("rlen, zlen, nsegs = %d %d %d\n", rlen, zlen, nsegs);

  if (wpot == NULL){//assuming xlen, ylen, zlen never change as for setup_efld
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

  /*now read the tables*/
  for (fno = segno = 0; segno < nsegs; segno++, fno++){
    /* data in files is for segments 0,1,2, etc*/
    if ((fp = fopen(fnames[fno], "r")) == NULL){
      error("failed to open file: %s\n", fnames[fno]);
      return -1;
    }
    tell(NORMAL, "Reading weighting potential from file: %s\n",
	 fnames[fno]);
    lineno = 0;

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
      i = (int)rint((cyl.r - rmin)/rstep);
      j = (int)rint((cyl.z - zmin)/zstep);
      if (i < 0 || i >= rlen || j < 0 || j >= zlen) continue;
      cyl.phi = 0;
      cart = cyl_to_cart(cyl);
      //if (segment_number(cart) < 0) continue;
      //printf("%d %d %.2f %.2f %.4f \n", i, j, cyl.r, cyl.z, wp); fflush(stdout);
      wpot[i][j][segno] = wp;
    }
    tell(NORMAL,"Read %d lines from file %s\n", lineno, fnames[fno]);
    fclose(fp);
  }
  tell(NORMAL, "Done reading weighting potential\n");
  return 0;
}


/* free malloc()'ed memory and do other cleanup*/
int fields_finalize(void){
  int i, j;

  for (i = 0; i < (int) rint(rmax - rmin)/rstep; i++){
    free(efld[i]);
    for (j=0; j < (int) rint(zmax - zmin)/zstep; j++) free(wpot[i][j]);
    free(wpot[i]);
  }
  free(efld);
  free(wpot);
  free(v_lookup);
  efld = NULL;
  wpot = NULL;
  v_lookup = NULL;

  return 1;
}




static float tcorr(float q, float abse){
  float dt, res, vc, f;
  int i;

  for (i = 0; i < v_lookup_len - 2 && abse > v_lookup[i+1].e; i++);
  f = (abse - v_lookup[i].e)/(v_lookup[i+1].e - v_lookup[i].e); 
  if (q > 0){
    vc = (v_lookup[i+1].hcorr - v_lookup[i].hcorr)*f 
      + v_lookup[i].hcorr;
  }else{
    vc = (v_lookup[i+1].ecorr - v_lookup[i].ecorr)*f 
      + v_lookup[i].ecorr;
  }

  dt = temperature - REF_TEMP;
  res = 1+vc*dt/100;

  return res;
}

void set_temp(float temp){
  if (temp < MIN_TEMP || temp > MAX_TEMP){
    printf("temperature out of range: %f\n", temp);
  }else{
    temperature = temp;
    printf("temperature set to %f\n", temperature);
  }
}

void set_fudge_E_r(float f1, float f2){
  fudge_radial_field_1 = f1;
  fudge_radial_field_2 = f2;
}
