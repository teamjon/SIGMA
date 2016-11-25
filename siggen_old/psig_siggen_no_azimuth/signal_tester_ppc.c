/* signal_tester.c
 *  Karin Lagergren
 *
 * a simple program that communicates with the signal 
 * calculation code. Given coordinate => spectrum, print-out of signal,
 * or segment number. Given two sets of coords=>rms distance
 * "cart" and "cyl" switch between cartesian and cylindrical coordinates
 * for input
 *
 * to compile: 
 *  gcc -o st signal_tester.c point.c cyl_point.c calc_signal.c\
 *    fields.c detector_geometry.c signal_calc_util.c -lm -lreadline
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
 double trunc(double x);
#include <readline/readline.h>
#include <readline/history.h>
#include <ctype.h>

#include "calc_signal.h"
#include "cyl_point.h"
#include "signal_calc_util.h"
#include "detector_geometry.h"
#include "fields.h"

#define PROMPT ">> "
#ifdef PPC
#define GEOMETRY_FILE "geometry_setup_ppc.dat"
#define FIELDS_DIR_FILE "field_setup_ppc.dat"
#define SIGNAL_PARS_FILE "calc_signal_setup_ppc.dat"
#elif defined WEIRD_COAXIAL_GEOMETRY
#define GEOMETRY_FILE "geometry_setup_wc.dat"
#define FIELDS_DIR_FILE "field_setup_wc.dat"
#define SIGNAL_PARS_FILE "calc_signal_setup_wc.dat"
#else
#define GEOMETRY_FILE "geometry_setup.dat"
#define FIELDS_DIR_FILE "field_setup.dat"
#define SIGNAL_PARS_FILE "calc_signal_setup.dat"
#endif

#define RZ_STEP 0.5
#define XYZ_STEP 0.5
#define P_STEP (M_PI/100)

#define MAX_LINE 512

#define MAX_SPE_CHS 8192

#define FREE(p) do{free(p); p = NULL; }while(0)

#define CYL 0
#define CART 1

#define C2_SIGNAL_MULT_F 10000
#define C2_THRESH 100

#define NON_LINEAR_CHISQ 1

static int coord_type = CYL;

static int rl_gets(char *line, int MAX_LEN);
static int write_spectrum(float *spec, int nchs, char *name);
#ifndef PPC
static int get_neighbors(int seg, int *neigh);
#endif

static float fit_time(int nsigs, int ntimes, float **s1, float **s2, 
		      int shift_s1);
static int rc_integrate(float **s_in, float **s_out, float tau);
int shift_signal(float **s, float dt, int nsigs, int ntimes);

static int set_cyl(char *cmd);
static int set_cart(char *cmd);
static int which_segment(char *cmd);
static int save_signal(char *cmd);
static int az_drift_time(char *cmd);
static int print_signal(char *cmd);
static int rms_dist(char *cmd);
static int geometrical_seg(char *cmd);
static int print_help(char *cmd);
static int drift_paths(char *cmd);
static int c2_dist(char *cmd);
static int fit_t0(char *cmd);
static int sett(char *cmd);
static int set_tau(char *cmd);
static int set_fEr(char *cmd);
void set_fudge_E_r(float f1, float f2); // FIXME

static int time_steps;
static int nsegments;
static float tau, tau7;

static struct{
  char cmd[MAX_LINE];
  int (*f)(char *);
  char help_str[MAX_LINE];
}cmds[] = {{"cyl", set_cyl, "cyl"},
           {"cart", set_cart, "cart"},
           {"seg", which_segment, "seg x y z or seg r p z"},
           {"sig", save_signal, "sig x y z file.spe or sig r p z file.spe"},
	   {"psig", print_signal, "psig x y z or psig r p z"},
           {"dt", az_drift_time, "dt r z fn.dat"},
	   {"rmsd", rms_dist, "rmsd x1 y1 z1 x2 y2 z2 or rmsd r1 p1 z2 r2 p2 z2"},
	   {"c2d", c2_dist, "c2d fix_coord val out_file.dat [rc_int tau]"},	   
	   {"ft0", fit_t0, "ft0 x1 y1 z1 x2 y2 z2 fn or ft0 r1 p1 z2 r2 p2 z2 fn"},	   
	   {"geo", geometrical_seg, "geo x y z or geo r p z"},
	   {"dp", drift_paths, "dp fn.dat"},
	   {"st", sett, "st %f"},
	   {"tau", set_tau, "tau %f"}, 
	   {"fer", set_fEr, "fer %f %f"}, 
           {"help", print_help, "help"}};

int main(void){
  char ans[MAX_LINE], *cp;
  int i, ncmds;
  
#ifdef CHATTY
  set_signal_calc_output(ANNOYINGLY_CHATTY, vprintf);
#endif
  if (signal_calc_init(GEOMETRY_FILE, FIELDS_DIR_FILE,
		       SIGNAL_PARS_FILE, &time_steps, &nsegments) != 0){
    return 1;
  }
  
  ncmds = sizeof(cmds)/sizeof(cmds[0]);
  for ( ; ; ){
    rl_gets(ans, MAX_LINE);
    for (cp = ans; isspace(*cp); cp++);
    if (strlen(cp) == 0) continue;
    if (strncmp(cp, "quit", 4) == 0) return 0;
    for (i = 0; i < ncmds; i++){
      if (strncmp(cp, cmds[i].cmd, strlen(cmds[i].cmd)) == 0){
	cp += strlen(cmds[i].cmd);
	cmds[i].f(cp);
	break;
      }
    } 
    if (i == ncmds){
      printf("unknown command: %s\n", ans);
    }
  }
  return 0;
}

/*readline wrapper*/
static int rl_gets(char *line, int MAX_LEN){
  static char *line_read = (char *) NULL;

  if (line_read != NULL){
    FREE(line_read);
    line_read = (char *) NULL;
  }

  line_read = readline(PROMPT);
  if (line_read != NULL && *line_read != '\0'){
    add_history(line_read);
  }
  strncpy(line,line_read,MAX_LEN);
  line[MAX_LEN -1] = '\0';
  FREE(line_read);

  return 1;
}

#ifndef PPC
static int get_neighbors(int seg, int *neigh){
  int i, n;

  neigh[0] = (seg % 6) == 0 ? seg + 5 : seg - 1;
  neigh[1] = (seg % 5) == 0 ? seg - 5 : seg + 1;
  neigh[2] = (seg < 6)      ? -1      : seg - 6;
  neigh[3] = (seg >= 30)    ? -1      : seg + 6;
  
  for (n = 0, i = 0; i < 4; i++)
    if (neigh[i] >=0) n++;
  return n;
}
#endif

struct spe_header{
  int reclA;            /* 24 */
  unsigned int title[2]; /*don't know why this is uint, but seems to
                           work, so...*/ 
  int dim;
  int a1;               /*  1 */
  int a2;               /*  1 */
  int a3;               /*  1 */
  int reclB;            /* 24 */
};

/*write_spectrum
 *
 * saves the spectrum pointed to by "spec". Returns 0 if unsuccessful,
 * 1 if successful.
 * The name is truncated by removing any trailing .spe (and any
 * subsequent characters...), as well as any leading directory names.
 * If the resulting string is longer than the maximum allowed 8
 * characters, only the first 8 are retained
 */
static int write_spectrum(float *spec, int nchs, char *name){  
  FILE *fp;
  int record_length;
  struct spe_header header;
  char *suffix_start;
  char *fname_start;

  header.reclA = header.reclB = 24; 
  header.title[0] = header.title[1] = 0;
  header.a1 = header.a2 = header.a3 = 1;

  fp = fopen(name,"w");
  if (fp == NULL){
    fprintf(stderr,"Error! Unable to open spectrum file %s \n",name);
    return 0;
  }
  header.dim = nchs;
  if ((suffix_start = strstr(name,".spe")) == NULL){
    suffix_start = name + strlen(name);
  }
  if ((fname_start = rindex(name,'/')) == NULL){
    fname_start = name;
  }else{
    fname_start++;/*get rid of the '/'*/
  }
  if (suffix_start - fname_start < 8){
    memcpy(header.title,"       ",8);/*blank the title*/
    memcpy(header.title,fname_start,suffix_start - fname_start);
  }else{ 
    memcpy(header.title,suffix_start - 8,8);
  }
  record_length = sizeof(float)*header.dim;

  fwrite(&header, sizeof(header), 1, fp);
  fwrite(&record_length, sizeof(record_length), 1, fp);
  fwrite(spec, sizeof(float), nchs, fp); 
  fwrite(&record_length, sizeof(record_length), 1,fp);
  fclose(fp);

  return 1;
}


static int set_cyl(char *cmd){
  coord_type = CYL;
  printf("coordinate system: cylindrical\n");
  
  return 0;
}

static int set_cart(char *cmd){
  coord_type = CART;
  printf("coordinate system: cartesian\n");

  return 0;
}

/*prints segment number (-1 for outside crystal)*/
static int which_segment(char *cmd){
  int seg;
  struct point cart;
  struct cyl_pt cyl;

  if (coord_type == CYL){
    if (sscanf(cmd, "%f %f %f", &cyl.r, &cyl.phi, &cyl.z) == 3
	|| sscanf(cmd, "%f, %f, %f", &cyl.r, &cyl.phi, &cyl.z) == 3){
      cart = cyl_to_cart(cyl);
    }else{
      printf("error parsing coordinate: %s\n", cmd);
      return 1;
    }
  }else{
    if (sscanf(cmd, "%f %f %f", &cart.x, &cart.y, &cart.z) == 3
	  || sscanf(cmd, "%f, %f, %f", &cart.x, &cart.y, &cart.z) == 3){
      ;//nothing
    }else{
      printf("error parsing coordinate: %s\n", cmd);
      return 1;
    }
  }
  seg = hit_segment(cart);
  printf("hit segment for (x = %.1f, y = %.1f, z = %.1f): %d\n", 
	 cart.x, cart.y, cart.z,seg);

  return 0;
}

static int save_signal(char *cmd){
  float spec[MAX_SPE_CHS];
  struct point cart;
  struct cyl_pt cyl;
  int seg;
  int comp_f;
  int i, j, k, t50;
  char *cp, *cp2;
  static float **s;

  if (s == NULL){//first call
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
  }

  if (coord_type == CYL){
    if (sscanf(cmd, "%f %f %f", &cyl.r, &cyl.phi, &cyl.z) == 3
	  || sscanf(cmd, "%f, %f, %f", &cyl.r, &cyl.phi, &cyl.z) == 3){
      cart = cyl_to_cart(cyl);
    }else{
      printf("error parsing coordinate: %s\n", cmd);
      return 1;
    }
  }else{
    if (sscanf(cmd, "%f %f %f", &cart.x, &cart.y, &cart.z) == 3
	  || sscanf(cmd, "%f, %f, %f", &cart.x, &cart.y, &cart.z) == 3){
      ;//nothing
    }else{
      printf("error parsing coordinate: %s\n", cmd);
      return 1;
    }
  }
  
  for (cp = cmd; isspace(*cp); cp++);
  for ( ; isdigit(*cp) || *cp == '.'; cp++);//skip r/x coord.
  for ( ; isspace(*cp) || *cp == ','; cp++); //skip whitespace, commas
  for ( ; isdigit(*cp) || *cp == '.'; cp++);//skip phi/y coord.
  for ( ; isspace(*cp) || *cp == ','; cp++);
  for ( ; isdigit(*cp) || *cp == '.'; cp++);//skip z coord

  for ( ; isspace(*cp); cp++); //skip whitespace
  for (cp2 = cp + strlen(cp); isspace(*cp2); cp2--)//remove whitespace
    *cp2 = '\0';
  if (strlen(cp) == 0){
    fprintf(stderr, "must supply file name\n");
    return 1;
  }
  seg = get_signal(cart, s);

  if (seg < 0) {
    printf("point not in crystal: (x = %.1f, y = %.1f, z = %.1f)\n", 
	   cart.x, cart.y, cart.z);
    return 1;
  }
  // FIXME
  // for (comp_f = 1; nsegments*time_steps/comp_f > MAX_SPE_CHS; comp_f *= 2) ;
  for (comp_f = 1; time_steps/comp_f > 400; comp_f *= 2) ;

  rc_integrate(s, s, tau);

  printf("spectrum will be compressed by factor %d\n", comp_f);
  for (i = 0; i < sizeof(spec)/sizeof(spec[0]); i++) spec[i] = 0;

  /* find 50% time of PC (seg 0) */
  for (t50 = 0; t50 < time_steps && s[0][t50]/s[0][time_steps-1] < 0.5f; t50++) ;
  tell(ANNOYINGLY_CHATTY, "t50 = %d; %f\n", t50, s[0][time_steps-2]);
  t50 -= 260; // align to ch 260, so shift left by new value of t50

  /*copy signal data to spectrum array*/
  for (i = 0; i < nsegments; i++){
    k = i*400;
    if (i>0) k += 2800;  // FIXME
    if (0) {
      /* old version, without time-alignment to PC t50 */
      for (j = 0; j < time_steps; j++){
	//spec[(i*time_steps+j)/comp_f] += s[i][j]*1000/comp_f;
	spec[(k+j)/comp_f] += s[i][j]*1000/comp_f;
      }
    } else {
      k -= t50;
      if (t50 > 0) {  // t50 positive; shift left
	for (j = t50; j < time_steps; j++){
	  spec[(k+j)/comp_f] += s[i][j]*1000/comp_f;
	}
      } else {  // t50 negative; shift right
	for (j = 0; j < time_steps && j - t50 < 400; j++){
	  spec[(k+j)/comp_f] += s[i][j]*1000/comp_f;
	}
      }
    }
  }
  /* spread seg 8 over segs 1-7 */  // FIXME
  /* multiply seg 19 by 4 */  // FIXME
  for (i = 0; i < 400; i++) {
    for (j=1; j<8; j++) {
      spec[i + j*400] = spec[i+3200]/8.0;
    }
    spec[i + 19*400] = spec[i + 19*400] * 4.0f;
  }
  /* fudge relative gain */  // FIXME
  float seg_chs[] = {1877, 1462, 1476, 1443, 1403, 1449, 1456, 1502, 1426, 1319,
		     1345, 1436, 1336, 1427, 1485, 1450, 1463, 1427, 1304, 247};
  float egam = 1332.5f;
  for (i = 0; i < 400; i++) {
    for (j=0; j<nsegments + 7; j++) {  // FIXME
      spec[i + j*400] *= -5.006f/1.409f * seg_chs[j]/egam;
    }
  }

  // write_spectrum(spec, nsegments*time_steps/comp_f, cp);
  // printf("%d channels of data saved in spectrum %s\n", nsegments*time_steps/comp_f, cp);
  write_spectrum(spec, MAX_SPE_CHS, cp);
  printf("%d channels of data saved in spectrum %s\n", MAX_SPE_CHS, cp);
  return 0;
}

static int print_signal(char *cmd){
  static float **s;
  struct point cart;
  struct cyl_pt cyl;
#ifndef PPC
  int j, neigh[4];
#endif
  int i, seg;

  if (s == NULL){//first call
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
  }

  if (coord_type == CYL){
    if (sscanf(cmd, "%f %f %f", &cyl.r, &cyl.phi, &cyl.z) == 3
	|| sscanf(cmd, "%f, %f, %f", &cyl.r, &cyl.phi, &cyl.z) == 3){
      cart = cyl_to_cart(cyl);
    }else{
      printf("error parsing coordinate: %s\n", cmd);
      return 1;
    }
  }else{
    if (sscanf(cmd, "%f %f %f", &cart.x, &cart.y, &cart.z) == 3
	|| sscanf(cmd, "%f, %f, %f", &cart.x, &cart.y, &cart.z) == 3){
      ;//nothing
    }else{
      printf("error parsing coordinate: %s\n", cmd);
      return 1;
    }
  }

  seg = get_signal(cart, s);

  if (seg < 0) {
    printf("point not in crystal: (x = %.1f, y = %.1f, z = %.1f)\n", 
	   cart.x, cart.y, cart.z);
    return 1;
  }

  printf("signal is for hit segment (no %d)\n", seg);
  for (i = 0; i < time_steps; i++){
    printf("%.3f ", s[seg][i]);
    if ((i+1)%10 == 0)
      printf("\n");
  }
  
#ifndef PPC
  get_neighbors(seg, neigh);
  for (j = 0; j < 4; j++){
    if (neigh[j] < 0) continue;
    printf("for neighboring segment %d\n",neigh[j]);
    for (i = 0; i < time_steps; i++){
      printf("%.3f ", s[neigh[j]][i]);
      if ((i+1)%10 == 0)
	printf("\n");
    }
  }
  printf("for central contact: \n");
  for (i = 0; i < time_steps; i++){
    printf("%.3f ", s[nsegments-1][i]);
    if ((i+1)%10 == 0)
      printf("\n");
  }
#endif
  return 0;
}

static int az_drift_time(char *cmd){
  FILE *fp;
  struct point cart;
  struct cyl_pt cyl;
  int dt;
  int i, phi;
  char *cp, *cp2;
  static float **s;

  if (s == NULL){//first call
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
  }

  if (sscanf(cmd, "%f %f", &cyl.r, &cyl.z) == 2 ||
      sscanf(cmd, "%f,%f", &cyl.r, &cyl.z) == 2) {
    cyl.phi = 0.0f;
    cart = cyl_to_cart(cyl);
  } else {
    printf("error parsing coordinate: %s\n", cmd);
    return 1;
  }
  
  for (cp = cmd; isspace(*cp); cp++);
  for ( ; isdigit(*cp) || *cp == '.'; cp++); //skip r coord.
  for ( ; isspace(*cp) || *cp == ','; cp++); //skip whitespace, commas
  for ( ; isdigit(*cp) || *cp == '.'; cp++); //skip z coord

  for ( ; isspace(*cp); cp++); //skip whitespace
  for (cp2 = cp + strlen(cp); isspace(*cp2); cp2--)//remove whitespace
    *cp2 = '\0';
  if (strlen(cp) == 0){
    fprintf(stderr, "must supply file name\n");
    return 1;
  }
  if ((fp = fopen(cp, "w")) == NULL){
    printf("failed to open output file: %s\n", cp);
    return 1;
  }
  fprintf(fp,
	  " 1 2\n"
	  "# r=%.2f z=%.2f\n"
	  "# ang(deg) drift_time(ps)\n", cyl.r, cyl.z);

  for (phi=0; phi<100; phi++) {
    cyl.phi = 3.141593 * (float) phi / 180.0;
    cart = cyl_to_cart(cyl);
    dt = get_signal(cart, s);
    if (dt < 0) {
      printf("point not in crystal: (r = %.1f, phi = %.1f, z = %.1f)\n", 
	     cyl.r, cyl.phi, cyl.z);
      fclose(fp);
      return 1;
    }
    printf("ang, dt = %3d %d\n", phi, dt);
    fprintf(fp, "%3d %d\n", phi, dt);
  }
  printf("Data saved in file %s\n", cp);
  fclose(fp);
  return 0;
}

static int rms_dist(char *cmd){
  struct point cart1, cart2;
  struct cyl_pt cyl1, cyl2;
  int seg1, seg2;
  static float **s1, **s2;
  int i,j;
  float d;
  float rmsd;

  if (s1 == NULL){
    if ((s1 = malloc(nsegments*sizeof(*s1))) == NULL
      || (s2 = malloc(nsegments*sizeof(*s2))) == NULL){
      printf("malloc failed\n");
    }
    for (i = 0; i < nsegments; i++){
      if ((s1[i] = malloc(time_steps*sizeof(*s1[i]))) == NULL
	|| (s2[i] = malloc(time_steps*sizeof(*s2[i]))) == NULL){
	printf("malloc failed\n");
      }

    }
  }

  if (coord_type == CYL){
    if (sscanf(cmd, "%f %f %f %f %f %f", 
	       &cyl1.r, &cyl1.phi, &cyl1.z,
	       &cyl2.r, &cyl2.phi, &cyl2.z) == 6
	||sscanf(cmd, "%f, %f, %f, %f, %f, %f", 
		 &cyl1.r, &cyl1.phi, &cyl1.z,
		 &cyl2.r, &cyl2.phi, &cyl2.z) == 6){
      
      cart1 = cyl_to_cart(cyl1);
      cart2 = cyl_to_cart(cyl2);
    }else{
      fprintf(stderr,"Failed to parse coordinates: %s\n", cmd);
      return 1;
    }
  }else{
    if (sscanf(cmd, "%f %f %f %f %f %f", 
	       &cart1.x, &cart1.y, &cart1.z,
	       &cart2.x, &cart2.y, &cart2.z) == 6
	||sscanf(cmd, "%f, %f, %f, %f, %f, %f", 
		 &cart1.x, &cart1.y, &cart1.z,
		 &cart2.x, &cart2.y, &cart2.z) == 6){
      ;//nothing
    }else{
      fprintf(stderr,"Failed to parse coordinates: %s\n", cmd);
      return 1;
    }
  }

  printf("finding signal for (x = %.1f, y = %.1f, z = %.1f)\n",
	 cart1.x, cart1.y, cart1.z);      
  seg1 = hit_segment(cart1);
  if (seg1 < 0) {
    printf("not a valid coordinate\n");
    return 1;
  }
  printf("signal 1 is for segment no %d\n", seg1);

  printf("finding signal for (x = %.1f, y = %.1f, z = %.1f)\n",
	 cart2.x, cart2.y, cart2.z);      
  seg2 = hit_segment(cart2);
  if (seg2 < 0) {
    printf("not a valid coordinate\n");
    return 1;
  }
  printf("signal 2 is for segment no %d\n", seg2);

  seg1 = get_signal(cart1, s1);
  seg2 = get_signal(cart2, s2);

  if (seg1 != seg2){
    printf("NB: points are in different segments\n");
  }

  rmsd = 0.0;
  for (i = 0; i < nsegments; i++){
    for (j = 0; j < time_steps; j++){
      d = s1[i][j] - s2[i][j];
      rmsd +=d*d;
    }
  }
  printf("RMS difference: %f\n", rmsd);
  
  return 0;
}

static int geometrical_seg(char *cmd){
  int seg;
  struct point cart;
  struct cyl_pt cyl;

  if (coord_type == CYL){
    if (sscanf(cmd, "%f %f %f", &cyl.r, &cyl.phi, &cyl.z) == 3
	|| sscanf(cmd, "%f, %f, %f", &cyl.r, &cyl.phi, &cyl.z) == 3){
      cart = cyl_to_cart(cyl);
    }else{
      printf("error parsing coordinate: %s\n", cmd);
      return 1;
    }
  }else{
    if (sscanf(cmd, "%f %f %f", &cart.x, &cart.y, &cart.z) == 3
	  || sscanf(cmd, "%f, %f, %f", &cart.x, &cart.y, &cart.z) == 3){
      ;//nothing
    }else{
      printf("error parsing coordinate: %s\n", cmd);
      return 1;
    }
  }
  seg = segment_number(cart);
  printf("geom. segment for (x = %.1f, y = %.1f, z = %.1f): %d\n", 
	 cart.x, cart.y, cart.z,seg);

  return 0;

}

static int sett(char *cmd){
  float t;
  char *endp;
  
  t = strtod(cmd, &endp);
  if (endp == cmd){
    printf("cannot parse temperature: %s\n", cmd);
    return 1;
  }
  set_temp(t);
  return 0;
}

static int set_tau(char *cmd){
  float t;
  char *endp, *e2;
  
  t = strtod(cmd, &endp);
  if (endp == cmd){
    printf("cannot parse tau: %s\n", cmd);
    return 1;
  }
  if (t >= 0) tau = t;
  tau7 = tau;
  t = strtod(endp, &e2);
  if (e2 != endp && t > 0.1) tau7 = t;
  tell(NORMAL, "Signals will be integrated with tau = %f %f\n", tau, tau7);

  return 0;
}

static int set_fEr(char *cmd){
  float f1, f2;
  char *endp, *e2;
  
  f1 = strtod(cmd, &endp);
  if (endp == cmd){
    printf("cannot parse f1: %s\n", cmd);
    f1 = f2 = 0;
    return 1;
  }
  f2 = strtod(endp, &e2);
  if (e2 == endp && fabs(f1) > 0.001){
    printf("cannot parse f2: %s\n", cmd);
    f1 = f2 = 0;
    return 1;
  }
  tell(NORMAL, "Fudging radial field with %f, %f\n", f1, f2);
  set_fudge_E_r(f1, f2);
  return 0;
}

static int print_help(char *cmd){
  int i;

  printf("available commands:\n");
  for (i = 0; i < sizeof(cmds)/sizeof(cmds[0]); i++){
    printf("%s\n",cmds[i].help_str);
  }
  return 0;
}


static int drift_paths(char *cmd){
  FILE *fp;
  point *dpe, *dph;
  int i, nt;
  char *cp, *cp2;

  for (cp = cmd; isspace(*cp); cp++);
  for (cp2 = cp + strlen(cp); isspace(*cp2); cp2--)//remove whitespace
    *cp2 = '\0';
  if (strlen(cp) == 0){
    fprintf(stderr, "must supply file name\n");
    return 1;
  }

  if ((fp = fopen(cp, "w")) == NULL){
    printf("failed to open output file: %s\n", cp);
    return 1;
  }
  
  nt = drift_path_e(&dpe);
  drift_path_h(&dph);
  fprintf(fp,"n ex ey ez hx hy hz\n");
  for (i = 0; i < nt; i++){
    fprintf(fp,"%d %.3f %.3f %.3f %.3f %.3f %.3f\n",
	   i, dpe[i].x, dpe[i].y, dpe[i].z, dph[i].x, dph[i].y, dph[i].z);
  }

  fclose(fp);
  return 0;
}

static int c2d_xy(float zval, FILE *fp, float rcint_tau);
static int c2d_xz(float yval, FILE *fp, float rcint_tau);
static int c2d_yz(float xval, FILE *fp, float rcint_tau);
static int c2d_rp(float zval, FILE *fp, float rcint_tau);
static int c2d_rz(float pval, FILE *fp, float rcint_tau);
static int c2d_pz(float rval, FILE *fp, float rcint_tau);


static int c2_dist(char *cmd){
  char dir;
  float fixval;
  char *cp, *cp2;
  FILE *fp;
  float rmax, zmax;
  float rcint_tau = -1;
  
  for (cp = cmd; isspace(*cp); cp++);
  dir = *cp;
  cp ++;
  if (sscanf(cp, " %f", &fixval) != 1){
    printf("failed to parse fix value for direction %c\n", dir);
    return -1;
  }
  rmax = rmax_detector();
  zmax = zmax_detector();
  if (dir == 'z' && (fixval < 0 || fixval > zmax)){
    printf(" z = %f is outside detector\n", fixval);
    return -1;
  }
  if (dir == 'r' && (fixval < 0 || fixval > rmax)){
    printf(" r = %f is outside detector\n", fixval);
    return -1;
  }
  if (dir == 'x' && (fixval < -rmax/2 || fixval > rmax/2)){
    printf(" x = %f is outside detector\n", fixval);
    return -1;
  }else if (dir == 'y' && (fixval < -rmax/2 || fixval > rmax/2)){
    printf(" y = %f is outside detector\n", fixval);
    return -1;
  }   
  while(isspace(*cp)) cp++;
  while(!isspace(*cp)) cp++;

  for (; isspace(*cp); cp++)
    ;//nada
  for (cp2 = cp; !isspace(*cp2) && *cp2 != '\0'; cp2++)
    ;//nada
  if (sscanf(cp2, " %f", &rcint_tau) == 1){
    if (rcint_tau > 0){
      printf("signals will be integrated with tau = %f\n", rcint_tau);
    }else{
      printf("invalid tau: %f\n", rcint_tau);
    }
  }
  *cp2 = '\0';
  if ((fp = fopen(cp, "w")) == NULL){
    printf("failed to open output file: %s\n", cp);
    return -1;
  }

  if (dir == 'r')
    c2d_pz(fixval, fp, rcint_tau);
  else if (dir == 'p')
    c2d_rz(fixval, fp, rcint_tau);
  else if (dir == 'z' && coord_type == CYL)
    c2d_rp(fixval, fp, rcint_tau);
  else if (dir == 'x')
    c2d_yz(fixval, fp, rcint_tau);
  else if (dir == 'y')
    c2d_xz(fixval, fp, rcint_tau);
  else if (dir == 'z')
    c2d_xy(fixval, fp, rcint_tau);

  fclose(fp);
  return 0;
}

#define TSTEPS_TOT (time_steps + 10)
static int c2d_xy(float zval, FILE *fp, float rcint_tau){
  printf("not implemented\n");
  return 0;
}
static int c2d_xz(float yval, FILE *fp, float rcint_tau){
  printf("not implemented\n");
  return 0;
}
static int c2d_yz(float xval, FILE *fp, float rcint_tau){
  printf("not implemented\n");
  return 0;
}



static int c2d_rp(float zval, FILE *fp, float rcint_tau){
  float ***s, ***s2;
  float **fr, **fphi, f, **fr2, **fphi2, f2;
  int nr, np;
  float rmax, pmax,  d;
  int i,j, k, l;
  cyl_pt cyl[2];
  point cart[2];
  int seg[2];
  int n1, n2;

  rmax = rmax_detector();
  pmax = 2*M_PI;
  nr = rmax/RZ_STEP + 1;
  np = pmax/P_STEP + 2;//?
  
  if ((s = malloc(2*sizeof(*s))) == NULL
      || (s2 = malloc(2*sizeof(*s2))) == NULL){
    printf("malloc failed\n");
    exit(1);
  }
  for (i = 0; i < 2; i++){
    if ((s[i] = malloc(nsegments*sizeof(*s[i]))) == NULL
	|| (s2[i] = malloc(nsegments*sizeof(*s2[i]))) == NULL){
      printf("malloc failed\n");
      exit(1);
    }    
    for (j = 0; j < nsegments; j++){
      if ((s2[i][j] = malloc(TSTEPS_TOT*sizeof(*s[i][j]))) == NULL){
	printf("malloc failed\n");
	exit(1);
      }    
      s[i][j] = s2[i][j] + (TSTEPS_TOT - time_steps);
      memset(s2[i][j], 0, sizeof(*s2[i][j])*TSTEPS_TOT);
    }
  }
  if ((fr = malloc(nr*sizeof(*fr))) == NULL
      || (fphi = malloc(nr*sizeof(*fphi))) == NULL
      || (fr2 = malloc(nr*sizeof(*fr2))) == NULL
      || (fphi2 = malloc(nr*sizeof(*fphi2))) == NULL){
    printf("malloc failed\n");
    exit(1);
  }
  for (i = 0; i < nr; i++){
    if ((fr[i] = malloc(np*sizeof(*fr[i]))) == NULL
	|| (fphi[i] = malloc(np*sizeof(*fphi[i]))) == NULL
	|| (fr2[i] = malloc(np*sizeof(*fr2[i]))) == NULL
	|| (fphi2[i] = malloc(np*sizeof(*fphi2[i]))) == NULL){
      printf("malloc failed\n");
      exit(1);
    }
    for (j = 0; j < np; j++){
      fr[i][j] = fphi[i][j] = fr2[i][j] = fphi2[i][j] = -2;
    }
  }
  n1 = 0; n2 = 1;
  cyl[n1].z = cyl[n2].z = zval;
  for (j = 0; j < np; j++){
    printf("\r%d/%d", j+1, np);fflush(stdout);
    cyl[n1].r = 0;
    cyl[n1].phi = cyl[n2].phi = 0 + j*P_STEP;
    cart[n1] = cyl_to_cart(cyl[n1]);
    seg[n1] = get_signal(cart[n1], s[n1]);
    if (rcint_tau > 0 && seg[n1] >= 0) rc_integrate(s[n1], s[n1], rcint_tau);
    for (i = 1; i < nr; i++){
      cyl[n2].r = 0 + i*RZ_STEP;
      cart[n2] = cyl_to_cart(cyl[n2]);
      seg[n2] = get_signal(cart[n2], s[n2]);
      if (rcint_tau > 0 && seg[n2] >= 0) rc_integrate(s[n2], s[n2], rcint_tau);
      if (seg[n2] < 0 || seg[n1] < 0){
	n1 = n2;
	n2 = !n1;
	continue;
      }
      f2 = fit_time(nsegments, TSTEPS_TOT, s2[n1], s2[n2], 0);
      if (f2 > C2_THRESH){
	f2 = 2.354*RZ_STEP*sqrt(C2_SIGNAL_MULT_F)/sqrt(f2);
	fr2[i-1][j] = f2;
      }else{
	fr2[i-1][j] = -1;
      }
      f = 0;
      for (k = 0; k < nsegments; k++){
	for (l = 0; l < time_steps; l++){
	  d = C2_SIGNAL_MULT_F*(s[n1][k][l] - s[n2][k][l]);
	  f += d*d;
	}
      }
      if (f > C2_THRESH){
	f = 2.354*RZ_STEP*sqrt(C2_SIGNAL_MULT_F)/sqrt(f);	
	fr[i-1][j] = f;
      }else{
	fr[i-1][j] = -1;
      }
      n1 = n2;
      n2 = !n1;
    }
  }
  printf("\n");
  for (i = 0; i < nr; i++){
    printf("\r%d/%d", i+1, nr);fflush(stdout);
    cyl[n1].r = cyl[n2].r = 0 + i*RZ_STEP;
    cyl[n1].phi = 0;
    cart[n1] = cyl_to_cart(cyl[n1]);
    seg[n1] = get_signal(cart[n1], s[n1]);
    if (rcint_tau > 0 && seg[n1] >= 0) rc_integrate(s[n1], s[n1], rcint_tau);
    for (j = 1; j < np; j++){
      cyl[n2].phi = 0 + j*P_STEP;
      cart[n2] = cyl_to_cart(cyl[n2]);
      seg[n2] = get_signal(cart[n2], s[n2]);
      if (rcint_tau > 0 && seg[n2]>= 0) rc_integrate(s[n2], s[n2], rcint_tau);
      if (seg[n2] < 0 || seg[n1] < 0){
	n1 = n2;
	n2 = !n1;
	continue;
      }
      f2 = fit_time(nsegments, TSTEPS_TOT, s2[n1],s2[n2], 0);
      if (f2 > C2_THRESH){
	f2 = 2.354*P_STEP*cyl[n1].r*sqrt(C2_SIGNAL_MULT_F)/sqrt(f2);
	fphi2[i][j-1] = f2;
      }else{
	fphi2[i][j-1] = -1;
      }
      f = 0;
      for (k = 0; k < nsegments; k++){
	for (l = 0; l < time_steps; l++){
	  d = C2_SIGNAL_MULT_F*(s[n1][k][l] - s[n2][k][l]);
	  f += d*d;
	}
      }
      if (f > C2_THRESH){
	f = 2.354*P_STEP*cyl[n1].r*sqrt(C2_SIGNAL_MULT_F)/sqrt(f);	
	fphi[i][j-1] = f;
      }else{
	fphi[i][j-1] = -1;
      }
      n1 = n2;
      n2 = !n1;
    }
  }
  printf("\n");
  fprintf(fp,"r p fr fp fr2 fp2--- z = %f\n", zval);
  for (i = 0; i < nr; i++){
    cyl[n1].r = RZ_STEP/2 + i*RZ_STEP;
    for (j = 0; j < np; j++){
      cyl[n1].phi = P_STEP/2 + j*P_STEP;
      fprintf(fp, "%.2f %.2f %f %f %f %f\n", cyl[n1].r, cyl[n1].phi, fr[i][j], fphi[i][j], fr2[i][j], fphi2[i][j]);
    }
    fprintf(fp, "\n");
  }
  for (i = 0; i < 2; i++){
    for (j = 0; j < nsegments; j++){
      free(s2[i][j]);    
    }
    free(s[i]);
    free(s2[i]);
  }
  free(s);
  free(s2);

  for (i = 0; i < nr; i++){
    free(fr[i]);
    free(fphi[i]);
    free(fr2[i]);
    free(fphi2[i]);
  }
  free(fr);
  free(fphi);
  free(fr2);
  free(fphi2);
  return 0;
}

static int c2d_rz(float pval, FILE *fp, float rcint_tau){
  float ***s, ***s2;
  float **fr, **fz, f, **fr2, **fz2, f2;
  int nr, nz;
  float rmax, zmax,  d;
  int i,j, k, l;
  cyl_pt cyl[2];
  point cart[2];
  int seg[2];
  int n1, n2;


  rmax = rmax_detector();
  zmax = zmax_detector();
  nr = rmax/RZ_STEP + 1;
  nz = zmax/RZ_STEP + 1;
  
  if ((s = malloc(2*sizeof(*s))) == NULL
      ||(s2 = malloc(2*sizeof(*s2))) == NULL){
    printf("malloc failed\n");
    exit(1);
  }
  for (i = 0; i < 2; i++){
    if ((s[i] = malloc(nsegments*sizeof(*s[i]))) == NULL
	|| (s2[i] = malloc(nsegments*sizeof(*s2[i]))) == NULL){
      printf("malloc failed\n");
      exit(1);
    }    
    for (j = 0; j < nsegments; j++){
      if ((s2[i][j] = malloc(TSTEPS_TOT *sizeof(*s[i][j]))) == NULL){
	printf("malloc failed\n");
	exit(1);
      }    
      s[i][j] = s2[i][j] + (TSTEPS_TOT - time_steps);
      memset(s2[i][j], 0, sizeof(*s2[i][j])*TSTEPS_TOT);
    }
  }
  if ((fr = malloc(nr*sizeof(*fr))) == NULL
      || (fz = malloc(nr*sizeof(*fz))) == NULL
      || (fr2 = malloc(nr*sizeof(*fr2))) == NULL
      || (fz2 = malloc(nr*sizeof(*fz2))) == NULL){
    printf("malloc failed\n");
    exit(1);
  }
  for (i = 0; i < nr; i++){
    if ((fr[i] = malloc(nz*sizeof(*fr[i]))) == NULL
	|| (fz[i] = malloc(nz*sizeof(*fz[i]))) == NULL
	|| (fr2[i] = malloc(nz*sizeof(*fr2[i]))) == NULL
	|| (fz2[i] = malloc(nz*sizeof(*fz2[i]))) == NULL){
      printf("malloc failed\n");
      exit(1);
    }
    for (j = 0; j < nz; j++){
      fr[i][j] = fz[i][j] = fr2[i][j] = fz2[i][j] = -2;
    }
  }
  n1 = 0; n2 = 1;
  cyl[n1].phi = cyl[n2].phi = pval;
  for (j = 0; j < nz; j++){
    printf("\r%d/%d", j+1, nz);fflush(stdout);
    cyl[n1].r = 0;
    cyl[n1].z = cyl[n2].z = 0 + j*RZ_STEP;
    cart[n1] = cyl_to_cart(cyl[n1]);
    seg[n1] = get_signal(cart[n1], s[n1]);
    if (rcint_tau > 0 && seg[n1] >= 0) rc_integrate(s[n1], s[n1], rcint_tau);
    for (i = 1; i < nr; i++){
      cyl[n2].r = 0 + i*RZ_STEP;
      cart[n2] = cyl_to_cart(cyl[n2]);
      seg[n2] = get_signal(cart[n2], s[n2]);
      if (rcint_tau > 0 && seg[n2] >= 0) rc_integrate(s[n2], s[n2], rcint_tau);
      if (seg[n2] < 0 || seg[n1] < 0){
	n1 = n2;
	n2 = !n1;
	continue;
      }
      f2 = fit_time(nsegments, TSTEPS_TOT, s2[n1], s2[n2], 0);
      if (f2 > C2_THRESH){
	f2 = 2.354*RZ_STEP*sqrt(C2_SIGNAL_MULT_F)/sqrt(f2);
	fr2[i-1][j] = f2;
      }else{
	fr2[i-1][j] = -1;
      }
      f = 0;
      for (k = 0; k < nsegments; k++){
	for (l = 0; l < time_steps; l++){
	  d = C2_SIGNAL_MULT_F*(s[n1][k][l] - s[n2][k][l]);
	  f += d*d;
	}
      }
      if (f > C2_THRESH){
	f = 2.354*RZ_STEP*sqrt(C2_SIGNAL_MULT_F)/sqrt(f);	
	fr[i-1][j] = f;
      }else{
	fr[i-1][j] = -1;
      }
      n1 = n2;
      n2 = !n1;
    }
  }
  printf("\n");
  for (i = 0; i < nr; i++){
    printf("\r%d/%d", i+1, nr);fflush(stdout);
    cyl[n1].r = cyl[n2].r = 0 + i*RZ_STEP;
    cyl[n1].z = 0;
    cart[n1] = cyl_to_cart(cyl[n1]);
    seg[n1] = get_signal(cart[n1], s[n1]);
    if (rcint_tau > 0 && seg[n1] >= 0) rc_integrate(s[n1], s[n1], rcint_tau);
    for (j = 1; j < nz; j++){
      cyl[n2].z = 0 + j*RZ_STEP;
      cart[n2] = cyl_to_cart(cyl[n2]);
      seg[n2] = get_signal(cart[n2], s[n2]);
      if (rcint_tau > 0 && seg[n2] >= 0) rc_integrate(s[n2], s[n2], rcint_tau);
      if (seg[n2] < 0 || seg[n1] < 0){
	n1 = n2;
	n2 = !n1;
	continue;
      }
      f2 = fit_time(nsegments, TSTEPS_TOT, s2[n1], s2[n2],0);
      if (f2 > C2_THRESH){
	f2 = 2.354*RZ_STEP*sqrt(C2_SIGNAL_MULT_F)/sqrt(f2);
	fz2[i][j-1] = f2;
      }else{
	fz2[i][j-1] = -1;
      }      f = 0;
      for (k = 0; k < nsegments; k++){
	for (l = 0; l < time_steps; l++){
	  d = C2_SIGNAL_MULT_F*(s[n1][k][l] - s[n2][k][l]);
	  f += d*d;
	}
      }
      if (f > C2_THRESH){
	f = 2.354*RZ_STEP*sqrt(C2_SIGNAL_MULT_F)/sqrt(f);	
	fz[i][j-1] = f;
      }else{
	fz[i][j-1] = -1;
      }
      n1 = n2;
      n2 = !n1;
    }
  }
  printf("\n");
  fprintf(fp,"r z fr fz  fr2 fz2--- phi = %f\n", pval);
  for (i = 0; i < nr; i++){
    cyl[n1].r = RZ_STEP/2 + i*RZ_STEP;
    for (j = 0; j < nz; j++){
      cyl[n1].z = RZ_STEP/2 + j*RZ_STEP;
      fprintf(fp, "%.2f %.2f %f %f %f %f\n", cyl[n1].r, cyl[n1].z, fr[i][j], 
	      fz[i][j], fr2[i][j], fz2[i][j]);
    }
    fprintf(fp, "\n");
  }
  for (i = 0; i < 2; i++){
    for (j = 0; j < nsegments; j++){
      free(s2[i][j]);    
    }
    free(s[i]);
    free(s2[i]);
  }
  free(s);
  free(s2);

  for (i = 0; i < nr; i++){
    free(fr[i]);
    free(fz[i]);
    free(fr2[i]);
    free(fz2[i]);
  }
  free(fr);
  free(fz);
  free(fr2);
  free(fz2);
  return 0;
}

static int c2d_pz(float rval, FILE *fp, float rcint_tau){
  float ***s, ***s2;
  float **fphi, **fz, f, **fphi2, **fz2, f2;
  int np, nz;
  float pmax, zmax,  d;
  int i,j, k, l;
  cyl_pt cyl[2];
  point cart[2];
  int seg[2];
  int n1, n2;

  pmax = 2*M_PI - P_STEP/2;
  zmax = zmax_detector();
  np = pmax/P_STEP + 2;
  nz = zmax/RZ_STEP + 1;
  
  if ((s = malloc(2*sizeof(*s))) == NULL
      ||(s2 = malloc(2*sizeof(*s2))) == NULL){
    printf("malloc failed\n");
    exit(1);
  }
  for (i = 0; i < 2; i++){
    if ((s[i] = malloc(nsegments*sizeof(*s[i]))) == NULL
	||(s2[i] = malloc(nsegments*sizeof(*s2[i]))) == NULL){
      printf("malloc failed\n");
      exit(1);
    }    
    for (j = 0; j < nsegments; j++){
      if ((s2[i][j] = malloc(TSTEPS_TOT*sizeof(*s[i][j]))) == NULL){
	printf("malloc failed\n");
	exit(1);
      }    
      s[i][j] = s2[i][j] + (TSTEPS_TOT - time_steps);
      memset(s2[i][j], 0, sizeof(*s2[i][j])*TSTEPS_TOT);
    }
  }
  if ((fphi = malloc(np*sizeof(*fphi))) == NULL
      || (fz = malloc(np*sizeof(*fz))) == NULL
      || (fphi2 = malloc(np*sizeof(*fphi2))) == NULL
      || (fz2 = malloc(np*sizeof(*fz2))) == NULL){
    printf("malloc failed\n");
    exit(1);
  }
  for (i = 0; i < np; i++){
    if ((fphi[i] = malloc(nz*sizeof(*fphi[i]))) == NULL
	|| (fz[i] = malloc(nz*sizeof(*fz[i]))) == NULL
	|| (fphi2[i] = malloc(nz*sizeof(*fphi2[i]))) == NULL
	|| (fz2[i] = malloc(nz*sizeof(*fz2[i]))) == NULL){
      printf("malloc failed\n");
      exit(1);
    }
    for (j = 0; j < nz; j++){
      fphi[i][j] = fz[i][j] = fphi2[i][j] = fz2[i][j] = -2;
    }
  }
  n1 = 0; n2 = 1;
  cyl[n1].r = cyl[n2].r = rval;
  for (j = 0; j < nz; j++){
    printf("\r%d/%d", j+1, nz);fflush(stdout);
    cyl[n1].phi = 0;
    cyl[n1].z = cyl[n2].z = 0 + j*RZ_STEP;
    cart[n1] = cyl_to_cart(cyl[n1]);
    seg[n1] = get_signal(cart[n1], s[n1]);
    if (rcint_tau > 0 && seg[n1] >= 0) rc_integrate(s[n1], s[n1], rcint_tau);
    for (i = 1; i < np; i++){
      cyl[n2].phi = 0 + i*P_STEP;
      cart[n2] = cyl_to_cart(cyl[n2]);
      seg[n2] = get_signal(cart[n2], s[n2]);
      if (rcint_tau > 0 && seg[n2] >= 0) rc_integrate(s[n2], s[n2], rcint_tau);
      if (seg[n2] < 0 || seg[n1] < 0){
	n1 = n2;
	n2 = !n1;
	continue;
      }
      f2 = fit_time(nsegments, TSTEPS_TOT, s2[n1], s2[n2],0);
      if (f2 > C2_THRESH){
	f2 = 2.354*P_STEP*cyl[n1].r*sqrt(C2_SIGNAL_MULT_F)/sqrt(f2);
	fphi2[i-1][j] = f2;
      }else{
	fphi2[i-1][j] = -1;
      }      f = 0;
      for (k = 0; k < nsegments; k++){
	for (l = 0; l < time_steps; l++){
	  d = C2_SIGNAL_MULT_F*(s[n1][k][l] - s[n2][k][l]);
	  f += d*d;
	}
      }
      if (f > C2_THRESH){
	f = 2.354*P_STEP*cyl[n1].r*sqrt(C2_SIGNAL_MULT_F)/sqrt(f);	
	fphi[i-1][j] = f;
      }else{
	fphi[i-1][j] = -1;
      }
      n1 = n2;
      n2 = !n1;
    }
  }
  printf("\n");
  for (i = 0; i < np; i++){
    printf("\r%d/%d", i+1, np);fflush(stdout);
    cyl[n1].phi = cyl[n2].phi = 0 + i*P_STEP;
    cyl[n1].z = 0;
    cart[n1] = cyl_to_cart(cyl[n1]);
    seg[n1] = get_signal(cart[n1], s[n1]);
    if (rcint_tau > 0 && seg[n1] >= 0) rc_integrate(s[n1], s[n1], rcint_tau);
    for (j = 1; j < nz; j++){
      cyl[n2].z = 0 + j*RZ_STEP;
      cart[n2] = cyl_to_cart(cyl[n2]);
      seg[n2] = get_signal(cart[n2], s[n2]);
      if (rcint_tau > 0 && seg[n2] >= 0) rc_integrate(s[n2], s[n2], rcint_tau);
      if (seg[n2] < 0 || seg[n1] < 0){
	n1 = n2;
	n2 = !n1;
	continue;
      }
      f2 = fit_time(nsegments, TSTEPS_TOT, s2[n1], s2[n2], 0);
      if (f2 > C2_THRESH){
	f2 = 2.354*RZ_STEP*sqrt(C2_SIGNAL_MULT_F)/sqrt(f2);
	fz2[i][j-1] = f2;
      }else{
	fz2[i][j-1] = -1;
      }   
      f = 0;
      for (k = 0; k < nsegments; k++){
	for (l = 0; l < time_steps; l++){
	  d = C2_SIGNAL_MULT_F*(s[n1][k][l] - s[n2][k][l]);
	  f += d*d;
	}
      }
      if (f > C2_THRESH){
	f = 2.354*RZ_STEP*sqrt(C2_SIGNAL_MULT_F)/sqrt(f);	
	fz[i][j-1] = f;
      }else{
	fz[i][j-1] = -1;
      }
      n1 = n2;
      n2 = !n1;
    }
  }
  printf("\n");
  fprintf(fp,"p z fp fz fp2 fz2 --- phi = %f\n", rval);
  for (i = 0; i < np; i++){
    cyl[n1].phi = P_STEP/2 + i*P_STEP;
    for (j = 0; j < nz; j++){
      cyl[n1].z = RZ_STEP/2 + j*RZ_STEP;
      fprintf(fp, "%.4f %.2f %f %f %f %f\n", cyl[n1].phi, cyl[n1].z, 
	      fphi[i][j], fz[i][j], fphi2[i][j], fz2[i][j]);
    }
    fprintf(fp, "\n");
  }
  for (i = 0; i < 2; i++){
    for (j = 0; j < nsegments; j++){
      free(s2[i][j]);    
    }
    free(s[i]);
    free(s2[i]);
  }
  free(s);
  free(s2);

  for (i = 0; i < np; i++){
    free(fphi[i]);
    free(fz[i]);
    free(fphi2[i]);
    free(fz2[i]);
  }
  free(fphi);
  free(fz);
  free(fphi2);
  free(fz2);
  return 0;

}

/* calculates min chisq for fitted t-zero,
   trying to match two signals s1 and s2.
   return min chisq value.
   from David's posres.c
*/
static float fit_time(int nsigs, int ntimes, float **s1, float **s2, int shift_s1)
{
  int   i, t, dt=0, lo, hi;
  double d, chisq, dsdt, fdt, chisq2, prev_chisq2, sum1, sum2, fdt2=0, prev_fdt2;
  float retval;


  chisq2 = prev_chisq2 = 100*nsigs*ntimes;

  while (1) {
    chisq = sum1 = sum2 = 0.0f;
    if (dt > 0) {
      lo = 0; hi = ntimes-dt-1;
    } else {
      lo = -dt; hi = ntimes-1;
    }
    for (i=0; i<nsigs; i++) {
      for (t=lo; t<hi; t++) {
	d = C2_SIGNAL_MULT_F*(s1[i][t+dt] - s2[i][t]);
	dsdt = C2_SIGNAL_MULT_F*(s1[i][t+dt+1] - s1[i][t+dt]);
	chisq += d*d;
	sum1 += d*dsdt;
	sum2 += dsdt*dsdt;
      }
    }
    if (chisq == 0.0f || sum2 == 0.0f) return chisq;   // deriv = 0 => constant signal

    prev_chisq2 = chisq2;
    prev_fdt2 = fdt2;
    fdt = -sum1/sum2;
    if (fdt >= 0.0f) {
      if (fdt <= 1.0f){ //return (chisq + (2.0*sum1 + sum2*fdt)*fdt);
	retval = (chisq + (2.0*sum1 + sum2*fdt)*fdt);
	fdt += dt;
	break;
      }
      if (dt > ntimes/4){ //return chisq;      // best-fit dt-zero is too big
	retval = chisq;
	fdt += dt;
	break;
      }
      chisq2 = chisq + 2.0*sum1 + sum2;
      fdt2 = dt + 0.9999;
    } else {
      if (dt < -1 - ntimes/4){// return chisq; // best-fit dt-zero is too big
	retval = chisq;
	fdt += dt;
	break;
      }
      chisq2 = chisq;
      fdt2 = dt;
    }
    fdt += (float) dt;
    dt = (int) fdt;
    if (fdt < 0.0f) dt--;

    if (dt > ntimes/4)      dt = ntimes/4;
    if (dt < -1 - ntimes/4) dt = -1 - ntimes/4;
    if (prev_chisq2 <= chisq2){
      retval = prev_chisq2;
      fdt = prev_fdt2;
      break;
    }
  } 

  if (shift_s1) shift_signal(s1, fdt, nsigs, ntimes);

  return retval;
}


int shift_signal(float **s, float dt, int nsigs, int ntimes){
  float *stmp, fdt;
  int idt;
  int i, j;
  
  idt = trunc(dt); // to nearest int not larger in absolute value
  fdt = fabs(dt - idt);

  if ((stmp = malloc(ntimes *sizeof(*stmp))) == NULL){
    printf("malloc failed\n");
    exit(1);
  }
  for (i = 0; i < nsigs; i++){
    for (j = 0; j < ntimes; j++){
      stmp[j] = s[i][j];
      s[i][j] = 0.0;
    }
    for (j = 0; j < ntimes; j++){
      if (j + idt < 0) continue;
      if (j + idt >= ntimes) continue;
      s[i][j] += stmp[j+idt]*(1-fdt);
      if (j +1 < ntimes) s[i][j+1] += stmp[j+idt]*fdt;
    }
  }
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

static int fit_t0(char *cmd){
  struct point cart1, cart2;
  struct cyl_pt cyl1, cyl2;
  int seg1, seg2;
  static float **s1, **s2;
  int i,j;
  float c2d;
  char fn[MAX_LINE], *cp, *cp2;
  float spec[MAX_SPE_CHS];
  int comp_f;

  if (s1 == NULL){
    if ((s1 = malloc(nsegments*sizeof(*s1))) == NULL
      || (s2 = malloc(nsegments*sizeof(*s2))) == NULL){
      printf("malloc failed\n");
    }
    for (i = 0; i < nsegments; i++){
      if ((s1[i] = malloc(time_steps*sizeof(*s1[i]))) == NULL
	|| (s2[i] = malloc(time_steps*sizeof(*s2[i]))) == NULL){
	printf("malloc failed\n");
      }

    }
  }

  if (coord_type == CYL){
    if (sscanf(cmd, "%f %f %f %f %f %f", 
	       &cyl1.r, &cyl1.phi, &cyl1.z,
	       &cyl2.r, &cyl2.phi, &cyl2.z) == 6
	||sscanf(cmd, "%f, %f, %f, %f, %f, %f", 
		 &cyl1.r, &cyl1.phi, &cyl1.z,
		 &cyl2.r, &cyl2.phi, &cyl2.z) == 6){
      
      cart1 = cyl_to_cart(cyl1);
      cart2 = cyl_to_cart(cyl2);
    }else{
      fprintf(stderr,"Failed to parse coordinates: %s\n", cmd);
      return 1;
    }
  }else{
    if (sscanf(cmd, "%f %f %f %f %f %f", 
	       &cart1.x, &cart1.y, &cart1.z,
	       &cart2.x, &cart2.y, &cart2.z) == 6
	||sscanf(cmd, "%f, %f, %f, %f, %f, %f", 
		 &cart1.x, &cart1.y, &cart1.z,
		 &cart2.x, &cart2.y, &cart2.z) == 6){
      ;//nothing
    }else{
      fprintf(stderr,"Failed to parse coordinates: %s\n", cmd);
      return 1;
    }
  }
  for (cp = cmd; isspace(*cp); cp++);
  for ( ; isdigit(*cp) || *cp == '.'; cp++);//skip r/x coord.
  for ( ; isspace(*cp) || *cp == ','; cp++); //skip whitespace, commas
  for ( ; isdigit(*cp) || *cp == '.'; cp++);//skip phi/y coord.
  for ( ; isspace(*cp) || *cp == ','; cp++);//skip whitespace, commas
  for ( ; isdigit(*cp) || *cp == '.'; cp++);//skip z coord
  for ( ; isspace(*cp)              ; cp++);//skip whitespace
  for ( ; isdigit(*cp) || *cp == '.'; cp++);//skip r/x coord.
  for ( ; isspace(*cp) || *cp == ','; cp++); //skip whitespace, commas
  for ( ; isdigit(*cp) || *cp == '.'; cp++);//skip phi/y coord.
  for ( ; isspace(*cp) || *cp == ','; cp++);//skip whitespace, commas
  for ( ; isdigit(*cp) || *cp == '.'; cp++);//skip z coord

  for ( ; isspace(*cp); cp++); //skip whitespace
  for (cp2 = cp + strlen(cp); isspace(*cp2); cp2--)//remove whitespace
    *cp2 = '\0';
  if (strlen(cp) == 0){
    fprintf(stderr, "must supply file name\n");
    return 1;
  }
  strcpy(fn, cp);

  printf("finding signal for (x = %.1f, y = %.1f, z = %.1f)\n",
	 cart1.x, cart1.y, cart1.z);      
  seg1 = get_signal(cart1, s1);
  if (seg1 < 0) {
    printf("not a valid coordinate\n");
    return 1;
  }
  printf("signal 1 is for segment no %d\n", seg1);

  printf("finding signal for (x = %.1f, y = %.1f, z = %.1f)\n",
	 cart2.x, cart2.y, cart2.z);      
  seg2 = get_signal(cart2, s2);
  if (seg2 < 0) {
    printf("not a valid coordinate\n");
    return 1;
  }
  printf("signal 2 is for segment no %d\n", seg2);

  c2d = fit_time(nsegments, time_steps, s1, s2, 1); 

  for (comp_f = 1; nsegments*time_steps/comp_f > MAX_SPE_CHS; comp_f *= 2)
    ;

  for (i = 0; i < sizeof(spec)/sizeof(spec[0]); i++)
    spec[i] = 0.0;
  /*copy signal data to spectrum array*/
  for (i = 0; i < nsegments; i++){
    for (j = 0; j < time_steps; j++){
      spec[(i*time_steps+j)/comp_f] += s1[i][j]*1000/comp_f;
    }
  }
  printf("Saving spectrum as %s, compression f. : %d\n", fn, comp_f);
  write_spectrum(spec, nsegments*time_steps/comp_f, fn);

  return 0;
}
