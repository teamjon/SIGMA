/* scan_ppc.c (Ren Cooper, Nov 09)

   Program to generate signal database from a PPC detector
   
   Uses signal calculation codes written by Karin Lagergren

   ** updated to include ability to sample solid angle

  to compile: 
   gcc -o st scan_ppc.c point.c cyl_point.c calc_signal.c\
     fields.c detector_geometry.c signal_calc_util.c -lm -lreadline
 

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double trunc(double x);
double atan2(double a, double b);
#include <readline/readline.h>
#include <readline/history.h>
#include <ctype.h>

#include "calc_signal.h"
#include "cyl_point.h"
#include "signal_calc_util.h"
#include "detector_geometry.h"

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

#define MAX_SPE_CHS 16384

#define FREE(p) do{free(p); p = NULL; }while(0)

#define CYL 0
#define CART 1

#define C2_SIGNAL_MULT_F 10000
#define C2_THRESH 100

#define NON_LINEAR_CHISQ 1

//static int coord_type = CYL;
static int time_steps;
static int nsegments;

static float rmin = 0.0;
static float rmax = 33.5;
static float rstep = 0.5;
static float zmin = 54.5; //20.0;
static float zmax = 57.5; //22.0;
static float zstep = 0.1;
static float phimin = 0.0;
static float phimax = 0.785398;   //  -> 45 degrees 2*3.14;
static float phistep = 0.0174533; // -> 1 degree steps
static float t = 31.51;


static int rc_integrate(float **s_in, float **s_out, float tau);
static int get_rt(int *sig);

int main(){

  struct point cart;
  struct cyl_pt cyl;
  int sgnl[4096];
  int seg;
  int i;
  int k;
  int posns = 0;
  int slow;
  int slow_sigs = 0;
  int fast_sigs = 0;
  int longrt = 0;
  int shortrt = 0;
  float r,z,phi;
  int rt;
  
  int X = 0.0;
  double sumangles = 0.0;
  double angles[35];
  int events[35];

  int ii;
  float total_solid_angle;
  float solid_angle[70];
  float sum_solid_angles = 0.0;
  int num_events = 0;
  int loop_counter;

  //char *cp, *cp2;
  static float **s;
  FILE *file_out;

  printf("\n");
  printf("===== scan_ppc.c ====\n");
  printf("\n");
  printf(" - code to generate database of calculated signals\n");
  printf(" setup info....\n");
  printf("\n");

#ifdef CHATTY
  set_signal_calc_output(ANNOYINGLY_CHATTY, vprintf);
#endif
  if (signal_calc_init(GEOMETRY_FILE, FIELDS_DIR_FILE,
		       SIGNAL_PARS_FILE, &time_steps, &nsegments) != 0){
    return 1;
  }


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



  /* --------- calculate solid angle distribution --------- */
  float Dist = 50.0; //50.0
  printf("\n");
  printf("- calculating solid angle\n");  

  total_solid_angle = 2.0 * 3.14 * (1.0 - (Dist/(sqrt(( Dist * Dist) + (33.5 * 33.5)))));
  printf("total solid angle = %f\n",total_solid_angle);

  float sum1 = 0.0;
  float sum2 = 0.0;
  /* calculate solid angle contained within each radius */
  i = 0;
  for (r=0.0; r<rmax+rstep; r=r+rstep){
    solid_angle[i] = 2.0 * 3.14 * (1.0 - (Dist/(sqrt(( Dist * Dist) + (r * r))))) - sum_solid_angles;
    sum_solid_angles = sum_solid_angles + solid_angle[i];
    
    if (r<5.0) sum1 = sum1 + solid_angle[i];
    if (r>=5.0) sum2 = sum2 + solid_angle[i];

    i++;
  }
  
  //for (i=0; i<50; i++) printf("%d %f\n",i, solid_angle[i]);
  //printf("%f\n",sum2/sum1);
  
/*
  i = 0;
  sumangles = 0;
  for (r=rmin; r<rmax+(2.0*rstep); r=r+(2.0*rstep)){
    if (i == 0) X = atan2 ((double)r,10.0) * 180.0/3.14;
    if (i != 0) X = ((atan2((double)r,10.0)) * 180.0/3.14) - sumangles;

    angles[i] = X;
    sumangles = sumangles + X;
    i++;
    printf("%.2f %f\n",r,X);
  }
  printf("%f\n",sumangles);
  
  for (i=0; i<35; i++){
    events[i] = (int) (1000.0f * angles[i]/sumangles);
    printf("%d %f %d\n",i,angles[i],events[i]);
  }
  */
  /* ------------------------------------------------------- */


  /* -------------- do signal generation/scan ---------------- */

  printf("\n");
  printf("- calculating signals with %d steps\n",time_steps);
  printf("\n");
  
  /* open file */
  file_out = fopen("sigs.spn","w+");

  printf(" - scanning r = %.2f to %.2f, z = %.2f to %.2f, phi = %.2f to %.2f\n",rmin,rmax,zmin,zmax,phimin,phimax);
  printf("\n");

  /* loop over z and r */
  loop_counter = 0;
  for (r=rmin+(0.5*rstep); r<rmax+rstep; r=r+rstep){
    
    /* events at each radius - ignore core region */
    if (r > 5.0) {
      num_events = (int)(((solid_angle[loop_counter]/total_solid_angle) * 5000.0f)/45.0f);
    }
    if (r <= 5.0) num_events = 0;

    //printf("radius = %.2f, %d events at each phi\n",r,num_events);
    loop_counter++;

    for (z=zmin; z<zmax+zstep; z=z+zstep){
      for (phi=phimin; phi<phimax+zstep; phi=phi+phistep){
      
	cyl.r = r;
	cyl.phi = phi;
	cyl.z = z;

       
      /* zero signal array */
      for (i=0; i<4096; i++) sgnl[i] = 0;

      /* convert to cartesian coords*/
      cart = cyl_to_cart(cyl);

      /* get signal */
      seg = get_signal(cart, s);
      if (seg < 0) {
	//printf("point not in crystal: (x = %.1f, y = %.1f, z = %.1f)\n", cart.x, cart.y, cart.z);
	//printf("%.1f  %.1f %.2f\n",(float)z,(float)r, 10.0f);
	//printf("%.1f  %.1f  %.1f\n",r,z,10.0f);
	//return 1;
	}

      if (seg >= 0) {
	/* preamp */
	rc_integrate(&s[seg], &s[seg], 46.0); //36.0

	/* copy to array */
	k = 0;
	for (i=0; i<time_steps; i++) sgnl[i] = s[seg][i] * 1000;
	//for (i=0; i<time_steps; i++) printf("%d\n",sgnl[i]);
	
	//sgnl[3000] = r;
	//sgnl[3001] = phi;
	//sgnl[3002] = z;

	  /* get risetime */
	  rt = get_rt(sgnl);
	  
	  /* print sig # and position */
	  //printf("%d %d %d\n",posns,r,z);

	  /* write to output file */
	for (ii=0; ii<num_events; ii++) {
	  X++;
	  //printf("%d %f %f\n",X,r,z);
	  fwrite(sgnl,sizeof(sgnl),1,file_out);
	}

	  /* increment counter */
	  posns++;

      }


      }
    }
    //printf("\n");
  }

  fclose(file_out);

  printf("\n");
  printf("%d signals written to file\n",X);
  printf("\n");
 
 return 0;
}


/* FUNCTIONS */

static int rc_integrate(float **s_in, float **s_out, float tau){
  int i, j;
  float s_in_old, s;  /* DCR: added so that it's okay to
			 call this function with s_out == s_in */
  
  for (i = 0; i < nsegments; i++){
    s_in_old = s_in[i][0];
    s_out[i][0] = 0.0;
    for (j = 1; j < time_steps; j++){
      s = s_out[i][j-1] + (s_in_old - s_out[i][j-1])/tau;//Apr 26, KL -- BUG
      s_in_old = s_in[i][j];
      s_out[i][j] = s;
    }
  }
  return 0;
}

static int get_rt(int *sig){
  int t50 = 0;
  int t95 = 0;

  while (sig[t50] < 500.0) t50++;

  return t50;
}
