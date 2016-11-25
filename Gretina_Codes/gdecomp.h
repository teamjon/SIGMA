/* Header file for GRETA decomposition routines
   Author:  D.C. Radford    Aug 2004
 */

#ifndef _gdecomp_h
#define _gdecomp_h

#define BASIS_FILE  "basis.dat"         /* file containing basis data */
#define EVENT_FILE  "input_signals.dat" /* CHANGEME file containing observed data */
#define SUMS_FILE   "precalc.dat"       /* file containing precalculated sums */

#define GRID_PTS   226993  /* CHANGEME number of grid points in the basis */
#define GRID_SEGS  37      /* number of signals calculated for each basis point */
#define MEAS_SEGS  37      /* number of signals measured for each event */
#define TIME_STEPS 50      /* number of time steps calculated/measured for each segment */
#define TOT_SEGS   36      /* total number of segments in crystal */
#define MAX_AGS    550     /* CHANGEME max. no. of points in coarse grid for AGS */

#define MAX_SEGS   8              /* max. number of segments to take in events */
#define MAX_PARS   (8*MAX_SEGS+1) /* max. number of parameters to fit */
#define MAX_INTPTS (2*MAX_SEGS)   /* max. number of interaction points */
#define SSEG       TOT_SEGS       /* range of seg in basis grid */
#define SRAD       36             /* range of r in basis grid */
#define SPHI       13             /* range of phi in basis grid */
#define SZZZ       20             /* range of z in basis grid */

#define COAL_DIST_DEFAULT 2.0


/* data structure for measured event signals */
typedef struct {
  float total_energy;                    /* total energy in crystal */
  float seg_energy[TOT_SEGS];            /* net energy in each segment */
  float signal[MEAS_SEGS][TIME_STEPS];   /* actual measured signals */
#ifdef TIMESTAMP
  long long int time;                    /* timestamp for this crystal event */
#endif
} Event_Signal;

/* data structure for calculated basis signals */
typedef struct {
  char  iseg, ir, ip, iz;             /* integer cylindrical coordinates of grid point */
  float x, y, z;                      /* cartesian coordinates of grid point */
  float signal[GRID_SEGS][50];        /* actual basis signals */
  int   lo_time[GRID_SEGS], hi_time[GRID_SEGS];     /* limits for non-zero signal */
} Basis_Point;

/* data structure for Adaptive Grid Search coarse grid */
typedef struct {
  int    npts;                /* number of AGS points for this segment */
  int    grid_pos[MAX_AGS];   /* pointer to basis-signal-ID for each AGS point */
  double *da;                 /* precalculated sums */
  float  *db;                 /* precalculated sums */
} Adaptive_Grid;

/* data structure for interactions */
typedef struct {
  int    seg;                 /* segment */
  int    pos;                 /* basis signal position, if applicable, else -1 */
  double r, p, z, e;          /* parameters */
  double dr, dp, dz, de;      /* uncertainties */
} Interaction;

extern Basis_Point   *basis;                 /* basis-signal data */
extern int           grid_pos_lu[SSEG][SRAD][SPHI][SZZZ];    /* basis-grid position look-up table */
extern int           maxir[SSEG], maxip[SSEG], maxiz[SSEG];  /* max. values of ir, ip, iz for each segment */
extern Adaptive_Grid ags1[SSEG];             /* Adaptive Grid Search coarse grid, 1 segment */

extern int           quiet;                  /* set to 1 to prevent output of diagnostic info */
extern int           average_sigs;           /* set to 1 to output averaged observed and 
						fitted signals for single-segment (net=1) events */
extern int           *bad_segs;              /* list of segment numbers that should be disabled */

/* ===================
   function prototypes
   =================== */

/* many of these modules could be combined... */

/* in read_basis.c: */
int    read_basis(char *basis_file_name);

/* in gdecomp.c: */
void   pclock(int print_flag);

/* in decompose.c: */
/* returns number of best-fit interactions found */
int decompose_1(Event_Signal asig,  /* observed signals */
		Event_Signal *bsig, /* fitted signals */
		int seg, Interaction *ints, double *t0, double *chisq_out, /* parameters */
		int grid2, int fit0, int fit1, int fit2, int fit3,         /* control */
		int fit4, int final_fit, int coalesce, double min_e_fraction);  /* control */
/* returns number of best-fit interactions found */
int decompose_n(Event_Signal asig,  /* observed signals*/
		Event_Signal *bsig, /* fitted signals*/
		int nseg, int *seg, int coalesce, 
		Interaction *ints, double *t0, double *chisq_out);  /* parameters */
/* returns final number of interactions found */
int postprocess_events(Interaction *ints, int nints, float total_e,
		       int ppflag, float coal_dist,
		       double *x, double *y, double *z, double *e);

/* in grid.c: */
int    grid_init(void);
double coarse_grid_1(Event_Signal asig,  /* observed signals */
		     int seg, Interaction *ints,
		     double *chisq0, double min_e_fraction);
double refine_grid_1(Event_Signal asig,  /* observed signals */
		     double chisq, double chisq0, double min_e_fraction,
		     Interaction *ints);

/* in fitter.c: */
double fitter(Event_Signal asig,  /* observed signals */
	      Event_Signal *bsig, /* fitted signals */
	      Interaction *ints, double *t0, int nints, int flag);

/* in eval.c: */
int    eval(Event_Signal asig,  /* observed signals */
	    Event_Signal *bsig, /* fitted signals */
	    Interaction *ints, double t0, int nints,
	    double *chisq_out, double beta[MAX_PARS],
	    double alpha[MAX_PARS][MAX_PARS], int calc_deriv, int *ssel);
double eval_int_pos(Event_Signal asig,  /* observed signals */
		    Event_Signal *bsig, /* fitted signals */
		    Interaction *ints, double t0, int nints);

/* in interpolate.c: */
int    nearest_grid_points(int seg, double rin, double pin, double zin,
			   float *rdiff, float *pdiff, float *zdiff, int pos_out[8]);
int    interpolate(int seg, double rin, double pin, double zin, Basis_Point *sig,
		   Basis_Point sigderiv[3], int calc_deriv, int *ssel);
int    cyl_to_cart(int seg, double *pars, double*x, double *y, double*z, double *e);

/* in matinv.c: */
int    matinv(double *array, int norder, int dim);

/* data structs added to allow multithreading */

struct decomp_errcnt {
  int nonet;
  int toomanynet;
  int sumener;
  int badchisq;
};

struct decomp_state {
  Event_Signal *bsig;
  Interaction  *ints;
  int   cnt;
  float coal_dist;	
  struct decomp_errcnt *err;
};

struct crys_intpts {
  int   num;
  float tot_e;
  float t0;
  float chisq;
  float norm_chisq;
  long long int timestamp;
  struct {
    float x, y, z;
    float e;
  } intpts[MAX_INTPTS];
};

/* in decompose.c: */
/* routine to do overall initialization of a decomposition thread;
   returns a struct handle into the state of the new decomposition thread*/
struct decomp_state *dl_decomp_init(char *basis_file_name, int set_quiet);
struct crys_intpts *dl_decomp(struct decomp_state *di, Event_Signal *asig);
char *dl_crys_intpts_2s(struct crys_intpts *x);
void dl_set_coal_dist(struct decomp_state *inst, float d);
int num_net(Event_Signal *asig);


#endif	/* _gdecomp_h */
