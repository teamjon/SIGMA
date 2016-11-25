/* Header file for GRETA decomposition routines
   Author:  D.C. Radford    Aug 2004

   Modified for prototype segmented inverted coax
   D.C. Radford  Feb 2015
 */

#ifndef _pdecomp_h
#define _pdecomp_h

#define BASIS_FILE  "basis.dat"         /* file containing basis data */
#define EVENT_FILE  "input_signals.dat" /* CHANGEME file containing observed data */
#define SUMS_FILE   "precalc.dat"       /* file containing precalculated sums */

#define GRID_PTS    33092// 35303   /* CHANGEME number of grid points in the basis */
#define NUM_SIGS    20       /* number of signals calculated for each basis point
				and measured for each event */
#define TIME_STEPS_C 200     /* number of time steps calculated for each segment (basis) */
#define TIME_STEPS_M 300     /* number of time steps measured for each segment (event) */
#define TOT_SEGS     19      /* total number of segments in crystal */
#define MAX_AGS      550     /* CHANGEME max. no. of points in coarse grid for AGS */

#define MAX_SEGS   8              /* max. number of segments to take in events */
#define MAX_PARS   (8*MAX_SEGS+1) /* max. number of parameters to fit */
#define MAX_INTPTS (2*MAX_SEGS)   /* max. number of interaction points */
#define SRAD       36             /* range of r in basis grid */
#define SPHI       16             /* range of phi in basis grid */
#define SZZZ       81             /* range of z in basis grid */

#define COAL_DIST_DEFAULT 2.0


/* data structure for measured event signals */
typedef struct {
  float total_energy;                    /* total energy in crystal */
  float seg_energy[TOT_SEGS];            /* net energy in each segment */
  float signal[NUM_SIGS][TIME_STEPS_M];  /* actual measured signals */
#ifdef TIMESTAMP
  long long int time;                    /* timestamp for this crystal event */
#endif
} Event_Signal;

/* data structure for calculated basis signals */
typedef struct {
  short iseg, ir, ip, iz;                /* integer cylindrical coordinates of grid point */
  float x, y, z;                         /* cartesian coordinates of grid point */
  float signal[NUM_SIGS][TIME_STEPS_C];  /* actual basis signals */
  short lo_time[NUM_SIGS], hi_time[NUM_SIGS];   /* limits for non-zero signal */
  short t80, t20;                        /* PC and segment signal times, in ns */
} Basis_Point;

/* data structure for Adaptive Grid Search coarse grid */
typedef struct {
  int     npts;                /* number of AGS points for this segment */
  int     grid_pos[MAX_AGS];   /* pointer to basis-signal-ID for each AGS point */
  double *da;                  /* precalculated sums */
  float  *db;                  /* precalculated sums */
} Adaptive_Grid;

/* data structure for interactions */
typedef struct {
  int    seg;                 /* segment */
  int    pos;                 /* basis signal position, if applicable, else -1 */
  double r, p, z, e;          /* parameters */
  double dr, dp, dz, de;      /* uncertainties */
} Interaction;

typedef struct {
  Basis_Point   *basis;                /* basis-signal data */
  int           grid_pos_lu[SRAD][SPHI][SZZZ];      /* basis-grid position look-up table */
  int           seg_basis_num[GRID_PTS][NUM_SIGS];  /* list of grid IDs inside each segment */
  int           t_order[GRID_PTS];                  /* list of grid IDs ordered by drift time */
  int           quiet;                 /* set to 1 to prevent output of diagnostic info */
  int           average_sigs;          /* set to 1 to output averaged observed and 
 				          fitted signals for single-segment (net=1) events */
  int           *bad_segs;             /* list of segment numbers that should be disabled */
  // Adaptive_Grid ags1[NUM_SIGS];     /* Adaptive Grid Search coarse grid, 1 segment */
} PDecomp;

/* ===================
   function prototypes
   =================== */

/* many of these modules could be combined... */

/* in read_basis.c: */
int    read_basis(char *basis_file_name, PDecomp *pdd);

/* in pdecomp.c: */
void   pclock(int print_flag);

/* in decompose.c: */
/* returns number of best-fit interactions found */
int decompose_1(PDecomp *pdd,
		Event_Signal asig,  /* observed signals */
		Event_Signal *bsig, /* fitted signals */
		int seg, Interaction *ints, double *t0, double *chisq_out, /* parameters */
		int grid2, int fit0, int fit1, int fit2, int fit3,         /* control */
		int fit4, int final_fit, int coalesce, double min_e_fraction);  /* control */
/* returns number of best-fit interactions found */
int decompose_n(PDecomp *pdd,
		Event_Signal asig,  /* observed signals*/
		Event_Signal *bsig, /* fitted signals*/
		int nseg, int *seg, int coalesce, 
		Interaction *ints, double *t0, double *chisq_out);  /* parameters */
/* returns final number of interactions found */
int postprocess_events(PDecomp *pdd,
		       Interaction *ints, int nints, float total_e,
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
double fitter(PDecomp *pdd,
	      Event_Signal asig,  /* observed signals */
	      Event_Signal *bsig, /* fitted signals */
	      Interaction *ints, double *t0, int nints, int flag);

/* in eval.c: */
int    eval(PDecomp *pdd,
	    Event_Signal asig,  /* observed signals */
	    Event_Signal *bsig, /* fitted signals */
	    Interaction *ints, double t0, int nints,
	    double *chisq_out, double beta[MAX_PARS],
	    double alpha[MAX_PARS][MAX_PARS], int calc_deriv, int *ssel);
double eval_int_pos(PDecomp *pdd,
		    Event_Signal asig,  /* observed signals */
		    Event_Signal *bsig, /* fitted signals */
		    Interaction *ints, double t0, int nints);

/* in interpolate.c: */
int    nearest_grid_points(PDecomp *pdd,
			   int seg, double rin, double pin, double zin,
			   float *rdiff, float *pdiff, float *zdiff, int pos_out[8]);
int    interpolate(PDecomp *pdd,
		   int seg, double rin, double pin, double zin, Basis_Point *sig,
		   Basis_Point sigderiv[3], int calc_deriv, int *ssel);

/* in matinv.c: */
int    matinv(double *array, int norder, int dim);

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
/* routine to do overall initialization */
int decomp_init(char *basis_file_name, int set_quiet);
//int decomp(Event_Signal *asig);
//char *crys_intpts_2s(struct crys_intpts *x);
//void set_coal_dist(float d);
//int num_net(Event_Signal *asig);

#endif	/* _pdecomp_h */
