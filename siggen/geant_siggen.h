#ifndef _GEANT_SIGGEN_H
#define _GEANT_SIGGEN_H

#include <stdarg.h>
#include "pdecomp.h"

#define BASIS "geant_siggen_basis.dat"
//#define MODE2 "Geant_Heather/ScanCs-24mmR-2mmColl.dat"		// Cs Geant 4 simulation data 2mm coll
//#define MODE2 "Geant_Heather/ScanCo-24mmR-2mmColl.dat"			// Co Geant 4 simulation data 2mm coll
#define MODE2 "Geant_Heather/ScanCs-24mmR-1mmColl.dat"			// Cs Geant 4 simulation data 1mm coll
//#define MODE2 "Geant_Heather/Scan.dat"				// Geant 4 test simulation data
#define NEVNTS 100

typedef struct {
  int type;				/* Interaction type id -> temp set to 25 */
  int num;				/* # of interaction points from decomp, or # of nets on decomp error */
  float tot_e;				/* CC energy for the central contact selected for use in decomp (calibrated, and for 10MeV channels, includes DNL correction) */
  int nseg;				/* Number of segments */
  int time_step;			/* Time step used for output pulses */
  struct {
    float x, y, z, e;				/* here e refers to the fraction */
    int seg;					/* segment number hit */
    float seg_ener;				/* energy of the hit segment */
    float t_drift;				/* PC drift time */
  }intpts[MAX_INTPTS];
  float signal_mult[NUM_SIGS][TIME_STEPS_C];
} SigGen_Output;

typedef struct {
  int type;			//as of June2012: abcd5678
  int crystal_id;
  int num;			//# of interaction points from decomp, or # of nets on decomp error
  float tot_e;			//CC energy for the central contact selected for use in decomp (calibrated, and for 10MeV channels, includes DNL correction)
  int core_e[4];		//4 raw core energies from FPGA filter (uncalibrated)
  long long int timestamp;
  long long trig_time;		//not yet implemented
  float t0;
  float cfd;
  float chisq;
  float norm_chisq;
  float baseline;
  float prestep;		//avg trace value before step (baseline)
  float poststep;		//avg traces value after step (flat-top)
  int pad;			//non-0 with a decomp error, value gives error type
  struct {
    float x, y, z, e;		//here e refers to the fraction
    int seg;			//segment number hit
    float seg_ener;		//energy of the hit segment
  }intpts[MAX_INTPTS];
} mode2_struct;

#define MAX_SIM_GAMMAS 10

typedef struct {
  float e;
  float x, y, z;
  float phi, theta;
  float beta;
} g4Sim_emittedGamma;

static struct GlobalHeaderStruct {
  int type;            // Integer type value, indicating type of data to follow/
                       // Type values are assigned and documented at LBNL.
  int length;          // Integer length, in bytes, of the data payload to follow.
  long long timestamp; // 50MHz timestamp corresponding to the data to follow.
} GHd;

typedef struct {
  short iseg;		                         /* hit segment & integer cylindrical coordinates of grid point */
  float x, y, z;                                 /* cartesian coordinates of grid point */
  float t_drift;                                 /* PC drift time */
  float signal[NUM_SIGS][TIME_STEPS_C];          /* actual basis signals */
} Basis_Struct;

typedef struct {
  mode2_struct   *basis;                /* basis-signal data */
} MDecomp;

int read_mode2(FILE *file, mode2_struct *m2);

#endif
