#ifndef _GEANT_SIGGEN_H
#define _GEANT_SIGGEN_H

#include <stdarg.h>
#include "pdecomp.h"

#define BASIS "geant_siggen_basis.dat"					// Basis file 
//#define MODE2 "Geant_Heather/ScanCs-24mmR-2mmColl.dat"		// Cs Geant 4 simulation data 2mm coll
//#define MODE2 "Geant_Heather/ScanCo-24mmR-2mmColl.dat"			// Co Geant 4 simulation data 2mm coll
#define MODE2 "Geant_Heather/ScanCs-24mmR-1mmColl.dat"			// Cs Geant 4 simulation data 1mm coll
//#define MODE2 "Geant_Heather/Scan.dat"				// Geant 4 test simulation data
#define NEVNTS 100

/*
	- For use with OakRidge/LBNL detector, mode2_struct defines the data structure for mode 2 GRETINA data i.e. the output of Heathers G4 code

	- For use with SIGMA, Marc_G4_Struct defines the struct for the output of Marc's G4 simulation, with SigGen_G4_Struct defining the parsed struct for use in SigGen

	- For both detectors, SigGen_Output defines the struct for the output of SigGen
*/

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

typedef struct {
  int SGeventm, SGnHits;					//event #, # of interactions in event
  int SGring, SGmodule, SGsector, SGslice, SGanneau;	//i.d. for ring, module, sector, slice and anneau
  float SGenergy;					//energy of interaction
  float SGx, SGy, SGz;					//xyz of interaction
  float SGtime;						//time of interaction
} Marc_G4_Struct;

typedef struct {
  int num;			//# of interactions in 1 event
  float tot_e;			//total energy of event
  struct {
    float x, y, z, e;		//xyz and energy of each interaction
  } intpts[MAX_INTPTS];
} SigGen_G4_Struct;

//int read_mode2(FILE *file, mode2_struct *m2);
int read_G4(FILE *file, SigGen_G4_Struct *m2);

#endif
