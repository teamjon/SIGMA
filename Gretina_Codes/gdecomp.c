/* program to decompose event signals for GRETINA / GRETA.

   Author:  D.C. Radford    Aug 2004
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef _USING_ROOT
#include "root_includes.h"
#endif // _USING_ROOT
#include "gdecomp.h"

#define NEVENTS 1000     /* number of events to process, normally 10000 for timing tests */

long int lrintf(float x);

int quiet;               /* set to 1 to prevent output of diagnostic info */

int main(int argc, char **argv)
{
  double xfinal[MAX_PARS], yfinal[MAX_PARS], zfinal[MAX_PARS], efinal[MAX_PARS];
  double chisq, chisq2, esum;
  int    i, ii, j, k, t, event_number, found, photopk, events_proc = 0, nints;
  int    post_process = 0;
  short  mat[4096];
  char   in_fn[256] = EVENT_FILE, out_fn[256] = "output/test_gdecomp.out";
  FILE   *event_file, *file_out = 0, *mat_file = 0;

  Event_Signal  asig;                   /* observed event data */
  Event_Signal  bsig;                   /* calculated/fitted event data */
  Interaction   ints[2*MAX_SEGS];       /* interactions */
  int           nseg;                   /* number of segments that fired in the event */
  int           seg[TOT_SEGS];          /* IDs of the segments that fired */
  double        t0;                     /* time-zero parameter */
  Event_Signal  avg_obs_sig[TOT_SEGS];  /* averaged net=1 observed event data */
  Event_Signal  avg_fit_sig[TOT_SEGS];  /* averaged net=1 fitted event data */
  int           average_sigs = 0;       /* set to 1 to prevent output averaged observed and
					   fitted signals for single-segment (net=1) events */

  /* process possible command line options:
     -q = quiet (no debug output, used for timing tests)
     -a = write out averaged signals for net=1 events
     -p = perform postprocessing of interactions
     -i <input_file_name>
     -o <output_file_name>
  */
  quiet = 0;
  for (i=1; i<argc; i++) {
    if (*argv[i] == '-') {
      if (strstr(argv[i], "q")) quiet = 1;
      if (strstr(argv[i], "a")) {
	average_sigs = 1;
	memset(&avg_obs_sig[0].total_energy, 0, sizeof(Event_Signal)*TOT_SEGS);
	memset(&avg_fit_sig[0].total_energy, 0, sizeof(Event_Signal)*TOT_SEGS);
      }
      if (strstr(argv[i], "p")) post_process = 1;
      if (strstr(argv[i], "pp")) post_process = 2;
      if (strstr(argv[i], "i")) {
	strncpy(in_fn, argv[++i], sizeof(in_fn));
      } else if (strstr(argv[i], "o")) {
	strncpy(out_fn, argv[++i], sizeof(out_fn));
      }
    }
  }

  /* open event data file for reading */
  if (!(event_file=fopen(in_fn, "r"))) {
    printf("ERROR -- Cannot open file %s for reading.\n", in_fn);
    return 1;
  }
  /* open results output file for writing */
  if (!(file_out=fopen(out_fn, "w+"))) {
    printf("ERROR -- Cannot open file %s for writing.\n", out_fn);
    return 1;
  }
  if (!quiet || average_sigs) {
    /* open matrix file for writing */
    if (!(mat_file=fopen("output/out.mat", "w+"))) {
      printf("ERROR -- Cannot open file output/out.mat for writing.\n");
      return 1;
    }
  }

#ifdef _USING_ROOT
  printf("USING ROOT\n");
  printf("opening root file for output!\n");
  TTree::SetMaxTreeSize(10*Long64_t(2000000000)); // e.g. up to ~20GB
  TFile *fout_gdecomp = new TFile("fout_gdecomp.root","RECREATE");
  TTree *gdtree = new TTree("tout","");

  gdtree->Branch("found",&found,"found/I");
  //  gdtree->Branch("merge",&merge,"merge/I");
  gdtree->Branch("photopk",&photopk,"photopk/I");
  gdtree->Branch("asig",&asig,"total_energy/F:seg_energy[36]/F:signal[37][50]");
  gdtree->Branch("bsig",&bsig,"total_energy/F:seg_energy[36]/F:signal[37][50]");
  gdtree->Branch("ints",ints,"seg/I:pos/I:r/F:p/F:z/F:e/F:dr/F:dp/F:dz/F:de/F");
  gdtree->Branch("t0",&t0,"t0/D");
  gdtree->Branch("chisq",&chisq,"chisq/D");
  gdtree->Branch("chisq2",&chisq2,"chisq2/D");
  gdtree->Branch("event_number",&event_number,"event_number/I");
  gdtree->Branch("xfinal",xfinal,"xfinal[found]/D");
  gdtree->Branch("yfinal",yfinal,"yfinal[found]/D");
  gdtree->Branch("zfinal",zfinal,"zfinal[found]/D");
  gdtree->Branch("efinal",efinal,"efinal[found]/D");
#endif // _USING_ROOT

  /* read in basis signal data */
  if (read_basis(BASIS_FILE)) return 1;

  /* initialize arrays for the grid-search algorithm */
  pclock(1);
  if (grid_init()) return 1;
  pclock(1);

  /* start of event decomposition loop over events */

  for (event_number = 0; event_number < NEVENTS; event_number++) {
    if (!quiet) printf("Event %d\n", event_number); fflush(stdout);

    /* read in an event */
    if (!fread(&asig, sizeof(Event_Signal), 1, event_file)) {
      printf("No more events in input file?\n"
	     "     %d events read, %d events processed.\n", event_number, events_proc);
      break;
    }
    if (asig.total_energy < 20.0) continue;

    /* set bad segment signals to zero */
    for (i=0; bad_segs[i]>=0; i++) {
      ii = bad_segs[i];
      for (t=0; t<TIME_STEPS; t++) asig.signal[ii][t] = 0.0f;
      asig.seg_energy[ii] = 0.0f;
    }

    /* determine which segments have net energy */
    nseg = 0;
    esum = 0.0;
    for (i=0; i<TOT_SEGS; i++) {
      seg[i] = 0;
      if (asig.seg_energy[i] > 30.0) {
	esum += asig.seg_energy[i];
	seg[nseg++] = i;
      }
    }
    if (nseg == 0) continue;

    if (!quiet) {
      printf("-----------\n"
	     " Event # %d,  total energy = %5.0fkeV\n"
	     "  %d segments hit: ",
	     event_number, asig.total_energy, nseg);
      for (i=0; i<nseg; i++) {
	printf(" %2d: %4.0fkeV ", seg[i], asig.seg_energy[seg[i]]);
      }
      printf("\n-----------\n");
    }

    if (fabs(esum - asig.total_energy) > 30.0) {
      if (!quiet) printf("***Error***  Sum energy (%.0f) != CC energy (%.0f)\n",
			 esum, asig.total_energy);
      continue;
    }

    if (nseg > MAX_SEGS) {
      if (!quiet) printf("ACKK! Too many segments! nseg =%d\n", nseg);
      continue;
    }

    /* clear out arrays */
    memset(&bsig, 0, sizeof(Event_Signal));
    memset(&ints, 0, 2*nseg*sizeof(Interaction));

    /* branch according to number of hit segments */
    if (nseg == 1) {
      //continue;
      nints = decompose_1(asig, &bsig, seg[0], ints, &t0, &chisq,
			  0, 1, 1, 1, 1, 1, 1, 1, 0.1);

    } else if (nseg == 2) {
      //continue;
      nints = decompose_n(asig, &bsig, nseg, seg, 1, ints, &t0, &chisq);

    } else {
      //continue;
      nints = decompose_n(asig, &bsig, nseg, seg, 1, ints, &t0, &chisq);
    }
    if (chisq > 99.9) {
      printf("Bad Chi-squared = %f for event number %d\n", chisq, event_number);
      continue;
    }
    events_proc++;

    if (post_process) {  /* post_process the results */
      found = postprocess_events(ints, nints, asig.total_energy, post_process,
				 COAL_DIST_DEFAULT, xfinal, yfinal, zfinal, efinal);

      /* write results to file_out */
      photopk = 0; // no longer used; flag is available for reuse?
      chisq2 = chisq / (float) (MEAS_SEGS * TIME_STEPS); // approx. degrees of freedom
      chisq2 /= (4.0/asig.total_energy)*(4.0/asig.total_energy);  // noise ~ 4 keV
      fprintf(file_out,
	      " %2d %2d %9.2f  %6.2f %24.6f %9.4f %4d\n",
	      found, photopk, asig.total_energy, t0, chisq,
	      chisq2, event_number);
      for (i=0; i<found; i++) {
	fprintf(file_out,
		" %7.2f %7.2f %7.2f %7.1f\n",
		xfinal[i], yfinal[i], zfinal[i], efinal[i]);
      }

#ifdef _USING_ROOT
      gdtree->Fill();
      if (event_number%1000 == 0){
        fout_gdecomp->cd();
	// printf("AutoSaving Tree....");
        gdtree->AutoSave("SaveSelf");
	// printf("ok\n");
      }
#endif // _USING_ROOT

      } else {  /* do *not* post_process the results */
      /* write results to file_out */
      chisq2 = chisq / (float) (MEAS_SEGS * TIME_STEPS); // approx. degrees of freedom
      chisq2 /= (4.0/asig.total_energy)*(4.0/asig.total_energy);  // noise ~ 4 keV
#ifdef W_SEGS
      fprintf(file_out,
	      " %2d %2d %9.2f  %6.2f %24.6f %9.4f %4d\n segs:",
	      nseg, nints, asig.total_energy, t0, chisq, chisq2, event_number);
      for (i=0; i<nseg; i++) {
	fprintf(file_out,
		" %3d %7.1f;", seg[i], asig.seg_energy[seg[i]]);
      }
      fprintf(file_out, "\n");
#else
      fprintf(file_out,
	      " %2d %2d %9.2f  %6.2f %24.6f %9.4f %4d\n",
	      nseg, nints, asig.total_energy, t0, chisq, chisq2, event_number);
#endif
      for (i=0; i<nints; i++) {
	fprintf(file_out,
		" %6d %7.2f %7.2f %7.2f %7.1f\n",
		ints[i].seg, ints[i].r, ints[i].p, ints[i].z, ints[i].e*asig.total_energy);
      }
    }

    if (!quiet) {
      /* write signals to matrix file */
      if (event_number < 4096) {
	for (i=0; i<4096; i++) {
	  mat[i] = 0;
	}
	for (i=0; i<MEAS_SEGS; i++) {
	  for (j=0; j<TIME_STEPS; j++) {
	    mat[50*i + j] = lrintf(asig.signal[i][j] * 10000.0f);
	    mat[2000 + 50*i + j] = lrintf(bsig.signal[i][j] * 10000.0f);
	  }
	}
	fwrite(mat, 8192, 1, mat_file);
      }
    }

    if (average_sigs && nseg == 1) {
      avg_obs_sig[seg[0]].total_energy += asig.total_energy;
      avg_fit_sig[seg[0]].total_energy += asig.total_energy;
      for (i=0; i<MEAS_SEGS; i++) {
	avg_obs_sig[seg[0]].seg_energy[i] += asig.seg_energy[i];
	avg_fit_sig[seg[0]].seg_energy[i] += asig.seg_energy[i];
	for (j=0; j<TIME_STEPS; j++) {
	  avg_obs_sig[seg[0]].signal[i][j] += asig.signal[i][j];
	  avg_fit_sig[seg[0]].signal[i][j] += bsig.signal[i][j];
	}
      }
    }

  }
#ifdef _USING_ROOT
  fout_gdecomp->cd();
  gdtree->Write();
  fout_gdecomp->Close();
#endif // _USING_ROOT

  /* pclock(3); */
  pclock(1);
  free(basis);
  for (i=0; i<SSEG; i++) {
    free(ags1[i].da);
    free(ags1[i].db);
  }

  if (average_sigs) {
      /* write averaged signals to matrix file */
    for (i=0; i<4096; i++) {
      mat[i] = 0;
    }
    for (k=0; k<TOT_SEGS; k++) {
      for (i=0; i<MEAS_SEGS; i++) {
	for (j=0; j<TIME_STEPS; j++) {
	  mat[       50*i + j] = lrintf(avg_obs_sig[k].signal[i][j] * 10000000.0f/avg_obs_sig[k].total_energy);
	  mat[2000 + 50*i + j] = lrintf(avg_fit_sig[k].signal[i][j] * 10000000.0f/avg_obs_sig[k].total_energy);
	}
      }
      fwrite(mat, 8192, 1, mat_file);
    }
  }

  return 0;
} /* main */

/* ================================================================ */

/* This is just a quick time profiling routine to see where the
   decomposition algorithm is spending all it's time.
   Author:  D.C. Radford    May 2003
*/

void pclock(int print_flag)
{
  static int s = 0, t1 = 0, t2;
  /* usage:
     pclock(0) on entry to timed code
     pclock(1) on exit from timed code, to print time interval
     pclock(2) on exit from timed code, to add time interval to sum
     pclock(3) to print out sum of time intervals

     Time increments seem to be in 1/100 seconds
     So for multiple short time intervals, statistical fluctuations may be important.
  */

  t2 = (int) clock();
  if (print_flag == 1) printf("*** clock, interval: %10d %10d\n", t2, t2 - t1);
  if (print_flag == 2) s += t2 - t1;
  if (print_flag == 3) {
    printf("*** clock sum: %10d\n", s);
    s = 0;
  }
  t1 = t2;
} /* pclock */
