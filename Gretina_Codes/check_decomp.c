/* check_decomp.c
 *
 * by John Pavan, David Radford
 *
 * checks to see if there are reasonable agreement between simulated
 * input and the corresponding output from the decomposition program
 */

/* --- includes --- */
#include <stdio.h>
#include <strings.h>
#include <math.h>

/* --- defines --- */
#define MAX_INTERACTIONS 10
#define MAX_EVENTS 5000

/* --- structs --- */

struct event_out{
  int n;     // number of interactions found
  int event; // which event is it?
  float tot_e; // total energy for the event
  float chisqr;  // chisqr of the event
  float chisqrdof; // chisqr per degree of freedom
  float x[MAX_INTERACTIONS];
  float y[MAX_INTERACTIONS];
  float z[MAX_INTERACTIONS];
  float e[MAX_INTERACTIONS];
  int pp; // photopeak flag
};

struct event_in{
  int n;
  float tot_e;
  float x[MAX_INTERACTIONS];
  float y[MAX_INTERACTIONS];
  float z[MAX_INTERACTIONS];
  float e[MAX_INTERACTIONS];
};

int events_in = 0;
int events_out = 0;


/* --- globals --- */

struct event_out oevents[MAX_EVENTS];
struct event_in ievents[MAX_EVENTS];
int histogram[4096];
int xhist[4096];
int yhist[4096];
int zhist[4096];
int rhist[4096];
int poshist[4096];
int ehist[4096];

/* --- fuction declarations --- */

int ReadInput(char *sfilename);
int ReadOutput(char *sfilename);
void DoComparison(char *sfilename,char *hisfilename);

/* --- functions --- */
void DoComparison(char *sfilename,char *hisfilename)
{
  float rms, rmsold, r1, r2, RMSpos, RMSx, RMSy, RMSz, RMSe, RMSr, d;
  int   i, j, k, l, m, n, *best_perturbation;
  int   perturbations1[1][1] = {{0}};
  int   perturbations2[2][2] = {{0,1},{1,0}};
  int   perturbations3[6][3] = {{0,1,2},{1,0,2},{1,2,0},{2,1,0},{2,0,1},{0,2,1}};
  int   perturbations4[24][4] = {{0,1,2,3},{1,0,2,3},{1,2,0,3},{2,1,0,3},{2,0,1,3},{0,2,1,3},
				 {0,1,3,2},{1,0,3,2},{1,2,3,0},{2,1,3,0},{2,0,3,1},{0,2,3,1},
				 {0,3,1,2},{1,3,0,2},{1,3,2,0},{2,3,1,0},{2,3,0,1},{0,3,2,1},
				 {3,0,1,2},{3,1,0,2},{3,1,2,0},{3,2,1,0},{3,2,0,1},{3,0,2,1}};
  FILE  *output;
  

  /* --- open output file for event listings --- */
  if (!(output = fopen(sfilename,"w"))) {
    printf("Could not open %s for writing\n", sfilename);
    return;
  }

  /* --- initialize histograms --- */
  for (i = 0; i < 4096; i++) {
    histogram[i] = xhist[i] = yhist[i] = zhist[i] = rhist[i] = poshist[i] = ehist[i] = 0;
  }

  /* --- scan through events --- */
  for (i = 0; i < events_out; i++) {
    n = oevents[i].chisqrdof * 100.0 + 0.5;
    if (n >= 0 && n < 4096) histogram[n]++;

    l = oevents[i].event;
    if (l >= events_in) continue;
    /* --- only compare if they have the same number of interactions --- */
    if (ievents[l].n == oevents[i].n) {
      rmsold = 10000000;
      best_perturbation = 0;
      /* --- now we want to try and match up the interactions so that
	 --- the rms difference is minimized --- */
      /* --- make perturbations of the output events to try and match 
	 --- up with the input events --- */
      switch (oevents[i].n) {
      case 1:
	/*--- 1 possible combination --- */
	best_perturbation = perturbations1[0];
	break;
      case 2:
	/* --- 2 possible combinations --- */
	for (k = 0; k < 2; k++) {
	  rms = 0;
	  for (j = 0; j < ievents[l].n; j++) {
	    rms += sqrt(pow(ievents[l].x[j] - oevents[i].x[perturbations2[k][j]],2) +
			pow(ievents[l].y[j] - oevents[i].y[perturbations2[k][j]],2) +
			pow(ievents[l].z[j] - oevents[i].z[perturbations2[k][j]],2));
	  }
	  rms = sqrt(rms);
	  if (rms < rmsold) {
	    rmsold = rms;
	    best_perturbation = perturbations2[k];
	  }
	}
	break;
      case 3:
	/* --- 6 possible combinations --- */
	for (k = 0; k < 6; k++) {
	  rms = 0;
	  for (j = 0; j < ievents[l].n; j++) {
	    rms += sqrt(pow(ievents[l].x[j] - oevents[i].x[perturbations3[k][j]],2) +
			pow(ievents[l].y[j] - oevents[i].y[perturbations3[k][j]],2) +
			pow(ievents[l].z[j] - oevents[i].z[perturbations3[k][j]],2));
	  }
	  rms = sqrt(rms);
	  if (rms < rmsold) {
	    rmsold = rms;
	    best_perturbation = perturbations3[k];
	  }
	}
	break;
      case 4:
	/* --- 24 possible combinations --- */
	for (k = 0; k < 24; k++) {
	  rms = 0;
	  for (j = 0; j < ievents[l].n; j++) {
	    rms += sqrt(pow(ievents[l].x[j] - oevents[i].x[perturbations4[k][j]],2) +
			pow(ievents[l].y[j] - oevents[i].y[perturbations4[k][j]],2) +
			pow(ievents[l].z[j] - oevents[i].z[perturbations4[k][j]],2));
	  }
	  rms = sqrt(rms);
	  if (rms < rmsold) {
	    rmsold = rms;
	    best_perturbation = perturbations4[k];
	  }
	}
	break;
      default: 

	/* -- do nothing --- */
	printf("Error, out of range of interactions \n");
	break;
      }

      /* --- store relevant information in histograms --- */
      RMSpos = RMSx = RMSy = RMSz = RMSe = RMSr = 0;
      for (k = 0; k < ievents[l].n; k++) {
	j = best_perturbation[k];
	/* --- x hist --- */
	d = ievents[l].x[k] - oevents[i].x[j];
	RMSx += d*d;
	n = 200.5 + (d * 10.0);
	if (n >= 0 && n < 4096) xhist[n]++;
	/* --- y hist --- */
	d = ievents[l].y[k] - oevents[i].y[j];
	RMSy += d*d;
	n = 200.5 + (d * 10.0);
	if (n >= 0 && n < 4096) yhist[n]++;
	/* --- z hist --- */
	d = ievents[l].z[k] - oevents[i].z[j];
	RMSz += d*d;
	n = 200.5 + (d * 10.0);
	if (n >= 0 && n < 4096) zhist[n]++;
	/* --- pos hist --- */
	r1 = sqrt(pow(ievents[l].x[k],2) + pow(ievents[l].y[k],2) + pow(ievents[l].z[k],2));
	r2 = sqrt(pow(oevents[i].x[j],2) + pow(oevents[i].y[j],2) + pow(oevents[i].z[j],2));
	d = r1 - r2;
	RMSpos += d*d;
	n = 200.5 + (d * 10.0);
	if (n >= 0 && n < 4096) poshist[n]++;
	/* --- r hist --- */
	r1 = sqrt(pow(ievents[l].x[k],2) + pow(ievents[l].y[k],2));
	r2 = sqrt(pow(oevents[i].x[j],2) + pow(oevents[i].y[j],2));
	d = r1 - r2;
	RMSr += d*d;
	n = 200.5 + (d * 10.0);
	if (n >= 0 && n < 4096) rhist[n]++;
	/* --- e hist --- */
	d = ievents[l].e[k] - oevents[i].e[j];
	RMSe += d*d;
	n = 200.5 + d;
	if (n >= 0 && n < 4096) ehist[n]++;
      }
      RMSpos = sqrt(RMSpos/(float)ievents[l].n);
      RMSx = sqrt(RMSx/(float)ievents[l].n);
      RMSy = sqrt(RMSy/(float)ievents[l].n);
      RMSz = sqrt(RMSz/(float)ievents[l].n);
      RMSe = sqrt(RMSe/(float)ievents[l].n);
      RMSr = sqrt(RMSr/(float)ievents[l].n);

      /* --- now we need to start keeping track of stuff --- */
      fprintf(output,"----result----\n");
      fprintf(output,"%7d %7d %7.1f %7.1f %7.1f\n",
	      oevents[i].event,oevents[i].n,ievents[l].tot_e,oevents[i].tot_e,RMSr);
      fprintf(output,"%7.3f %7.3f %7.3f %7.3f %7.3f\n",
	      RMSx,RMSy,RMSz,RMSpos,RMSe);
      fprintf(output,"----input----\n");
      fprintf(output,"%d %7.3f\n",ievents[l].n,ievents[l].tot_e);
      for (m = 0; m < ievents[l].n; m++) {
	fprintf(output,"%7.3f %7.3f, %7.3f %7.3f\n",ievents[l].x[m],
		ievents[l].y[m],ievents[l].z[m],ievents[l].e[m]);
      }
      fprintf(output,"----output----\n");
      fprintf(output,"%d %7.3f\n",oevents[i].n,oevents[i].tot_e);
      for (m = 0; m < ievents[l].n; m++) {
	fprintf(output,"%7.3f %7.3f, %7.3f %7.3f\n",oevents[i].x[m],
		oevents[i].y[m],oevents[i].z[m],oevents[i].e[m]);
      }
      fprintf(output,"\n\n");
    }
  }
  fclose(output);

  /* --- open output file for histograms --- */
  if ((output = fopen(hisfilename, "w")) != NULL) {
    fwrite(histogram, sizeof(int), 4096, output);
    fwrite(xhist,     sizeof(int), 4096, output);
    fwrite(yhist,     sizeof(int), 4096, output);
    fwrite(zhist,     sizeof(int), 4096, output);
    fwrite(rhist,     sizeof(int), 4096, output);
    fwrite(poshist,   sizeof(int), 4096, output);
    fwrite(ehist,     sizeof(int), 4096, output);
    fclose(output);
    printf("Done writting histograms to %s.\n",hisfilename);
  } else {
    printf("Error opening %s for output.\n",hisfilename);
  }
}

int ReadOutput(char *sfilename)
{
  float a, b, c, d, e;
  int   i, j, k, l;
  char  dummystr[200];
  FILE  *input;

  if ((input = fopen(sfilename,"r")) != NULL) {
    while (fgets(dummystr,200,input)) {
      if ((sscanf(dummystr,"%d %d %f %d %f %f %d", &i, &j, &a, &k,&b,&c,&l)) == 7) {
	oevents[events_out].n = i;
	oevents[events_out].pp = j;
	oevents[events_out].tot_e = a;
	oevents[events_out].chisqr = b;
	oevents[events_out].chisqrdof = c;
	oevents[events_out].event = l;
	for (j = 0; j < i; j++) {
	  fgets(dummystr,200,input);
	  if (sscanf(dummystr,"%f %f %f %f", &b,&c,&d,&e) == 4) {
	    oevents[events_out].x[j] = b;
	    oevents[events_out].y[j] = c;
	    oevents[events_out].z[j] = d;
	    oevents[events_out].e[j] = e * a;
	  }
	}
	events_out++;
      }
    }
    
    fclose(input);
    return(1);
  } else {
    printf("Error opening file %s\n",sfilename);
    return(0);
  }
}

int ReadInput(char *sfilename)
{
  float a, b, c, d, e, f;
  int   i, j;
  char  dummystr[200];
  FILE *input;

  if ((input = fopen(sfilename,"r")) != NULL) {
    while (fgets(dummystr,200,input)) {
      if (sscanf(dummystr,"%d %f",&i,&a) == 2) {
	/* --- it is a new event --- */
	/* --- otherwise it is a comment or something is wrong --- */
	ievents[events_in].n = i;
	ievents[events_in].tot_e = a;
	for (j = 0; j < i; j++) {
	  if (fgets(dummystr,200,input)) {
	    if ((sscanf(dummystr,"%f %f %f %f %f %f",
			&a,&b,&c,&d,&e,&f)) ==6) {
	      if (j >= MAX_INTERACTIONS) continue;
	      if (d < 10.0) {
		j--;
		ievents[events_in].n = --i;
		continue;
	      }
	      ievents[events_in].x[j] = a;
	      ievents[events_in].y[j] = b;
	      ievents[events_in].z[j] = c;
	      ievents[events_in].e[j] = d;
	    }
	  } else {
	    /* --- we are out of file --- */
	    printf("Ran out of input file in an unexpected way\n");
	    fclose(input);
	    return(0);
	  }
	}
	events_in++;
      }
    }
    fclose(input);
    return(1);
  } else {
    printf("Error reading %s\n",sfilename);
    return(0);
  }
}

int main(int argn, char **argc)
{
  if (argn < 5) {
    printf("Usage: check_decomp sim_input_fn decomp_output_fn checker_output_fn output_spn_fn\n");
    return 1;
  }
  if (ReadInput(argc[1]) && ReadOutput(argc[2])) {
    DoComparison(argc[3],argc[4]);
  } else {
    printf("error reading names of files\n");
  }
  return 0;
}
