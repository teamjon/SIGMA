/* program to extract proper portion of measured signals to
   type expected by gdecomp

   Author:  D.C. Radford    Oct 2006
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "gdecomp.h"

#define SIG_FILE "/home/users/gretina/evtBuild/inBeam28.r1_eb.dat"

float gain[37] =
  {1.322816, 1.324088, 1.309852, 1.341286, 1.259293, 1.306283,
   1.346097, 0.652350, 1.323821, 1.255738, 1.285552, 1.297279,
   1.391064, 1.312160, 1.352038, 1.343364, 1.320638, 1.307491,
   1.364578, 1.326912, 1.345889, 1.310774, 1.282102, 1.324368,
   1.351193, 1.333868, 1.378939, 1.315090, 1.328648, 1.354061,
   1.320814, 1.319716, 1.305407, 1.303908, 1.299641, 1.292464, 1.529};

float delay0[37] =
  {19.64,21.04,18.01,19.22,19.55,19.14,19.14,21.78,18.40,18.45,21.18,18.48,
   19.59,20.96,18.80,19.71,20.81,20.40,21.80,23.19,19.31,21.80,19.88,19.88,
   21.98,22.77,19.20,19.81,19.16,20.45,17.63,18.75,17.52,17.97,16.68,17.37,
   0.0};
float delay1[37] =
  {0.28, 0.26, 1.51, 0.68, 0.68, 1.07, 1.47, 0.17, 1.58, 1.32,-0.16, 1.21,
   0.76,-0.12, 1.33, 0.85, 0.16, 0.03,-0.28,-0.58, 0.37,-0.57, 0.06, 0.19,
  -1.15,-0.96,-0.11,-0.55, 0.32,-0.14, 2.00, 1.50, 1.92, 1.46, 3.48, 2.23,
   0.0};

int t_cfd(int *buf)
{
  int i, int_length = 4, imax = 0, max_deriv = 0, delay = 4;
  int deriv[110], cfd[110];

  deriv[0] = 0;
  for (i = 0; i < int_length; i++) {
    deriv[0] += buf[i+int_length] - buf[i];
  }
  for (i = 1; i < 105 - 2*int_length; i++) {
    deriv[i] = deriv[i-1] +
       buf[i+2*int_length] - 2 * buf[i+int_length] + buf[i];
    if (max_deriv < deriv[i]) {
      /* imax = time step where current pulse has its maximum value */
      max_deriv = deriv[i];
      imax = i;
    }
  }
  for (i = 0; i < 105 - 2*int_length - delay; i++) {
    cfd[i] = deriv[i] - deriv[i+delay]/4;
  }
  for (i = imax + delay; i > 0; i--) {
    if (cfd[i] <= 0 && cfd[i+1] > 0) {
      /* interpolate zero crossing and return time to the nearest ns */
      return 10*i - 10*cfd[i]/(cfd[i+1] - cfd[i]);
      break;
    }
  }

  return 0;
}

int align_cfd(int buf[37][120])
{
  int   i, j, k = 1000, n = 0, m = 0, it;
  float s0 = 0.0f, s1 = 0.0f, t, tt;

  if (buf[36][107 ] < 400 || buf[36][107] > 630) return 0;  // CHANGEME: should put in terms of CC energy
  for (j=0; j<37; j++) {
    if (buf[j][107] > 100) {
      // s1 += ((float)t_cfd(buf[j]))/10.0 + delay1[j];
      // m++;
      i = t_cfd(buf[j]) + (int) (10.0f*delay1[j]);
      if (k > i) k = i;
      if (j<36) {
	s0 += delay0[j];
	n++;
      }
    }
  }
  if (n<1) return 0;
  t = ((float) k)/10.0f + s0/((float) n) - 16.0;
  // t = s1/((float) m) + s0/((float) n) - 16.0;
  if (t < 0) return 0;
  for (j=0; j<37; j++) {
    tt = t + delay1[j];
    if (tt <= 0.0f) continue;
    it = tt;
    tt -= (float) it;
    for (k=0; k<50; k++) {
      buf[j][k] = (1.0f - tt) * buf[j][k+it] + tt * buf[j][k+it+1];
    }
  }
  return 1;
}

int main(int argc, char **argv)
{
  unsigned short head[14], ibuf[4096], mbuf[4096], spec[37][4096];
  float fact;
  int   obuf[37][120];
  int   evt, i, j, k, nout=0;
  int   id, dsize, time = 0, tled, e, s;
  FILE  *file_in, *file_out, *file_mat;
  Event_Signal event;


  printf("Program convert_signals\n"
	 "Converts measured signal files to form expected by gdecomp\n\n");

  /* open the input and output files */
  printf("Reading signals from %s\n", SIG_FILE);
  if (!(file_in=fopen(SIG_FILE, "r"))) {
    printf("ERROR -- Cannot open file %s for reading.\n", SIG_FILE);
    return 1;
  }
  if (!(file_out=fopen(EVENT_FILE, "w"))) {
    printf("ERROR -- Cannot open file %s for writing.\n", EVENT_FILE);
    return 1;
  }
  if (!(file_mat=fopen("sig.mat", "w"))) {
    printf("ERROR -- Cannot open file sig.mat for writing.\n");
    return 1;
  }

  /* take first 50000 events */
  for (evt=0; evt<50000; evt++) {
    if (!(evt%100)) {printf(" %d\r", evt); fflush(stdout);}

    if (fread(head, sizeof(head), 1, file_in) != 1) {
      printf("ERROR - cannot read event number %d from file %s\n",
	     evt, SIG_FILE);
      return 1;
    }
    if (head[0] != 43690 || head[1] != 43690) {
      printf("\nBad header: %d %d\n", head[0], head[1]);
      break;
    }

    id = head[3] - 512;
    dsize = head[2] * 2 - 12;
    if (2*dsize > sizeof(ibuf)) {
      printf("\nBad buffer size: %d\n", dsize);
      break;
    }
    if (fread(ibuf, 2*dsize, 1, file_in) != 1) {
      printf("\nCannot read data: %d %d\n", head[2], dsize);
      break;
    }

    /* throw away data from detectors B and C for now */
    if (id > 40) continue;
    if (id > 35 && id < 39) continue;
    if (id == 39) id -= 3;
    // swap segments 28 & 34 for crystal A:
    if (id == 28) {
      id = 34;
    } else if (id == 34) {
      id = 28;
    }

    tled = head[4];
    tled = (tled << 16) + head[5];

    e = (head[9] & 63);
    e = (e << 16) + head[6];
    if (head[9] & (1<<6)) e -= (1<<22);
    if (id < 36) e = -e;

    if (abs(time - tled) > 2 && time != 0) {
      // new event is beginning
      if (align_cfd(obuf)) {
	memset(mbuf, 0, sizeof(mbuf));
	/* 3.4111 = CC energy/signal max,  obuf[36][110] = CC energy */
	fact = (float) obuf[36][110]/ 3.4111;
	for (j=0; j<36; j++) {
	  event.seg_energy[j] = (float) obuf[j][110] / (gain[j] * gain[36]);
	  for (k=0; k<50; k++) {
	    event.signal[j][k] = (float) obuf[j][k] / (gain[j] * fact);
	    mbuf[100*j + k] = 10000.0*event.signal[j][k];
	  }
	}
	event.total_energy = obuf[36][110] / gain[36];
	for (k=0; k<50; k++) {
	  event.signal[36][k] = -(float) obuf[36][k] / fact;
	  mbuf[100*36 + k] = 10000.0*event.signal[36][k];
	}
	fwrite(&event, sizeof(event), 1, file_out);
	fwrite(mbuf, sizeof(mbuf), 1, file_mat);
	nout++;
      }
      memset(obuf[0], 0, sizeof(obuf));
    }
    time = tled;

    for (i=0; i<108; i+=2) {
      obuf[id][i] = -(ibuf[i+1] & 2047);
      if (ibuf[i+1] & 2048) obuf[id][i] += 2048;
      if (id == 36) obuf[id][i] = -obuf[id][i];
      obuf[id][i+1] = -(ibuf[i] & 2047);
      if (ibuf[i] & 2048) obuf[id][i+1] += 2048;
      if (id == 36) obuf[id][i+1] = -obuf[id][i+1];
    }
    e = e >> 7;
    obuf[id][110] = e;
    // printf("id, e = %d %d\n", id, e);
    if (e > 100 && e < 4095) spec[id][e]++;

    s = 0;
    for (i=0; i<15; i++) {
      s += obuf[id][i];
    }
    s = (s + 7)/15;
    for (i=0; i<108; i++) {
      obuf[id][i] -= s;
    }

  }

  printf("\n  %d events processed, %d events output.\n\n", evt, nout);
 
  fclose(file_in);
  fclose(file_out);
  fclose(file_mat);
  if (!(file_mat=fopen("e.mat", "w"))) {
    printf("ERROR -- Cannot open file e.mat for writing.\n");
    return 1;
  }
  fwrite(spec[0], sizeof(spec), 1, file_mat);
  fclose(file_mat);  

  return 0;
}
