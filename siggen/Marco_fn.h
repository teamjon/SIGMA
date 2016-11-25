#include <stdlib.h>
#include <math.h>

/*=== Function ===============================================*/

double **Allocate2DSamples(int ny, int nx)

/*--- Description --------------------------------------------//

Allocates a 2 dim array of doubles.
Returns the pointers to the arrays on success, 0 on error.
The array of all double is consecutive.

//------------------------------------------------------------*/
{
  int i;
  double **yp=(double**)calloc(ny,sizeof(void*));
  double *dp=(double*)calloc(nx*ny,sizeof(double));
  for(i=0; i<ny; i++) yp[i]= &dp[i*nx];
  return yp;
}


/*=== Function ===============================================*/

int Free2DSamples(double **data)

/*--- Description --------------------------------------------//

Free the 2dim. array of doubles.
Returns 0 if data is a 0 pointer or 1 if data was freed.

//------------------------------------------------------------*/
{
  if(!data) return 0;
  free(data[0]);
  free(data);
  return 1;
}

//------------------------------------------------------------*/

/*=== Function ===============================================*/

double IntegrateSamples(int n, double *in, int start, int stop)

/*--- Description --------------------------------------------//

Return the sum of array elements within of the range start
and and stop.

//------------------------------------------------------------*/
{
  int i;
  double sum=0;
  if(start<0) start=0;
  if(stop>=n) stop=n-1;
  for(i=start;i<=stop;i++) sum+=in[i];
  return sum;
}


/*=== Function ===============================================*/

int DifferenceOfSamples(int n, double *in, double *out)

/*--- Description --------------------------------------------//

Differentiate the array in (size n) and store it in out.

WARNING: The first element of out is always the first element
of in, hence contains the information about the offset and
should be ignored for most applications.

The function returns n.

//------------------------------------------------------------*/
{
  int i;
  double last=0;
  double next=0;
  for(i=0; i<n; i++) next=in[i],out[i]=next-last,last=next;
  return n;
}


/*=== Function ===============================================*/

double FindLowerBoundOfSamples(int n, double *in, double low, double max)

/*--- Description --------------------------------------------//

Finds the position in the array in of size n, where the integral
reaches low as a fraction of the total integral of the array (no
normalization required). It requires that this position is bellow
max. It works with positive and negative pulses.

The exact position is interpolated linearly.

If array in is the derivative of a continously increasing
function then the return value is the time where the pulse crosses
low as a fraction of the total pulse height.

//------------------------------------------------------------*/
{
  if(max<0) max=0; if(max>n) max=n;
  int kk, sig=1;
  double ll=0,delta,a,v1,v2, root, pp;
  double igr=IntegrateSamples(n,in,0,n);
  if(igr==0){ return 0; };
  if(igr<0){ sig=-1; igr=-igr; }
  low*=igr;

  //lets go backwards, first
  for(kk=n-1; kk>max; kk--) igr-=sig*in[kk];
  for(;kk>=0;kk--){ igr-=sig*in[kk]; if(igr<low ) break; }
  if(kk>max-1) return max;
  if(kk<1) return 0;
  igr+=in[kk]; //move one up again for simpler interpolation algorithm

  v1=sig*in[kk-1]; v2=sig*in[kk];

  //interpolation (required condition: root>0, fabs(pp)<1, delta!=0 or v1!=0)
  delta=v2-v1; a=low-igr+v2; root=v1*v1+2.0*a*delta; pp=(sqrt(root)-v1)/delta;
  if( (delta==0 || fabs(pp)>1) && v1>a) ll=kk-1+a/v1; //approx
  else if(delta!=0 && root>=0 && fabs(pp)<1 ) ll=kk-1+pp; //normal value
  else ll=kk; //no other option
  //printf("kk %d ll %g min %g delta %g v1 %g v2 %g a %g\n",kk,ll,max,delta,v1,v2,a);
  return ll;
}


/*=== Function ===============================================*/

double FindUpperBoundOfSamples(int n,double *in, double hig, double min)

/*--- Description --------------------------------------------//

Finds the position in the array in of size n, where the integral
reaches hig as a fraction of the total integral of the array (no
normalization required). It requires that this position is above
min. It works with positive and negative pulses.

The exact position is interpolated linearly.

If array in is the derivative of a continously increasing
function then the return value is the time where the pulse crosses
hig as a fraction of total the pulse height.

//------------------------------------------------------------*/
{
  int kk, sig=1;
  if(min<0) min=0; if(min>n) min=n;
  double ll=0,delta,a,v1,v2, root, pp;
  double igr=IntegrateSamples(n,in,0,n);
  if(igr==0){ return 0; }
  if(igr<0){ sig=-1; igr=-igr; }
  hig*=igr; igr=0;

  //lets go forward
  for(kk=0; kk<min; kk++) igr+=sig*in[kk];
  for(; kk<n; kk++){ igr+=sig*in[kk]; if(igr>hig) break; }
  if(kk<min+1) return min;
  if(kk>n-1) return n-1;
  v1=sig*in[kk-1]; v2=sig*in[kk];

  //interpolation (required condition: root>0, fabs(pp)<1, delta!=0 or v1!=0)
  delta=v2-v1; a=hig-igr+v2; root=v1*v1+2.0*a*delta; pp=(sqrt(root)-v1)/delta;
  if( (delta==0 || fabs(pp)>1) && v1>=a) ll=kk-1+a/v1; //approx
  else if(delta!=0 && root>=0 && fabs(pp)<1 ) ll=kk-1+pp; //normal value

  else ll=kk; //no other option
  //printf("kk %d ll %g min %g delta %g v1 %g v2 %g a %g\n",kk,ll,min,delta,v1,v2,a);
  return ll;
}


/*=== Function ===============================================*/

double CenterOfGravityOfSamples(int n, double *in, int start, int stop)

/*--- Description --------------------------------------------//

Returns the center of gravity (cog) of the range start and stop.
All element must be greater than 0. there is no check within
of the function that this is fulfilled.

If the sum of the given range is 0 the function returns 0.


//------------------------------------------------------------*/

// data  must be >0

{
  int i;
  double isum=0;
  double sum=0;
  for(i=start;i<=stop;i++) isum+=in[i]*i,sum+=in[i];
  if(sum!=0) return isum/sum;
  return 0;
}


/*=== Function ===================================================*/

double GConCalculateDriftTime(int slen, double *psig, double *ssig,
  double low, double hig, double delta, double *ostart)

/*--- Description ------------------------------------------------//

Find the drift time between the point contact signal (psig) and the
segment signal (ssig) both of length (slen), with the limits low and
hig (in percent) using a shaping of sha and interpolate around delta
(in percent too).

ostart is filled with the start value if provided. It returns the
drift time.

Parameters used for SIGMA work low=0.02, hig=0.92, delta=0.03.

//----------------------------------------------------------------*/
{
  static int gcon_CDTlen=0;      //static: only initialized at first call
  static double **gcon_CDTarr=0; //static: only initialized at first call
  if(slen>gcon_CDTlen) //we only allocate when necessary
  {
    if(gcon_CDTarr) Free2DSamples(gcon_CDTarr);
    gcon_CDTarr=Allocate2DSamples(2,slen);
    gcon_CDTlen=slen;
  }
  int j=0;
  double start1=0, start2=0, start=0, stop=0, stop1=0, stop2=0;

  DifferenceOfSamples(slen,psig,gcon_CDTarr[0]);//central contact
  DifferenceOfSamples(slen,ssig,gcon_CDTarr[1]);//max segment contact
  stop2=FindUpperBoundOfSamples(slen,gcon_CDTarr[0],hig-delta,0);
  stop1=FindUpperBoundOfSamples(slen,gcon_CDTarr[0],hig,stop2);
  stop=stop1+(stop1-stop2)*(1-hig)/delta;
  start2=FindLowerBoundOfSamples(slen,gcon_CDTarr[1],low+delta, stop2);
  start1=FindLowerBoundOfSamples(slen,gcon_CDTarr[1],low, start2);
  start=start1-(start2-start1)*low/delta;

  if(ostart) *ostart=start;

  return CenterOfGravityOfSamples(slen, gcon_CDTarr[0], start, stop)-start;
}
