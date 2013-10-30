#ifndef CLANG
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cmath>
#endif
#ifdef CLANG
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#endif

#include <SpiceUsr.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv2.h>

#ifndef CLANG
using namespace std;
#endif

//////////////////////////////////////////////////////////////////////////////////
//SPICE INFO & MACROS
//////////////////////////////////////////////////////////////////////////////////
#define ABCORR "NONE"
#define FRAME  "ECLIPJ2000"

//////////////////////////////////////////////////////////////////////////////////
//MACROS
//////////////////////////////////////////////////////////////////////////////////
//PHYSICAL AND ASTRONOMICAL CONSTANTS
#define MSUN (GSL_CONST_MKSA_SOLAR_MASS) //kg
#define RSUN (6.96342E8) //m
#define LSUN (3.846E26) //W
#define AU (GSL_CONST_MKSA_ASTRONOMICAL_UNIT) //m
#define HOURS (3600.0) //s
#define YEAR (365.25*GSL_CONST_MKSA_DAY) //s
#define MEARTH (5.9736E24) //kg
#define REARTH (6.371E6) //m
#define GCONST (GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT) // m^3 / (kg s^2)
#define CLIGHT (GSL_CONST_MKSA_SPEED_OF_LIGHT) // m / s

#define QPR 1.0 //Radiation pressure coefficient

#define PI M_PI
#define PI2 (PI*PI)
#define D2R PI/180

//////////////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES
//////////////////////////////////////////////////////////////////////////////////
typedef double real;
typedef void* params;
typedef FILE* file;

real UL,UM,UT,GPROG;
real ALPHA;
real RHODUST;

//ELEMENTS
/*
  0 - RP: Distance at periapse
  1 - ECC: Eccentricity (rad)
  2 - INC: Inclination (rad)
  3 - LNODE: Longitude of the ascending node
  4 - ARGP: Argument of the periapse
  5 - M0: Initial mean anomaly
  6 - T0: Epoch
  7 - MU: G*Mcentral
  8 - PER: Time scale (period in case of elliptic orbits)
  9 - TPER: Time of periapse passage
  10 - AXISZ: Angle inclination of axis respect to orbit normal
  11 - AXISXY: Acimutal angle of body at periapse (x-axis: sun-periapse)
 */
#define NELEMS 12

enum ElementsEnum {RP,ECC,INC,LNODE,ARGP,M0,T0,MU,
		   TPER,PER,AXISZ,AXISXY,ENDELEM};
char Elements[][10]={"RP","ECC","INC","LNODE","ARGP","M0","T0","MU",
		     "TPER","PER","AXISZ","AXISXY"};

//STATE VECTOR
enum StateEnum {XPOS,YPOS,ZPOS,XVEL,YVEL,ZVEL,ENDSTATE};
char State[][10]={"XPOS","YPOS","ZPOS","XVEL","YVEL","ZVEL"};

//ROTATION MATRIX
enum Axis {AXIS,XAXIS,YAXIS,ZAXIS};

//////////////////////////////////////////////////////////////////////////////////
//ROUTINES
//////////////////////////////////////////////////////////////////////////////////
void fprintf_vec(FILE* stream,const char* fmt,const double x[],int end,bool qlabel=true,int ini=0)
{
  for(int i=ini;i<end;i++){
    if(qlabel) fprintf(stream,"%d:",i);
    fprintf(stream,fmt,x[i]);
  }
  fprintf(stream,"\n");
}

void fprintf_mat(FILE* stream,const char* fmt,const double x[][3])
{
  for(int i=0;i<3;i++){
    fprintf(stream,"|");
    for(int j=0;j<3;j++)
      fprintf(stream,fmt,x[i][j]);
    fprintf(stream,"|\n");
  }
}

void fprintf_elem(FILE *stream,const char* fmt,const double e[],bool qname=false)
{
  double elem;
  ElementsEnum ie;
  for(int i=RP;i!=ENDELEM;i++){
    ie=(ElementsEnum)i;
    if(qname) fprintf(stream,"%s:",Elements[ie]);
    fprintf(stream,fmt,e[ie]);
  }
  fprintf(stream,"\n");
}

void fprintf_state(FILE *stream,const char* fmt,const double x[],bool qname=false)
{
  double state;
  StateEnum ix;
  for(int i=XPOS;i!=ENDSTATE;i++){
    ix=(StateEnum)i;
    if(qname) fprintf(stream,"%s:",State[ix]);
    fprintf(stream,fmt,x[ix]);
  }
  fprintf(stream,"\n");
}

int setUnits(real ul,real um,real ut,real G)
{
  real ratio=(G/GCONST);
  if(ut==0){
    UL=ul;UM=um;GPROG=G;
    UT=sqrt(ratio*ul*ul*ul/um);
  }
  else if(ul==0){
    UM=um;UT=ut;GPROG=G;
    UL=1/pow(ratio/(um*ut*ut),1.0/3);
  }
  else if(um==0){
    UL=ul;UT=ut;GPROG=G;
    UM=ratio*ul*ul*ul/(ut*ut);
  }
  else if(G==0){
    UL=ul;UT=ut;UM=um;
    GPROG=GCONST/(ul*ul*ul/(um*ut*ut));
  }

  //CALCULATE DERIVED CONSTANTS
  //Radiation pressure: beta = alpha Qpr / (rho R) (Pater & Lissauer)
  ALPHA=3*LSUN/(4*PI*CLIGHT*GCONST*MSUN)/4;
  //fprintf(stdout,"SI: %e\n",ALPHA);
  ALPHA/=(UM/(UL*UL));
  //fprintf(stdout,"Program Units %e\n",ALPHA);
  //exit(0);
  
  return 0;
}

double tOrbit(double e[])
{
  double t;
  double a;
  if(e[ECC]>=1) a=e[RP];
  else
    a=e[RP]/(1-e[ECC]);

  t=sqrt(4*PI*PI*a*a*a/e[MU]);
  
  return t;
}

void rotateState(double R[3][3],double x[],double xr[])
{
  double *r,*v,rr[3],vr[3];
  r=x;
  v=x+3;
  mxv_c(R,r,rr);
  mxv_c(R,v,vr);
  memcpy(xr,rr,3*sizeof(double));
  memcpy(xr+3,vr,3*sizeof(double));
}

void copyState(double dest[],double origin[])
{
  memcpy(dest,origin,6*sizeof(double));
}

void stateAdd(double x1[],double x2[],double s[])
{
  vadd_c(x1,x2,s);
  vadd_c(x1+3,x2+3,s+3);
}

void stateOrbit(double t,double e[],double x[])
{
  e[T0]=(e[TPER]-t);
  conics_c(e,0,x);
}

void orbitNormal(double e[],double n[])
{
  double x1[6],x2[6];
  double P=tOrbit(e);
  stateOrbit(e[TPER],e,x1);
  stateOrbit(e[TPER]+P/10.0,e,x2);
  vcrss_c(x1,x2,n);
  unorm_c(n,n,&P);
}

void stateUnit(double x[],double ux[])
{
  double norm;
  unorm_c(x,ux,&norm);
  unorm_c(x+3,ux+3,&norm);
}

double* realAlloc(int n)
{
  double* v;
  v=(double*)calloc(n,sizeof(double));
  return v;
}

double** stateAlloc(int n)
{
  int i;
  double** x;
  x=(double**)calloc(n,sizeof(double*));
  for(i=0;i<n;i++)
    x[i]=(double*)calloc(6,sizeof(double));
  return x;
}

double volumeElement(double r1,double r2,double t1,double t2,double df)
{
  return 1./3*(r2*r2*r2-r1*r1*r1)*(cos(t1)-cos(t2))*df;
}

void divisionSphere(int nfrag,double* Rf)
{
  double R=1;
  double V;
  int it,ir,ip;
  int nv;
  double r,dr,t,dt,f,df,v;
  double factor;

  nv=(int)pow(nfrag,1./3);
  dr=R/nv;
  dt=PI/nv;
  df=2*PI/nv;
  V=4*PI/3*R*R*R;

  printf("Divisions per dimension: %d\n",nv);
  printf("Total volume: %e\n",V);

  v=volumeElement(R-dr,R,0,dt,df);
  printf("Volume fragment: %e (%e/%d = %e)\n",v,V,nfrag,V/nfrag);

  it=1;
  t=dt;
  printf("Initial dt = %e\n",dt);
  while(t<=PI/2){
    factor=sin(t)/(sin(t+dt));
    printf("Factor at t = %e = %e\n",t,factor);
    dt*=factor;
    v=volumeElement(R-dr,R,t,t+dt,df);
    printf("\tdt = %e\n",dt);
    printf("\tVolume = %e\n",v);
    t+=dt;
    it++;
  }
  printf("Final value of t: %e (x pi/2)\n",t/(PI/2));
  printf("Number of slices: %d\n",it);
  exit(0);

  t=0;
  dr*=1.0;
  printf("Initial r space: %e\n",dr);
  r=R-dr;
  ir=1;
  while(r>dr){
    dr*=(r*r)/((r-dr)*(r-dr));
    printf("r space at r = %e: %e\n",r,dr);
    r-=dr;
    ir++;
  }
  printf("Final value of r: %e\n",r);
  printf("Number of slices: %d\n",ir);

  exit(0);
}

gsl_rng *RANGEN;
void randomInit()
{
  RANGEN=gsl_rng_alloc(gsl_rng_taus2);

  #ifdef TIMESEED
  gsl_rng_set(RANGEN,time(NULL));
  #else
  gsl_rng_set(RANGEN,RANSEED);
  #endif
}

double randReal(void)
{
  return gsl_rng_uniform(RANGEN);
}

int gravSystem(double t,const double y[],double dydt[],void* params)
{
  int i,j,k;
  double *ps=(double*) params;
  double *Rs=ps+1;
  int nsys=(int)ps[0];
  int nfrag=nsys/6-1;

  //////////////////////////////////////////
  //ZEROED dydt
  //////////////////////////////////////////
  memset(dydt,0,nsys*sizeof(dydt[0]));
  
  //////////////////////////////////////////
  //drdt
  //////////////////////////////////////////
  for(i=0;i<=nfrag;i++){
    k=6*i;
    memcpy(dydt+k,y+k+3,3*sizeof(double));
  }

  //////////////////////////////////////////
  //dvdt
  //////////////////////////////////////////
  double beta,D,D3,M,R;
  double gs[3],rad[3],evap[3],gb[3];
  for(i=0;i<=nfrag;i++){
    k=6*i;
    D=vnorm_c(y+k);D3=D*D*D;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //PROPERTIES OF THE PARTICLE
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    beta=ALPHA*QPR/(RHODUST*Rs[i]);
    fprintf(stdout,"i = %d:\n",i);
    fprintf(stdout,"\tRHODUST = %e UM/UL^3 = %e kg/m^3\n",
	    RHODUST,RHODUST*UM/(UL*UL*UL));
    fprintf(stdout,"\tRs = %e UL = %e m\n",
	    Rs[i],Rs[i]*UL);
    fprintf(stdout,"\tbeta = %e\n",
	    beta);
    exit(0);
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //GRAVITATIONAL FORCE FROM THE SUN
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gs[0]=-GPROG*y[0+k]/D3;
    gs[1]=-GPROG*y[1+k]/D3;
    gs[2]=-GPROG*y[2+k]/D3;
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //RADIATION FORCES
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //See de Pater & Lissauer, Sec. 2.7
    //Radiation Pressure
    rad[0]=-beta*gs[0];
    rad[1]=-beta*gs[1];
    rad[2]=-beta*gs[2];

    //Poynting-Robertson

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //EVAPORATION RECOIL
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    evap[0]=0.0;
    evap[1]=0.0;
    evap[2]=0.0;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //GRAVITATIONAL FORCE FROM OTHER BODIES
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gb[0]=0.0;
    gb[1]=0.0;
    gb[2]=0.0;
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //TOTAL
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dydt[3+k]=gs[0]+rad[0]+gb[0]+evap[0];
    dydt[4+k]=gs[1]+rad[1]+gb[1]+evap[1];
    dydt[5+k]=gs[2]+rad[2]+gb[2]+evap[2];
  }

  return 0;
}

void earthObservations(double t,double eE[],double eC[],double y[],double x[])
{
  int i,k;
  int nfrag;
  double RE[3][3];
  double xE[6],xC[6],xCE[6];
  double xobs[3],yobs[3],zobs[3];
  double zaxis[]={0,0,1};

  //CALCULATE EARTH POSITION
  stateOrbit(t,eE,xE);

  //CALCULATE COMET CM POSITION
  stateOrbit(t,eC,xC);

  //ROTATION MATRIX
  vsubg_c(xC,xE,6,xCE);
  vpack_c(-xCE[0],-xCE[1],-xCE[2],zobs);
  vcrss_c(zaxis,zobs,xobs);
  twovec_c(xobs,XAXIS,zobs,ZAXIS,RE);
  
  //ROTATE POSITIONS
  vsubg_c(y,xE,6,xCE);
  rotateState(RE,xCE,x);
}

void vxk_c(double v[],double k,double vout[])
{
  vout[0]=k*v[0];
  vout[1]=k*v[1];
  vout[2]=k*v[2];
}
