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
#include <gsl/gsl_roots.h>

#ifndef CLANG
using namespace std;
#endif

//////////////////////////////////////////////////////////////////////////////////
// GLOBAL BEHAVIOR
//////////////////////////////////////////////////////////////////////////////////
#define ALLOW_COLISION /*ALLOW COLISION*/
//#define SHOW_COLISION
//#define TIMESEED
#define RANSEED 1


//////////////////////////////////////////////////////////////////////////////////
//SPICE INFO, MACROS & KERNELS
//////////////////////////////////////////////////////////////////////////////////
//PLANETS
/*
  1 MERCURY BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
  2 VENUS BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
  3 EARTH BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
  4 MARS BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
  5 JUPITER BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
  6 SATURN BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
  7 URANUS BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
  8 NEPTUNE BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
  9 PLUTO BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
  10 SUN w.r.t. 0 SOLAR SYSTEM BARYCENTER
  199 MERCURY w.r.t. 1 MERCURY BARYCENTER
  299 VENUS w.r.t. 2 VENUS BARYCENTER
  301 MOON w.r.t. 3 EARTH BARYCENTER
  399 EARTH w.r.t. 3 EARTH BARYCENTER
  499 MARS w.r.t. 4 MARS BARYCENTER
*/
#define ABCORR "NONE"
#define FRAME  "ECLIPJ2000"
#define ABSKERNEL HOME
#define DIRKERNEL "util/kernels"

//GENERIC KERNELS
//OBTAIN KERNELS AND INFO HERE: http://naif.jpl.nasa.gov/naif/data.html
#define PLSPK "de421.bsp" /* PLANETARY POSITIONS */
#define PLDATA "pck00010.tpc" /* PLANETARY DATA */
#define PLMASS "de-403-masses.tpc" /* PLANETARY MASSES */
#define LEAP "naif0010.tls" /* LEAP SECONDS */

//SPECIFIC KERNELS
//GENERATE NEW KERNELS HERE: http://ssd.jpl.nasa.gov/x/spk.html
/*
Object Name        = ISON (C/2012 S1)
Primary SPKID      = 1003203
Primary designation= C/2012 S1
Aliases            = K12S010
Kernel generated with : http://ssd.jpl.nasa.gov/x/spk.html
 */
#define ISONID "1003203" /*SPKID */
#define ISONSPK "ison-perihelion.bsp" /* ISON */
#define ISONTP "11/28/2013 18:51:11.52" /*PERIHELION*/
#define ISONFT "11/10/2013 00:00:00" /*FRAGMENTATION*/

/*
Object Name        = ENSOR (C/1925 X1)
Primary SPKID      = 9905348
Primary designation= C/1925 X1
Aliases            = 
Kernel generated with : http://ssd.jpl.nasa.gov/x/spk.html
 */
#define ENSORID "9905348" /*SPKID */
#define ENSORSPK "ensor.bsp" /*ENSOR (C/1925 X1)*/
#define ENSORTP "02/11/1926 23:20:49.92" /*PERIHELION*/
#define ENSORFT "02/11/1926 23:20:49.92" /*FRAGMENTATION*/

/*
Object Name        = LINEAR (C/1999 S4)
Primary SPKID      = 1000273
Primary designation= C/1999 S4
Aliases            = J99S040 
Kernel generated with : http://ssd.jpl.nasa.gov/x/spk.html
Desintegration report: 
http://www.ing.iac.es/PR/press/ing300.html
http://apod.nasa.gov/apod/ap000811.html
 */
#define LINEARID "1000273" /*SPKID*/
#define LINEARSPK "linear.bsp" /*LINEAR (C/1999 S4)*/
#define LINEARTP "07/26/2000 03:55:35.0" /*PERIHELION*/
#define LINEARFT "07/25/2000 00:00:00" /*FRAGMENTATION*/

/*
Object Name        = Elenin (C/2010 X1)
Primary SPKID      = 1003113
Primary designation= C/2010 X1
Aliases            = K10X010
Kernel generated with : http://ssd.jpl.nasa.gov/x/spk.html
Desintegration report: http://www.universetoday.com/88494,
http://www.thunderbolts.info/wp/2011/10/06/comet-elenin%E2%80%94the-debate-that-never-happened/
 */
#define ELENINID "1003113" /*SPKID*/
#define ELENINSPK "elenin.bsp" /*ELENIN (C/2010 X1)*/
#define ELENINTP "09/10/2011 17:15:51.8" /*PERIHELION*/
#define ELENINFT "08/19/2011 00:00:00" /*FRAGMENTATION*/

/*
Object Name        = Schwassmann-Wachmann 3
Primary SPKID      = 1000394
Primary designation= 73P
Aliases            = 1979 P1, 1930 J1, 1994w, 1990 VIII, 1989d1, 1979g, 1979 VIII, 1930d, 1930 VI, 4000073
Kernel generated with : http://ssd.jpl.nasa.gov/x/spk.html
Desintegration report: 
 */
#define SW3ID "1000394" /*SPKID*/
#define SW3SPK "sw3.bsp" /*73P/SCHWASSMANN-WACHMANN 3*/
#define SW3TP "09/22/1995 21:21:23.9" /*PERIHELION*/
#define SW3FT "09/22/1995 00:00:00" /*FRAGMENTATION*/

/*
Object Name        = Schwassmann-Wachmann 3, Fragment B
Primary SPKID      = 1000320
Primary designation= 73P-B
Aliases            = Kernel generated with : http://ssd.jpl.nasa.gov/x/spk.html
Desintegration report: 
http://arxiv.org/abs/0904.4733
 */
#define SW3BID "1000320" /*SPKID*/
#define SW3BSPK "sw3B.bsp" /*73P-B*/
#define SW3BTP "09/22/1995 21:21:23.9" /*PERIHELION*/
#define SW3BFT "04/01/2006 00:00:00" /*FRAGMENTATION*/

//COMETARY KERNELS
//GET OTHER COMETARY KERNELS HERE: http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/comets/
//FOR A LIST OF SPKIDS: util/kernels/pcomets_dastcom_v1.txt
#define COMETSPK "pcomets_v1_ieee.bsp" /* COMETS */

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
#define DAY (GSL_CONST_MKSA_DAY) //s
#define MEARTH (5.9736E24) //kg
#define REARTH (6.371E6) //m
#define GCONST (GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT) // m^3 / (kg s^2)
#define CLIGHT (GSL_CONST_MKSA_SPEED_OF_LIGHT) // m / s
#define KBOLTZ (GSL_CONST_MKSA_BOLTZMANN) //j / K
#define MH2O 2.98897E-26 //kg

#define QPR 1.0 //Radiation pressure coefficient
#define ZO 3E17/(1E-4)  //Molecules / m^2 s

#define PI M_PI
#define PI2 (PI*PI)
#define D2R (PI/180)

//////////////////////////////////////////////////////////////////////////////////
//GLOBAL VARIABLES
//////////////////////////////////////////////////////////////////////////////////
typedef ConstSpiceChar* date;
typedef double real;
typedef void* params;
typedef FILE* file;

real UL,UM,UT,GPROG;
real ALPHA;
real RHODUST,RHOCOMET,RHOVOLATILES,FR;
int NPARTICLES,NLARGE,NDEBRIS;
double AR,AT,AN,AMREF;

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
  11 - AXISXY: Azimutal angle of body at periapse (x-axis: sun-periapse)
  12 - MASS: Mass
 */
#define NELEMS 13

enum ElementsEnum {RP,ECC,INC,LNODE,ARGP,M0,T0,MU,
		   TPER,PER,AXISZ,AXISXY,EMASS,ENDELEM};
char Elements[][10]={"RP","ECC","INC","LNODE","ARGP","M0","T0","MU",
		     "TPER","PER","AXISZ","AXISXY","EMASS"};

//STATE VECTOR
#define NSTATE 7

enum StateEnum {XPOS,YPOS,ZPOS,XVEL,YVEL,ZVEL,MASS,ENDSTATE};
char State[][10]={"XPOS","YPOS","ZPOS","XVEL","YVEL","ZVEL","MASS"};

//ROTATION MATRIX
enum Axis {_AXIS_,XAXIS,YAXIS,ZAXIS};

//TYPE OF FRAGMENT
enum Fragment {_FRAGMENT_,LARGE,BOULDER,DUST};

//////////////////////////////////////////////////////////////////////////////////
//INITIALIZE
//////////////////////////////////////////////////////////////////////////////////
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
  ALPHA/=(UM/(UL*UL));
  
  return 0;
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

char* HOME;
#define MAXSTR 500
char* kernelLocation(const char* kernel)
{
  char *kerneloc;
  kerneloc=(char*)malloc(MAXSTR*sizeof(char));
  if(HOME==NULL)
    HOME=getenv("PWD");
  strcat(kerneloc,ABSKERNEL);
  strcat(kerneloc,"/");
  strcat(kerneloc,DIRKERNEL);
  strcat(kerneloc,"/");
  strcat(kerneloc,kernel);
  return kerneloc;
}

void cometsynInit(double ul,double um,double ut,double uG)
{
  //INITIALIZE UNITS
  setUnits(ul,um,ut,uG);

  //INITIALIZE RANDOM NUMBERS GENERATORS
  randomInit();
  
  //LOAD GENERIC SPICE KERNELS
  furnsh_c(kernelLocation(PLDATA));
  furnsh_c(kernelLocation(PLMASS));
  furnsh_c(kernelLocation(PLSPK));
  furnsh_c(kernelLocation(LEAP));

  //LOAD COMETARY KERNELS
  furnsh_c(kernelLocation(ISONSPK));
  furnsh_c(kernelLocation(COMETSPK));
}

//////////////////////////////////////////////////////////////////////////////////
//ROUTINES
//////////////////////////////////////////////////////////////////////////////////
void vxk_c(double v[],double k,double vout[])
{
  vout[0]=k*v[0];
  vout[1]=k*v[1];
  vout[2]=k*v[2];
}

void fprintf_vec(FILE* stream,const char* fmt,const double x[],int end,bool qlabel=true,int ini=0,bool endline=true)
{
  for(int i=ini;i<end;i++){
    if(qlabel) fprintf(stream,"%d:",i);
    fprintf(stream,fmt,x[i]);
  }
  if(endline)
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
  for(int i=0;i<NSTATE;i++){
    if(qname) fprintf(stream,"%s:",State[i]);
    fprintf(stream,fmt,x[i]);
  }
  fprintf(stream,"\n");
}

double orbitTimeScale(double q,double e,double mu)
{
  double t;
  double a;
  if(e>=1) a=10*q;
  else a=q/(1-e);
  t=sqrt(4*PI*PI*a*a*a/mu);
  return t;
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
  memcpy(dest,origin,NSTATE*sizeof(double));
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
  x[MASS]=e[EMASS];
}

void orbitNormal2(double e[],double n[])
{
  double x1[NSTATE],x2[NSTATE];
  double P=tOrbit(e);
  stateOrbit(e[TPER],e,x1);
  stateOrbit(e[TPER]+P/10.0,e,x2);
  vcrss_c(x1,x2,n);
  unorm_c(n,n,&P);
}

void setunitsState(double state[NSTATE])
{
  vxk_c(state,1E3/UL,state);
  vxk_c(state+3,1E3/(UL/UT),state+3);
  state[6]/=UM;
}

void stateObject(const char* spkid,double t,double x[])
{
  double xsun[NSTATE],light;
  spkezr_c("SUN",t,"ECLIPJ2000","NONE","SOLAR SYSTEM BARYCENTER",xsun,&light);
  spkezr_c(spkid,t,"ECLIPJ2000","NONE","SOLAR SYSTEM BARYCENTER",x,&light);
  vsubg_c(x,xsun,NSTATE-1,x);
  setunitsState(x);
}

void orbitDirector(const char* spkid,double t,double r[],double d[],double n[])
{
  double x1[NSTATE],x2[NSTATE],elements[NELEMS],P;

  stateObject(spkid,t,x1);
  oscelt_c(x1,t,GPROG*MSUN/UM,elements);
  P=orbitTimeScale(elements[RP],elements[ECC],elements[MU]);
  stateObject(spkid,t+P/10,x2);

  memcpy(r,x1,3*sizeof(x1[0]));
  vhat_c(r,r);

  vcrss_c(x1,x2,n);
  vhat_c(n,n);

  vcrss_c(r,n,d);
  vhat_c(d,d);
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
    x[i]=(double*)calloc(NSTATE,sizeof(double));
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

double randReal(void)
{
  return gsl_rng_uniform(RANGEN);
}

/*
  Ishiguro et al. (2006) originally from Marsden & Zekanina (1973)

  Z = Zo g(r)

  Where Zo = 3 x 10^17 molecules / cm^2 s 
*/
double sublimationRate(double r /*m*/,
		       /*WATER SUBLIMATION*/
		       double ao=0.111262/*m/s^2*/,
		       double ro=2.808/*AU*/,
		       double mo=2.15,
		       double no=5.093,
		       double ko=4.6142
		       )
{
  double Zo,Z;
  double g;

  Zo=ZO;
  g=ao/pow((r/AU/ro),mo)/pow((1+pow((r/AU/ro),no)),ko);
  Z=Zo*g;

  return Z;
}

/*
  Drying function
 */
double dryingFunction(double r)
{
  double dF=1;
  double Z;

  #ifdef DRYING_CONSTANT
  //==================================================
  //CONSTANT SUBLIMATION RATE AFTER 1 AU
  //==================================================
  if(r<AU){
    Z=sublimationRate(r);
    dF=ZO/Z;
  }
  #endif

  #ifdef DRYING_EXPONENTIAL
  //==================================================
  //DECAYING SUBLIMATION
  //==================================================
  if(r<AU){
    Z=sublimationRate(r);
    dF=ZO/Z*pow(r/AU,1.58);
  }
  #endif

  #ifdef DRYING_NULL
  //==================================================
  //NO DYRING
  //==================================================
  dF=1.0;
  #endif

  return dF;
}

/*
  See Sosa & Fernandez (2008)
 */
double outflowVelocity(double r,
		       /*WATER*/
		       double To=323.0,
		       double nT=0.5,
		       double csi=0.5)
{
  double T;
  double vth,u;

  T=To/pow(r/AU,nT);

  vth=sqrt(8*KBOLTZ*T/(PI*MH2O));

  /*
  fprintf(stdout,"\tTemperature: %e K / s\n",T);
  fprintf(stdout,"\tThermal velocity: %e m/s / s\n",vth);
  //*/

  u=csi*vth;
  return u;
}

void rocketAcceleration(double r,
			  double M,
			  double R,
			  /*Drying function*/
			  double (*dFunc)(double),
			  double *J,double *dMdt)
{
  double Z,A,Q,u,dF;

  A=2*PI*R*R;/*m^2*/

  Z=sublimationRate(r);/*molecules/m^2 s*/

  Q=A*Z;/*molecules/s*/

  dF=dFunc(r);/*Adimensional: drying function*/
  Q*=dF;

  u=outflowVelocity(r);/*m/s*/

  *dMdt=Q*MH2O; /*kg*/

  *J=Q*MH2O*u/M;/*m/s^2*/

  /*
  fprintf(stdout,"\tDistance = %e AU\n",r/AU);
  fprintf(stdout,"\tArea for radius R = %e m: %e m^2\n",R,A);
  fprintf(stdout,"\tSublimation rate: %e molec / m^2 s\n",Z);
  fprintf(stdout,"\tSublimation flux Q: %e molec\n",Q);
  fprintf(stdout,"\tDrying function: %e\n",dF);
  fprintf(stdout,"\tOutflow velocity: %e m/s\n",u);
  fprintf(stdout,"\tTotal outflow momentum: %e m/s\n",Q*MH2O*u);
  fprintf(stdout,"\tComet mass : %e kg\n",M);
  fprintf(stdout,"\tMass loss: %e kg/s\n",*dMdt);
  fprintf(stdout,"\tAcceleration : %e m/s^2\n",*J);
  //*/

}

int gravSystem(double t,const double y[],double dydt[],void* params)
{
  int i,j,k,kj;

  //INPUT PARAMETERS
  double *ps=(double*) params;
  int nsys=(int)ps[0];
  int nlarge=(int)ps[1];
  int nfrag=nsys/NSTATE-1;
  int ndebris=nfrag-nlarge;

  //RADIUS AND MASSES OF FRAGMENTS
  double *Rs=ps+2;
  double *Ms=ps+2+nfrag;

  double Mp,Rf,Rp,Mc,dMdt,Ar;
  double frp,rhom,Md;
  const double *r,*v;
  const double *ri,*rj;
  double *dyjdt,*dyidt;
  double Rij[3];
  double beta,D,Dsoft,epsoft,D3,M,R;
  double gs[]={0,0,0},rad[]={0,0,0},rocket[]={0,0,0},gb[]={0,0,0};
  double ur[3],ut[3],vr[3],vt[3],pf[3],vrmag=0,vtmag=0,fr=1,ft=0,gmag=0;
  bool qcolision;
  double Mf;

  //////////////////////////////////////////
  //INITIALIZE: dx/dt = 0
  //////////////////////////////////////////
  memset(dydt,0,nsys*sizeof(dydt[0]));

  //////////////////////////////////////////
  //dr_i/dt = v_i
  //////////////////////////////////////////
  for(i=0;i<=nfrag;i++){
    k=NSTATE*i;
    memcpy(dydt+k,y+k+3,3*sizeof(double));
  }

  //////////////////////////////////////////
  //dv_i/dt
  //////////////////////////////////////////
 init:
  for(i=0;i<=nfrag;i++){

    if(i==0) Mf=ps[2*nfrag+3];
    else Mf=Ms[i-1];
    //fprintf(stdout,"\nIntegrating particle: %d (Mf = %.11e, Mp = %.11e)\n",i,Mf,y[NSTATE*i+6]);

    //AVOID INTEGRATION OF PARTICLES ALREADY ABSORBED
    if(i>0&&Ms[i-1]==0){
      continue;
    }
    
    //FRESH RADIUS
    if(i>0) Rf=Rs[i-1];
    else Rf=ps[2*nfrag+2];

    qcolision=false;
    k=NSTATE*i;

    //ACTUAL RADIUS
    Md=FR*Mf;
    Mp=y[6+k];
    frp=Md/Mp;
    rhom=1/((1-frp)/RHOVOLATILES+frp/RHODUST);
    Rp=pow((Mp/(4*PI/3*rhom)),1./3);
    /*
    fprintf(stdout,"\tMd = %.11e, Mp = %.11e, frp = %.11e, rhom = %.11e, Rp = %.11e, Rf = %.11e\n",
	    Md,Mp,frp,rhom*UM/(UL*UL*UL),Rp,Rf);
    //*/

    //DISTANCE TO SUN
    D=vnorm_c(y+k);D3=D*D*D;

    //COLISION INTO THE SUN
    if(D<RSUN/UL){
      #ifdef SHOW_COLISION
      fprintf(stdout,"\t\tParticle %d has collided with the Sun\n",i-1);
      #endif
      //RESET GRADIENT FOR COLLIDED PARTICLE
      memset(dydt+k,0,NSTATE*sizeof(dydt[0]));
      //SET TO ZERO MASS AND RADIUS OF IMPACTOR
      Rs[i-1]=Ms[i-1]=0;
      //REDUCE NUMBER OF PARTICLES
      NPARTICLES--;
      if(i<=nlarge) NLARGE--;
      else NDEBRIS--;
      //RESET CALCULATION
      goto init;
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //GRAVITATIONAL FORCE FROM THE SUN
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gs[0]=0.0;
    gs[1]=0.0;
    gs[2]=0.0;
    #ifdef SUN_FORCE
    gs[0]=-GPROG*y[0+k]/D3;
    gs[1]=-GPROG*y[1+k]/D3;
    gs[2]=-GPROG*y[2+k]/D3;
    #else
    memset(gs,0,3*sizeof(gs[0]));
    #endif

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //ROCKET FORCE: FORCE COMING FROM SUBLIMATION
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #ifdef ROCKET_FORCE
    double go,g,J;

    //ONLY FOR LARGE PARTICLES
    Mp=y[6+k];
    if(i>0) Mf=Ms[i-1];
    else Mf=ps[2*nfrag+3];
    
    if(Mp<=FR*Mf && FR<1){
      /*
      fprintf(stdout,"\t\tPARTICLE %d HAS DRY OUT (Mp = %e < Mr = %e [=%e x %e])\n",
	      i,Mp,FR*Mf,FR,Mf);
      */
      Ms[i-1]=FR*Mf;
    }
    if(Rp>RMIN_ROCKET/UL /*Minimum size of evaporable rubble*/ 
       && 
       Mp>FR*Mf /*Rubble still has enough ice to evaporate*/
       ){
      rocketAcceleration(D*UL,Mp*UM,Rp*UL,dryingFunction,&J,&dMdt);
      g=sublimationRate(D*UL)/(3E17/1E-4);
      Ar=J/(AU/(DAY*DAY))/g;

      J/=(UL/(UT*UT));
      dMdt/=(UM/UT);

      rocket[0]=J*y[0+k]/D3;
      rocket[1]=J*y[1+k]/D3;
      rocket[2]=J*y[2+k]/D3;

      /*
      fprintf(stdout,"****************************************\n");
      fprintf(stdout,"ROCKET ACCELERATION: PARTICLE %d\n",i);
      fprintf(stdout,"Present fragment mass: %.11e UM = %.11e kg\n",Mp,Mp*UM);
      fprintf(stdout,"Fresh fragment mass: %.11e UM = %.11e kg\n",Mf,Mf*UM);
      fprintf(stdout,"Rocky fragment mass: %.11e UM = %.11e kg\n",FR*Mf,FR*Mf*UM);
      fprintf(stdout,"Rocket acceleration: %e UL/UT^2\n",J);
      fprintf(stdout,"Equivalent Ar: %e AU/day^2\n",Ar);
      fprintf(stdout,"Mass loss: %e UM/UT\n",dMdt);
      fprintf(stdout,"Rocket force (Sosa & Fernandez):");
      fprintf_vec(stdout,"%e ",rocket,3);
      fprintf(stdout,"****************************************\n");
      //*/
    }else{ 
      /*
      fprintf(stdout,"****************************************\n");
      fprintf(stdout,"ROCKET ACCELERATION: PARTICLE %d\n",i);
      fprintf(stdout,"No rocket effect\n");
      fprintf(stdout,"Rp = %e, Rpmin = %e\n",Rp,RMIN_ROCKET/UL);
      fprintf(stdout,"Mp = %e, FR Mf = %e\n",Mp,FR*Mf);
      fprintf(stdout,"****************************************\n");
      //*/
      dMdt=0;
      memset(rocket,0,3*sizeof(rocket[0]));
    }      
    #else
    dMdt=0;
    memset(rocket,0,3*sizeof(rocket[0]));
    #endif
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //GRAVITATIONAL FORCE FROM OTHER BODIES
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //fprintf(stdout,"Particle %d testing colisions\n",i);
    
    #ifdef FRAG_FORCE
    gb[0]=0.0;
    gb[1]=0.0;
    gb[2]=0.0;
    for(j=1;j<=nlarge&&i>0;j++){
      if(i==j || Rs[j-1]==0) continue;
      kj=NSTATE*j;
      ri=y+k;
      rj=y+kj;
      Mp=Ms[j-1];
      epsoft=0.0;

      vsub_c(ri,rj,Rij);
      D=vnorm_c(Rij);

      #ifdef ALLOW_COLISION
      //COLISION!
      if(D<(Rs[j-1]+Rs[i-1])){

	//fprintf(stdout,"Particle %d colision: %e < %e + %e\n",i,D,Rs[j-1],Rs[i-1]);

	//ANNOUNCE IT
        #ifdef SHOW_COLISION
	fprintf(stdout,"\t\tParticle %d has collided with particle %d\n",i-1,j-1);
	#endif
	/*
	fprintf(stdout,"\t\t\tAbsorbed particle: %d, Ri = %e m, Mi = %e kg\n",i,Rs[i-1]*UL,Ms[i-1]*UM);
	fprintf(stdout,"\t\t\tAbsorber fragment: %d, Rj = %e m, Mj = %e kg\n",j,Rs[j-1]*UL,Ms[j-1]*UM);
	fprintf(stdout,"\t\t\tDistance: %e m < %e m\n",D*UL,(Rs[i-1]+Rs[j-1])*UL);
	*/

	//RESET GRADIENT FOR COLLIDED PARTICLE
	memset(dydt+k,0,NSTATE*sizeof(dydt[0]));

	//CHANGE VELOCITY OF HIT PARTICLE
	dyidt=dydt+k;
	dyjdt=dydt+kj;

	/*
	if((i-1)<nlarge&&(j-1)<nlarge){

	  fprintf(stdout,"Velocity of impactor(%d):",k);fprintf_vec(stdout,"%e ",ri+3,3);
	  fprintf(stdout,"Velocity of hit before(%d):",kj);fprintf_vec(stdout,"%e ",rj+3,3);

	  fprintf(stdout,"Gradient of impactor:");fprintf_vec(stdout,"%e ",dyidt,NSTATE);
	  fprintf(stdout,"Gradient of hit before:");fprintf_vec(stdout,"%e ",dyjdt,NSTATE);
	}
	//*/

	dyjdt[0]=(Ms[j-1]*rj[3]+Ms[i-1]*ri[3])/(Ms[j-1]+Ms[i-1]);
	dyjdt[1]=(Ms[j-1]*rj[4]+Ms[i-1]*ri[4])/(Ms[j-1]+Ms[i-1]);
	dyjdt[2]=(Ms[j-1]*rj[5]+Ms[i-1]*ri[5])/(Ms[j-1]+Ms[i-1]);
	
	//INCREASE MASS AND RADIUS CORRESPONDINGLY
	Ms[j-1]+=Ms[i-1];
	Rs[j-1]=pow(Ms[j-1]/(4*PI/3*RHODUST),1./3);

	//SET TO ZERO MASS AND RADIUS OF IMPACTOR
	Rs[i-1]=Ms[i-1]=0;

	/*
	if((i-1)<nlarge&&(j-1)<nlarge){
	  fprintf(stdout,"Velocity of hit after:");fprintf_vec(stdout,"%e ",dyjdt,NSTATE);
	  fprintf(stdout,"Radius after of particle hit: %e m\n",Rs[j-1]*UL);
	  exit(0);
	}
	//*/

	//REDUCE NUMBER OF PARTICLES
	NPARTICLES--;
	if(i<=nlarge) NLARGE--;
	else NDEBRIS--;

	//RESET CALCULATION
	goto init;
      }
      #endif
      Dsoft=sqrt(D*D+epsoft*epsoft);
      D3=Dsoft*Dsoft*Dsoft;
      gb[0]+=-GPROG*Mp*Rij[0]/D3;
      gb[1]+=-GPROG*Mp*Rij[1]/D3;
      gb[2]+=-GPROG*Mp*Rij[2]/D3;

      /*
      fprintf(stdout,"Force on %d from %d...\n",i,j);
      fprintf(stdout,"Rs:");fprintf(stdout,"%e\n",Rs[j-1]);
      fprintf(stdout,"Ms:");fprintf(stdout,"%e\n",Mp);
      fprintf(stdout,"ri:");fprintf_vec(stdout,"%-+23.17e ",ri,3);
      fprintf(stdout,"rj:");fprintf_vec(stdout,"%-+23.17e ",rj,3);
      fprintf(stdout,"D:");fprintf(stdout,"%e km\n",D*UL);
      fprintf(stdout,"Dsoft:");fprintf(stdout,"%e\n",Dsoft);
      fprintf(stdout,"D3:");fprintf(stdout,"%e\n",D3);
      fprintf(stdout,"Force:");fprintf_vec(stdout,"%e ",gb,3);
      //*/
    }
    #else
    memset(gb,0,3*sizeof(gb[0]));
    #endif

    //fprintf(stdout,"Particle %d pass colisions\n",i);
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //RADIATION FORCES
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //See de Pater & Lissauer, Sec. 2.7

    #ifdef RADIATION_FORCE
    if(i>0)
      beta=ALPHA*QPR/(RHODUST*Rs[i-1]);
    else
      beta=0;
    #ifdef PR_CORRECTION
    //&&&&&&&&&&&&&&&&&&&&
    //POYNTING-ROBERTSON
    //&&&&&&&&&&&&&&&&&&&&
    //Eqs. 2.49a,b, Pater-Lissauer
    r=y+k;//Radius
    v=y+k+3;//Velocity

    //Radial velocity component
    vhat_c(r,ur);
    vrmag=vdot_c(v,ur);
    vxk_c(ur,vrmag,vr);

    //Tangential velocity component
    vsub_c(v,vr,vt);
    vhat_c(vt,ut);
    vtmag=vnorm_c(vt);

    //Radial correction
    fr=(1-2*(fabs(vrmag)*UL/UT)/CLIGHT);
    //Tangential correction
    ft=(vtmag*UL/UT)/CLIGHT;
    
    //PR-correction
    gmag=vnorm_c(gs);
    rad[0]=beta*gmag*(fr*ur[0]-ft*ut[0]);
    rad[1]=beta*gmag*(fr*ur[1]-ft*ut[1]);
    rad[2]=beta*gmag*(fr*ur[2]-ft*ut[2]);

    /*
    if(i>nlarge){
      fprintf(stdout,"r=");fprintf_vec(stdout,"%e ",r,3);
      fprintf(stdout,"v=");fprintf_vec(stdout,"%e ",v,3);
      fprintf(stdout,"ur=");fprintf_vec(stdout,"%e ",ur,3);
      fprintf(stdout,"vrmag=");fprintf(stdout,"%e\n",vrmag);
      fprintf(stdout,"vr=");fprintf_vec(stdout,"%e ",vr,3);
      fprintf(stdout,"vt=");fprintf_vec(stdout,"%e ",vt,3);
      fprintf(stdout,"vtmag=");fprintf(stdout,"%e\n",vtmag);
      fprintf(stdout,"ut=");fprintf_vec(stdout,"%e ",ut,3);
      fprintf(stdout,"fr=");fprintf(stdout,"%e\n",fr);
      fprintf(stdout,"ft=");fprintf(stdout,"%e\n",ft);
      fprintf(stdout,"rad=");fprintf_vec(stdout,"%e ",rad,3);
      exit(0);
    }
    //*/
    #else /*NO POYNTING-ROBERTSON CORRECTION*/
    //&&&&&&&&&&&&&&&&&&&&
    //RADIATION PRESSURE
    //&&&&&&&&&&&&&&&&&&&&
    rad[0]=-beta*gs[0];
    rad[1]=-beta*gs[1];
    rad[2]=-beta*gs[2];
    #endif
    #else
    memset(rad,0,3*sizeof(rad[0]));
    #endif
    /*
    if(i>nlarge){
      fprintf(stdout,"rad=");fprintf_vec(stdout,"%e ",rad,3);
      exit(0);
    }
    //*/

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //TOTAL
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dydt[3+k]=gs[0]+rad[0]+gb[0]+rocket[0];
    dydt[4+k]=gs[1]+rad[1]+gb[1]+rocket[1];
    dydt[5+k]=gs[2]+rad[2]+gb[2]+rocket[2];

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //MASS-LOSS
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dydt[6+k]=-dMdt;

    /*
    if(i>nlarge || true){
      if(i==0) Mf=ps[2*nfrag+3];
      else Mf=Ms[i-1];
      fprintf(stdout,"Particle: %d\n",i);
      fprintf(stdout,"Fresh mass: %.11e\n",Mf);
      fprintf(stdout,"Present mass: %.11e\n",y[6+k]);
      fprintf(stdout,"Sun gravitatinal force:");
      fprintf_vec(stdout,"%e ",gs,3);
      fprintf(stdout,"Radiative force:");
      fprintf_vec(stdout,"%e ",rad,3);
      fprintf(stdout,"Fragments gravitational force:");
      fprintf_vec(stdout,"%e ",gb,3);
      fprintf(stdout,"Rocket acceleration:");
      fprintf_vec(stdout,"%e ",rocket,3);
      fprintf(stdout,"Total force:");
      fprintf_vec(stdout,"%e ",dydt+k,NSTATE);
      fprintf(stdout,"Mass-loss: %e\n",dydt[6+k]);
      //exit(0);
    }
    //*/
  }
  //exit(0);

  return 0;
}

void earthObservations2(double t,double eE[],double eC[],double y[],double x[])
{
  int i,k;
  int nfrag;
  double RE[3][3];
  double xE[NSTATE],xC[NSTATE],xCE[NSTATE];
  double xobs[3],yobs[3],zobs[3];
  double zaxis[]={0,0,1};

  //CALCULATE EARTH POSITION
  stateOrbit(t,eE,xE);

  //CALCULATE COMET CM POSITION
  stateOrbit(t,eC,xC);

  //ROTATION MATRIX
  vsubg_c(xC,xE,NSTATE,xCE);
  vpack_c(-xCE[0],-xCE[1],-xCE[2],zobs);
  vcrss_c(zaxis,zobs,xobs);
  twovec_c(xobs,XAXIS,zobs,ZAXIS,RE);
  
  //ROTATE POSITIONS
  vsubg_c(y,xE,NSTATE,xCE);
  rotateState(RE,xCE,x);

  //MASS
  x[NSTATE-1]=y[NSTATE-1];
  /*
  fprintf(stdout,"Input state:");fprintf_vec(stdout,"%e ",y,NSTATE);
  fprintf(stdout,"Output state:");fprintf_vec(stdout,"%e ",x,NSTATE);
  exit(0);
  //*/
}

void earthObservations(const char* spkid,double t,double y[],double x[])
{
  int i,k;
  int nfrag;
  double RE[3][3];
  double xE[NSTATE],xC[NSTATE],xCE[NSTATE];
  double xobs[3],yobs[3],zobs[3];
  double zaxis[]={0,0,1};

  //CALCULATE EARTH POSITION
  stateObject("EARTH",t,xE);

  //CALCULATE COMET CM POSITION
  stateObject(spkid,t,xC);

  //ROTATION MATRIX
  vsubg_c(xC,xE,NSTATE,xCE);
  vpack_c(-xCE[0],-xCE[1],-xCE[2],zobs);
  vcrss_c(zaxis,zobs,xobs);
  twovec_c(xobs,XAXIS,zobs,ZAXIS,RE);
  
  //ROTATE POSITIONS
  vsubg_c(y,xE,NSTATE,xCE);
  rotateState(RE,xCE,x);

  //MASS
  x[NSTATE-1]=y[NSTATE-1];
  /*
  fprintf(stdout,"Input state:");fprintf_vec(stdout,"%e ",y,NSTATE);
  fprintf(stdout,"Output state:");fprintf_vec(stdout,"%e ",x,NSTATE);
  exit(0);
  //*/
}

double radiusGenerate(double rmin,double rmax,double q)
{
  double u,r,e,eo,p,lambda;
  p=-q;
  do{
    u=randReal();
    if(p>1){
      lambda=p-1;
      e=-1/lambda*log(1-u);
    }else{
      eo=log10(rmax/rmin);
      e=eo*u;
    }
    r=rmin*pow(10,e);
  }while(r>=rmax);
  return r;
}
       
double radiusNumber(double Rmin,double Rmax,double no,double q)
{
  double N,alpha;
  
  alpha=-(q+1);
  N=no/alpha*(1/pow(Rmin,alpha)-1/pow(Rmax,alpha));
  
  return N;
}

double radiusOveri(int i,double no,double q,double Rmax)
{
  double alpha=-(1+q);
  double uR,R;
  
  uR=(i+1)*alpha/no+1/pow(Rmax,alpha);
  R=pow(1/uR,1/alpha);

  return R;
}

double massOverN(double no,void *params)
{
  double *ps=(double*)params;
  
  int N=(int)ps[0];
  double q=ps[1];
  double Rmin=ps[2];
  double Rmax=ps[3];

  int i;
  double Ri,RN,EN,m=0;
  double beta=-(3+q);

  for(i=0;i<N;i++){
    Ri=radiusOveri(i,no,q,Rmax);
    m+=Ri*Ri*Ri;
    //fprintf(stdout,"Fragment %d: Ri = %e, Mc = %e (cum. %e)\n",i,Ri,Ri*Ri*Ri,m);
  }
  RN=radiusOveri(N,no,q,Rmax);
  EN=no/(1-beta)*(pow(RN,1-beta)-pow(Rmin,1-beta));
  ps[4]=EN;
  m+=EN;

  return m-1;
}

void noFromfragments(int N,double q,double Rmin,double Rmax,double *no,double *mlarge)
{
  double nomin,nomax,deltam;
  double pars[]={N,q,Rmin,Rmax,0};
  gsl_root_fsolver *s=gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
  gsl_function function={&massOverN,&pars};
  gsl_root_fsolver_set(s,&function,0,10);
 
  do{
    gsl_root_fsolver_iterate(s);
    gsl_root_fsolver_root(s);
    nomin=gsl_root_fsolver_x_lower(s);
    nomax=gsl_root_fsolver_x_upper(s);
    *no=gsl_root_fsolver_root(s);
    deltam=massOverN(*no,pars);
    //fprintf(stdout,"nomin = %e, nomax = %e, deltam = %e\n",nomin,nomax,deltam);
  }while(gsl_root_test_residual(deltam,1E-5)!=GSL_SUCCESS);
  deltam=massOverN(*no,pars);
  *mlarge=pars[4];
}
