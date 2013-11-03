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
#define DAY (GSL_CONST_MKSA_DAY) //s
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
real RHODUST,RHOCOMET,FR;
int NPARTICLES,NLARGE,NDEBRIS;

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
enum Fragment {_FRAGMENT_,LARGE,DEBRIS};

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

void orbitNormal(double e[],double n[])
{
  double x1[NSTATE],x2[NSTATE];
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

  double Mp;
  const double *r,*v;
  const double *ri,*rj;
  double *dyjdt,*dyidt;
  double Rij[3];
  double beta,D,Dsoft,epsoft,D3,M,R;
  double gs[]={0,0,0},rad[]={0,0,0},evap[]={0,0,0},gb[]={0,0,0};
  double ur[3],ut[3],vr[3],vt[3],pf[3],vrmag=0,vtmag=0,fr=1,ft=0,gmag=0;
  bool qcolision;

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

    //AVOID INTEGRATION OF PARTICLES ALREADY ABSORBED
    if(i>0&&Ms[i-1]==0){
      continue;
    }

    qcolision=false;
    k=NSTATE*i;

    //DISTANCE TO SUN
    D=vnorm_c(y+k);D3=D*D*D;

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
    //GRAVITATIONAL FORCE FROM OTHER BODIES
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

	//ANNOUNCE IT
	fprintf(stdout,"\t\tParticle %d has collided with particle %d\n",i-1,j-1);
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
    //EVAPORATION RECOIL
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #ifdef EVAP_FORCE
    evap[0]=0.0;
    evap[1]=0.0;
    evap[2]=0.0;
    #else
    memset(evap,0,3*sizeof(evap[0]));
    #endif

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //TOTAL
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dydt[3+k]=gs[0]+rad[0]+gb[0]+evap[0];
    dydt[4+k]=gs[1]+rad[1]+gb[1]+evap[1];
    dydt[5+k]=gs[2]+rad[2]+gb[2]+evap[2];

    /*
    if(i>nlarge){
      fprintf(stdout,"Particle: %d\n",i);
      fprintf(stdout,"Sun gravitatinal force:");
      fprintf_vec(stdout,"%e ",gs,3);
      fprintf(stdout,"Radiative force:");
      fprintf_vec(stdout,"%e ",rad,3);
      fprintf(stdout,"Fragments gravitational force:");
      fprintf_vec(stdout,"%e ",gb,3);
      fprintf(stdout,"Evaporation recoil force:");
      fprintf_vec(stdout,"%e ",evap,3);
      fprintf(stdout,"Total force:");
      fprintf_vec(stdout,"%e ",dydt+k,NSTATE);
      exit(0);
    }
    //*/
  }

  return 0;
}

void earthObservations(double t,double eE[],double eC[],double y[],double x[])
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

