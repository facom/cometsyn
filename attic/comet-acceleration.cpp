#include <cometsyn.cpp>

int main(int argc,char *argv[])
{
  double u,Z,Q;
  double r,A,R,M,v;
  double rho,I,MoI;
  double tau,Pini,Pend,omegaini,alpha,omegaend,omega;
  double deltat,deltaw;
  FILE *fl;
  date Tini_Date,Tend_Date;
  double tini,tend,t,dt;
  double rc[NSTATE],rmed;

  //==================================================
  //INITIALIZATION
  //==================================================
  double ul=1.49597870691E11; //1 AU (m)
  double um=1.98892E30; //1 Msun (kg)
  double ut=365.25*DAY; //1 year (s)
  double uG=0.0; //Compute
  cometsynInit(ul,um,ut,uG);

  Tini_Date="11/22/2013 00:43:00.000 UTC";
  Tend_Date="11/28/2013 00:43:00.000 UTC";
  str2et_c(Tini_Date,&tini);
  str2et_c(Tend_Date,&tend);
  deltat=tend-tini;

  t=tini;
  stateObject(ISONID,t,rc);
  r=vnorm_c(rc)*UL;
  v=vnorm_c(rc+3)*UL/UT;

  printf("Date: %s\n",Tini_Date);
  printf("Distance [AU]: %e\n",r/AU);
  printf("Velocity [km/s]: %e\n",v/1E3);
  
  //==================================================
  //COMET PROPERTIES
  //==================================================
  r=0.95*AU;
  rho=2.6E3;
  R=0.7 /*km*/*1E3;
  MoI=2./5;
  Pini=4*HOURS;

  M=4./3*PI*R*R*R*rho;
  omegaini=2*PI/Pini;
  I=MoI*M*R*R;

  printf("Comet properties:\n");
  printf("r [m] = %e\n",r);
  printf("rho [kg/m^3] = %e\n",rho);
  printf("R [m] = %e\n",R);
  printf("M [kg] = %e\n",M);
  printf("I [kg m^2] = %e\n",I);
  printf("Pini [s] = %e\n",Pini);
  printf("Omegaini [1/s] = %e\n",omegaini);
 
  //==================================================
  //SIMPLE APPROXIMATION
  //==================================================
  u=outflowVelocity(r);
  A=4*PI*R*R/4;
  Z=sublimationRate(r);
  Q=A*Z;
  tau=Q*MH2O*u*R;
  alpha=tau/I;

  printf("\nOutflow properties:\n");
  printf("Q [molec/s] = %e\n",Q);
  printf("u [m/s] = %e\n",u);
  printf("tau [N m] = %e\n",tau);
  printf("alpha [1/s^2] = %e\n",alpha);

  deltaw=alpha*deltat;
  omegaend=omegaini+deltaw;
  Pend=2*PI/omegaend;

  printf("\nSpin-up:\n");
  printf("deltat [s] = %e\n",deltat);
  printf("deltaw [1/s] = %e\n",deltaw);
  printf("Omegaend [1/s] = %e\n",omegaend);
  printf("Pend [h] = %e\n",Pend/HOURS);

  //==================================================
  //ORBITAL INTEGRATION
  //==================================================
  printf("\nDynamical:\n");
  dt=deltat/100;
  printf("Deltat [s] = %e = %e days\n",deltat,deltat/DAY);
  omega=omegaini;

  fl=fopen("outflow.dat","w");
  rmed=0;
  int i=0;
  for(t=tini;t<=tend;t+=dt){
    stateObject(ISONID,t,rc);
    r=vnorm_c(rc)*AU;
    rmed+=r;
    u=outflowVelocity(r);
    A=4*PI*R*R/4;
    Z=sublimationRate(r);
    Q=A*Z;
    tau=Q*MH2O*u*R;
    alpha=tau/I;
    deltaw=alpha*dt;
    omega+=deltaw;
    printf("t = %e s, r = %e AU, omega = %e 1/s, tend = %e, Z = %e, Q = %e\n",t,r/AU,omega,tend,Z,Q);
    fprintf(fl,"%e %e %e %e\n",t,r/AU,Q,omega);
    i++;
  }
  fclose(fl);
  rmed/=(i+1);
  omegaend=omega;
  Pend=2*PI/omegaend;
  printf("Omegaend [1/s] = %e\n",omegaend);
  printf("Pend [h] = %e\n",Pend/HOURS);
  printf("<r> [AU] = %e\n",rmed/AU);
}
