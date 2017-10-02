//////////////////////////////////////////////////////////////////////////////////
//    ___                     _   __             
//   / __\___  _ __ ___   ___| |_/ _\_   _ _ __  
//  / /  / _ \| '_ ` _ \ / _ \ __\ \| | | | '_ \ 
// / /__| (_) | | | | | |  __/ |__\ \ |_| | | | |
// \____/\___/|_| |_| |_|\___|\__\__/\__, |_| |_|
//                                   |___/       
// 2013 (C) Jorge Zuluaha, zuluagajorge@gmail.com
// Instituto de Fisica / FCEN - Universidad de Antioquia
//////////////////////////////////////////////////////////////////////////////////
//DISINTEGRATE CALCULATION
//////////////////////////////////////////////////////////////////////////////////
#include <cometsyn.cpp>

int main(int argc,char *argv[])
{
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INITIALIZATION AND CONFIGURATION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #include<fragments.hpp>
  #include<simulate.hpp>

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //FRAGMENT INFORMATION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //COUNT NUMBER OF FRAGMENTS
  int nfrag=numLines("fragments.dat");
  printf("\nREADING INFORMATION ABOUT %D FRAGMENTS\n",nfrag);

  //PROPERTY VECTORS
  double* Ms=realAlloc(nfrag);
  double* Rs=realAlloc(nfrag);
  double* rhos=realAlloc(nfrag);
  int* type=intAlloc(nfrag);//Type of fragment: 1-large, 2-boulders, 3-dust
  double** xs=stateAlloc(nfrag);

  //READ PROPERTIES OF FRAGMENTS
  char linea[MAXSTR];
  char tmp[MAXSTR];
  FILE *ffrag=fopen("fragments.dat","r");
  int i=0;
  nlarge=0;
  int ndebris=0;
  while(fgets(linea,sizeof linea,ffrag)!=NULL){
    if(strstr(linea,"#")!=NULL) continue;
    printf("Fragment %d:\n",i);
    sscanf(linea,"%s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	   tmp,&type[i],&rhos[i],&Ms[i],&Rs[i],
	   &xs[i][0],&xs[i][1],&xs[i][2],
	   &xs[i][3],&xs[i][4],&xs[i][5]
	   );
    printf("\tType:%d\n",type[i]);
    printf("\tDensity,Mass,Radius:%.3e,%.3e,%.3e\n",rhos[i],Ms[i],Rs[i]);
    printf("\tState vector (m,m/s): %s\n",sprintf_vec("%.3e ",xs[i],6));
    //TRANSFORM TO PROGRAM UNITS
    xs[i][0]/=UL;xs[i][1]/=UL;xs[i][2]/=UL;
    xs[i][3]/=(UL/UT);xs[i][4]/=(UL/UT);xs[i][5]/=(UL/UT);
    xs[i][6]=Ms[i]/UM;
    printf("\tState vector (UL,UL/UT): %s\n",sprintf_vec("%.3e ",xs[i],NSTATE));
    switch(type[i]){
    case 1:
      nlarge++;
      break;
    default:
      ndebris++;
    }
    i++;
  }
  fclose(ffrag);
  comstate[6]=Mc;

  //COUNTING FRAGMENTS
  NPARTICLES=nfrag;
  NLARGE=nlarge;
  NDEBRIS=ndebris;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //TRANSFORM TO HELIOCENTRIC COORDINATES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  printf("\nTRANSFORMING TO HELIOCENTRIC COORDINATES:\n");

  //==================================================
  //COMET ROTATION MATRICES
  //==================================================
  fprintf(stdout,"Orientation parameters:\n");
  fprintf(stdout,"\tState vector (UL, UL/UT) = %s\n",sprintf_vec("%.14e",comstate,6,false,0,false));
  fprintf(stdout,"\tradC = %s\n",sprintf_vec("%.2e",radC,3,false,0,false));
  fprintf(stdout,"\ttangC = %s\n",sprintf_vec("%.2e",tangC,3,false,0,false));
  fprintf(stdout,"\tnormC = %s\n",sprintf_vec("%.2e",normC,3,false,0,false));
  fprintf(stdout,"\tRotation Helio->Orbit = \n%s",sprintf_mat("%+5.2e",RHO,"\t\t"));
  fprintf(stdout,"\tRotation Orbit->Axis = \n%s",sprintf_mat("%+5.2e",ROA,"\t\t"));
  fprintf(stdout,"\tRotation Helio->Axis = \n%s",sprintf_mat("%+5.2e",RHA,"\t\t"));

  //==================================================
  //TRANSFORM TO HELIOCENTRIC
  //==================================================
  ffrag=fopen("fragmentshe.txt","a");
  fprintf(ffrag,"\nParticles in Heliocentric RF:\n");
  for(i=0;i<nfrag;i++){
    //HELIOCENTRIC RF
    rotateState(IRHA,xs[i],xs[i]);
    //TRANSLATION
    stateAdd(xs[i],comstate,xs[i]);
    fprintf(ffrag,"Particle %d:",i);
    fprintf_vec(ffrag,"%-+23.17e ",xs[i],NSTATE);
  }
  fclose(ffrag);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //PREPARE INTEGRATION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  printf("\nPREPARING INTEGRATION:\n");
  double tint=tend-tini;
  double dt=0.01*DAY/UT;
  double dtout=0.5*DAY/UT;
  int nsteps=(int)(tint/dtout);
  int nscreen=ceil(nsteps/10);
  fprintf(stdout,"Time configuration:\n");
  fprintf(stdout,"\tInitial date: %s = %.13e s after J2000.0\n",Tini_Date,tini*UT);
  fprintf(stdout,"\tFinal date: %s = %.13e s after J2000.0\n",Tend_Date,tend*UT);
  fprintf(stdout,"\tIntegration time: %e UT = %e s = %e days\n",tint,tint*UT,tint*UT/DAY);
  fprintf(stdout,"\tInitial integration step size: %e UT = %e s = %e days\n",dt,dt*UT,dt*UT/DAY);
  fprintf(stdout,"\tOutput step size: %e UT = %e s = %e days\n",dtout,dtout*UT,dtout*UT/DAY);
  fprintf(stdout,"\tNumber of output points: %d\n",nsteps);
  fprintf(stdout,"\tFrequency of screen report: %d\n",nscreen);

  //==================================================
  //PREPARE SYSTEM
  //==================================================
  int nsys=NSTATE*(nfrag+1);
  double* y=realAlloc(nsys);
  double* dydt=realAlloc(nsys);

  //INITIAL CONDITIONS
  memcpy(y,comstate,NSTATE*sizeof(double));
  int k;
  for(i=0;i<nfrag;i++){
    k=NSTATE*(i+1);
    memcpy(y+k,xs[i],NSTATE*sizeof(double));
  }
  ffrag=fopen("fragmentsin.txt","w");
  fprintf(ffrag,"Initial state (nsys = %d):\n",nsys);
  fprintf_vec(ffrag,"%e ",y,nsys);
  fclose(ffrag);

  //==================================================
  //INITIALIZE INTEGRATOR
  //==================================================
  double* params=realAlloc(2*nfrag+2/*nsys&nlarge*/+2/*Rc&Mc*/);
  params[0]=nsys;
  params[1]=nlarge;
  memcpy(params+2,Rs,nfrag*sizeof(double));
  memcpy(params+2+nfrag,Ms,nfrag*sizeof(double));
  params[2*nfrag+2]=Rc;
  params[2*nfrag+3]=Mc;
  gsl_odeiv2_system sys={gravSystem,NULL,nsys,params};
  gsl_odeiv2_driver* driver=
    gsl_odeiv2_driver_alloc_y_new(&sys,
				  gsl_odeiv2_step_rk4,
				  dt,0.0,1E-6);

  //==================================================
  //OUTPUT FILES
  //==================================================
  FILE* fint=fopen("comet-orbit.dat","w");
  FILE* fdate=fopen("dates.dat","w");
  

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INTEGRATE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  printf("\nINTEGRATE:\n");
  fprintf(stdout,"\nStarting Integration at et = %e UT = %e s\n",tini,tini*UT);

  //==================================================
  //INTEGRATE
  //==================================================
  int iout=1;
  double Mp=0;
  double* xf;

  //Auxiliary variables

  //Counters
  int n=0;
  double t=0;
  char date[100];

  //Conditions
  bool qfinal=false;
  bool lastreport=false;

  do{
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //REPORT
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if((n%nscreen)==0){
      fprintf(stdout,"\tStep %d/%d: t = (%e - %e) UT = (%e - %e) s = (%e - %e) days (nparticles = %d, nlarge = %d, ndebris = %d)\n",
	      n,nsteps,t,t+dtout,t*UT,(t+dtout)*UT,t*UT/DAY,(t+dtout)*UT/DAY,NPARTICLES,NLARGE,NDEBRIS);
    }
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //UPDATE MASSES
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for(i=0;i<nfrag;i++){
      k=NSTATE*(i+1);
      xf=y+k;
      Mp=params[2+nfrag+i];
      if(Mp==0){
        #ifdef SHOW_COLISION
	fprintf(stdout,"\t\t\tSetting particle %d mass to zero\n",i);
	#endif
	xf[6]=0;
      }
    }

  report:
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //SAVE DATE
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    etcal_c((t+tini)*UT,48,date);
    fprintf(fdate,"%d:%-+23.17e=%s\n",n,t,date);

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //SAVE HELIOCENTRIC ORBIT
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fprintf(fint,"%-23.17e ",t);
    fprintf_vec(fint,"%-+23.17e ",y,nsys,false);
    if(qfinal){
      lastreport=true;
      break;
    }
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //INTEGRATE
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    gsl_odeiv2_driver_apply(driver,&t,t+dtout,y);
    n++;
  }while(t<tint-dtout);

  qfinal=true;
  if(!lastreport) goto report;

  fprintf(stdout,"Integration ended at t = %e UT = %e s = %e days after %d steps\n",
	  t,t*UT,t*UT/DAY,n);

  //////////////////////////////////////////
  //CLOSING COMMANDS
  //////////////////////////////////////////
  fclose(fint);
  fclose(fdate);

  return 0;
}
