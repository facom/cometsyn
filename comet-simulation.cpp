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
//////////////////////////////////////////////////////////////////////////////////
// BEHAVIOR
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//**************************************************
//BODY ID
//**************************************************
#define COMETID ISONID

//**************************************************
//BEHAVIOR
//**************************************************
//#define COMET_PROPERTIES_ONLY
//#define SHOW_COLISION

//DRYING FUNCTION
//#define DRYING_NULL
//#define DRYING_CONSTANT
#define DRYING_EXPONENTIAL

//**************************************************
//FORCE
//**************************************************
#define SUN_FORCE /*FORCE FROM THE SUN*/
#define FRAG_FORCE /*FORCE FROM LARGE FRAGMENTS*/
//#define RADIATION_FORCE /*RADIATION FORCE*/
#define PR_CORRECTION /*POYINTING-ROBERTSON CORRECTION*/
#define ROCKET_FORCE /*EVAPORATION RECOIL*/
#define RMIN_ROCKET 1.0 /*m*/ /*MINIMUM SIZE FOR ROCKET FORCE*/

//**************************************************
//INITIAL CONDITIONS
//**************************************************
#define RANSEED 1
//#define RADIALMODE /*DEBRIS HAVE A RADIAL VELOCITY*/

//**************************************************
//DYNAMICS
//**************************************************
#define ALLOW_COLISION /*ALLOW COLISION*/

//**************************************************
//OUTPUT
//**************************************************
//#define SAVE_TRAJECTORIES

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// PROGRAM
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
#include <cometsyn.cpp>
int main(int argc,char *argv[])
{
  //////////////////////////////////////////
  //VARIABLES
  //////////////////////////////////////////
  FILE** ftraj;
  double rmax,vmax;
  double vesc,vtan,Pmin,fc;
  int i,j,k,l;
  double xaxis[]={1,0,0},yaxis[]={0,1,0},zaxis[]={0,0,1};
  double xE[NSTATE],xCE[NSTATE];
  date Tini_Date,Tend_Date,Tpresent;
  double t,tini,tend,tint,dt,dtout;
  int iout,nsteps,nscreen;
  double xcm[NSTATE],axis[NSTATE];
  file faxis,ffrag,fint,fearth,fobs,fdir,fcom,fpfrag,fdate;
  double ul,um,ut,uG;
  double eC[NELEMS],normC[3],tangC[3],radC[3];
  double RHO[3][3],IRHO[3][3];
  double ROA[3][3],IROA[3][3];
  double RHA[3][3],IRHA[3][3];
  double eE[NELEMS],normE[3];
  double q,no,Rmin,Rmax;
  double Rlarge_min,Rlarge_max,Rboulder_min,Rboulder_max,Rdust_min,Rdust_max;
  double flarge,fboulder,fdust;
  double Mlarge,Mboulder,Mdust,Mtot;
  double Rc,Prot,rho,Mc;

  int nfrag;
  int nlarge;
  int ndebris;
  int nbould;
  int ndust;

  double mp,Mp,Rp,Rpmin,Rpmax;
  double dx[3];
  double* type;
  double* Ms;
  double* Rs;
  double* rhos;
  double** xs;
  double r,rp,phi,theta,v;
  int ilarge,iboulder,idust;
  bool qfrag;
  int nsys;
  double* y;
  double* xobs;
  double* dy;
  double* dydt;
  double* dirs=realAlloc(2*NSTATE);
  double dth,Rdbmin;
  double ur[3],vr[3];
  double* xf;
  double* rcom=realAlloc(NSTATE);
  double* rcm=realAlloc(NSTATE);
  double rsun,rfearth[3],rearth,laxis;
  bool lastreport=false;
  int status;
  double *params;
  double rpmin,rpmax,alpha;
  int n;
  bool qfinal;
  char date[100];
  double AxisZ,AxisXY;
  double rhom;

  //SET AXIS
  memset(axis,0,NSTATE*sizeof(axis[0]));
  axis[5]=1.0;

  //////////////////////////////////////////
  //INITIALIZE
  //////////////////////////////////////////
  ul=1.49597870691E11; //1 AU (m)
  um=1.98892E30; //1 Msun (kg)
  ut=3.1557600E7; //1 year (s)
  uG=0.0; //Compute

  cometsynInit(ul,um,ut,uG);

  //////////////////////////////////////////
  //INTEGRATION PARAMETERS
  //////////////////////////////////////////

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //TIME
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Tini_Date="11/08/2013 00:00:00.000 UTC";
  Tend_Date="12/28/2013 00:00:00.000 UTC";
  Tend_Date="11/28/2013 00:00:00.000 UTC";

  str2et_c(Tini_Date,&tini);
  str2et_c(Tend_Date,&tend);

  tini/=UT;
  tend/=UT;
  tint=tend-tini;

  dt=0.01*DAY/UT;
  dtout=0.5*DAY/UT;
  nsteps=(int)(tint/dtout);
  nscreen=ceil(nsteps/10);
  
  fprintf(stdout,"Time configuration:\n");
  fprintf(stdout,"\tInitial date: %s = %.13e s after J2000.0\n",Tini_Date,tini*UT);
  fprintf(stdout,"\tFinal date: %s = %.13e s after J2000.0\n",Tend_Date,tend*UT);
  fprintf(stdout,"\tIntegration time: %e UT = %e s = %e days\n",tint,tint*UT,tint*UT/DAY);
  fprintf(stdout,"\tInitial integration step size: %e UT = %e s = %e days\n",dt,dt*UT,dt*UT/DAY);
  fprintf(stdout,"\tOutput step size: %e UT = %e s = %e days\n",dtout,dtout*UT,dtout*UT/DAY);
  fprintf(stdout,"\tNumber of output points: %d\n",nsteps);
  fprintf(stdout,"\tFrequency of screen report: %d\n",nscreen);

  //exit(0);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //COMETARY PHYSICAL PROPERTIES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //ORIENTATION ROTATIONAL AXIS
  AxisZ=0.0*D2R;
  AxisXY=0.0*D2R;

  //SIZE DISTRIBUTION OF DEBRIS
  rpmin=100E-6 /*m*//UL;
  rpmax=1.0E0 /*m*//UL;
  alpha=2.0; /*P(r)~r^(-alpha)*/

  //DENSITY OF ROCKY FRAGMENTS
  RHODUST=2E3 /*kg/m^3*//(UM/(UL*UL*UL));
  RHOCOMET=0.6E3 /*kg/m^3*//(UM/(UL*UL*UL));

  //FRACTION OF HEALTHY MASS IN ROCK+DUST
  /*See Greenberg (1998) "Making a cometary nucleus"*/
  FR=0.26;

  //DENSITY OF VOLATILES
  RHOVOLATILES=(1-FR)/(1/RHOCOMET-FR/RHODUST);

  //CORE TOTAL RADIUS 
  Rc=2 /*km*/ *1E3/UL;

  //SIZE OF INITIAL REGION OF FRAGMENTS AND DEBRIS
  fc=1.5;

  //PERIOD OF ROTATION
  Prot=3*HOURS /*secs*/ /UT;

  //NUMBER OF LARGE FRAGMENTS
  nlarge=2;
  flarge=1E0;

  //EXPONENT OF FRAGMENT DISTRIBUION
  q=-3.5;

  //RANGE OF DEBRIS AND FRAGMENT SIZES
  Rmin=1E-6 /*m*//UL;
  Rmax=Rc;

  //MINIMUM SIZE OF BOULDERS
  Rboulder_min=10 /*m*//UL;
  fboulder=1E-3;

  //NUMBER OF DUST TEST PARTICLES
  ndust=0;
  Rdust_min=1E-6 /*m*//UL;
  Rdust_max=1E-2 /*m*//UL;

  //THICK OF DEBRIS CRUST 
  Rdbmin=1.0*fc*Rc; /* Minimum radius for debris origin */
  dth=0.1; /* fc Rc */

  //HEALTHY CORE MASS
  Mc=4*PI/3*Rc*Rc*Rc*RHOCOMET; 

  //AVERAGE MASS AND RADIUS OF ROCKY+DUST FRAGMENTS
  if(nlarge>0) mp=Mc/nlarge; //Average mass of large fragments
  else mp=0;
  Rp=pow((mp/RHOCOMET)/(4*PI/3),1./3);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //REPORT
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(stdout,"Cometary properties:\n");
  fprintf(stdout,"\tFraction of refractory: %e\n",FR);
  fprintf(stdout,"\tDensity of rock and dust: %e UM/UL^3 = %e kg/m^3\n",
	  RHODUST,RHODUST*(UM/(UL*UL*UL)));
  fprintf(stdout,"\tDensity of volatiles: %e UM/UL^3 = %e kg/m^3\n",
	  RHOVOLATILES,RHOVOLATILES*(UM/(UL*UL*UL)));
  fprintf(stdout,"\tAverage density comet: %e UM/UL^3 = %e kg/m^3\n",RHOCOMET,RHOCOMET*UM/(UL*UL*UL));
  fprintf(stdout,"\tRotation period: %e UT = %e h\n",Prot,Prot*UT/3600);
  fprintf(stdout,"\tCore radius: %e UL = %e km\n",Rc,Rc*UL/1E3);
  fprintf(stdout,"\tMass: %e UM = %e ton\n",Mc,Mc*UM/1E3);
  fprintf(stdout,"\tLarge fragments: %d\n",nlarge);
  fprintf(stdout,"\tAverage mass of fragments: %e UM = %e ton\n",mp,mp*UM/1E3);
  fprintf(stdout,"\tAverage radius of fragments: %e UL = %e km\n",Rp,Rp*UL/1E3);
  fprintf(stdout,"\tMinimum radius of large fragments: %e UL = %e km\n",
	  Rlarge_min,Rlarge_min*UL/1E3);
  fprintf(stdout,"\tMaximum radius of large fragments: %e UL = %e km\n",
	  Rlarge_max,Rlarge_max*UL/1E3);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //DERIVED PHYSICAL PROPERTIES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  vesc=sqrt(2*GPROG*FR*Mc/Rc);
  vtan=2*PI*Rc/Prot;
  double Omegamin=sqrt(4*PI*GPROG/3*RHODUST);
  Pmin=2*PI/Omegamin;

  fprintf(stdout,"Cometary derived properties:\n");
  fprintf(stdout,"\tMaximum rotation rate for breakup: %e 1/UT = %e rad/s\n",Omegamin,Omegamin/UT);
  fprintf(stdout,"\tMinimum rotation period for breakup: %e UT = %e h\n",Pmin,Pmin*UT/3600);
  fprintf(stdout,"\tComet escape velocity: %e UL/UT = %e km/s\n",vesc,vesc*UL/UT/1E3);
  fprintf(stdout,"\tTangential velocity: %e UL/UT = %e km/s\n",vtan,vtan*UL/UT/1E3);
  
  //////////////////////////////////////////
  //COMETARY ORBIT ROTATION MATRICES
  //////////////////////////////////////////

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //NORMAL VECTOR
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  orbitDirector(COMETID,tini*UT,radC,tangC,normC);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //ROTATION MATRICES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //HELIO->ORBIT
  twovec_c(radC,XAXIS,normC,ZAXIS,RHO);
  //twovec_c(radC,ZAXIS,normC,XAXIS,RHO);
  //twovec_c(radC,XAXIS,normC,YAXIS,RHO);
  invert_c(RHO,IRHO);
  //ORBIT->AXIS
  eul2m_c(0.0*D2R,0.0*D2R,0.0*D2R,YAXIS,ZAXIS,YAXIS,ROA);
  invert_c(ROA,IROA);
  //HELIO->AXIS
  mxm_c(RHO,ROA,RHA);
  invert_c(RHA,IRHA);

  //////////////////////////////////////////
  //FRAGMENT PROPERTIES
  //////////////////////////////////////////
  noFromfragments(nlarge,q,Rmin/Rc,Rmax/Rc,&no,&Mlarge);
  Mlarge*=Mc;
  Rpmax=radiusOveri(0,no,q,Rmax/Rc)*Rc;
  Rpmin=radiusOveri(nlarge,no,q,Rmax/Rc)*Rc;
  fprintf(stdout,"\tFragment distribution parameters: no = %e, q = %e\n",no,q);
  fprintf(stdout,"\tMass in large fragments: %e UM = %e kg = %e Mc\n",Mlarge,Mlarge*UM,Mlarge/Mc);
  fprintf(stdout,"\tLargest fragment: %e UL = %e m = %e Rc\n",
	  Rpmax,Rpmax*UL,Rpmax/Rc);
  fprintf(stdout,"\tSmallest fragment: %e UL = %e m = %e Rc\n",
	  Rpmin,Rpmin*UL,Rpmin/Rc);

  //NUMBER OF BOULDERS
  Rboulder_max=Rpmin;
  nbould=fboulder*radiusNumber(Rboulder_min/Rc,Rboulder_max/Rc,no,q); 
  nbould=1;

  //NUMBER OF DEBRIS
  ndebris=nbould+ndust;

  //COUNTING FRAGMENTS
  NPARTICLES=nfrag=nlarge+ndebris;
  NLARGE=nlarge;
  NDEBRIS=ndebris;

  //INFORMATION ABOUT FRAGMENTS
  //Mass
  Ms=realAlloc(nfrag);
  //Radii
  Rs=realAlloc(nfrag);
  //Radii
  rhos=realAlloc(nfrag);
  //Type
  type=realAlloc(nfrag);//Type of fragment: 1-large, 2-boulders, 3-dust
  //State
  xs=stateAlloc(nfrag);

  fprintf(stdout,"Fragments:\n");
  fprintf(stdout,"\tLarge Fragments: %d\n",nlarge);
  fprintf(stdout,"\tDebris (boulders): %d (real / %.1f)\n",nbould,1/fboulder);
  fprintf(stdout,"\tDebris (dust): %d (real = %e)\n",ndust,
	  radiusNumber(Rdust_min/Rc,Rdust_max/Rc,no,q));
  fprintf(stdout,"\tTotal: %d\n",nfrag);

  ffrag=fopen("fragments.out","w");
  //ffrag=stdout;
  fprintf(stdout,"Storing fragment information in fragments.out\n");
  fprintf(ffrag,"Fragments:\n");
  fprintf(ffrag,"\tNumber: %d\n",nfrag);
  fprintf(ffrag,"\tAverage mass of large fragments: %e\n",mp);
  fprintf(ffrag,"\n");

  vmax=rmax=0;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //CENTRAL FRAGMENT
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Mlarge=Mboulder=Mdust=Mtot=0;
  i=0;
  type[i]=1;
  Rs[i]=Rpmax;
  Ms[i]=4*PI/3*Rs[i]*Rs[i]*Rs[i]*RHOCOMET;
  //SAVE STATE
  memset(xs[i],0,NSTATE*sizeof(xs[0][0]));
  xs[i][6]=Ms[i];
  Mlarge+=Ms[i];
  Mtot+=Ms[i];
  fprintf(ffrag,"\tFragment %d:\n",i);
  fprintf(ffrag,"\t\tType (1:large,2:debris): %d\n",(int)type[i]);
  fprintf(ffrag,"\t\tMass: %e UM = %e kg = %e Mc\n",Ms[i],Ms[i]*UM,Ms[i]/Mc);
  fprintf(ffrag,"\t\tCumulative mass: %e UM = %e kg = %e Mc\n",Mtot,Mtot*UM,Mtot/Mc);
  fprintf(ffrag,"\t\tRadius: %e UL = %e m = %e Rc\n",Rs[i],Rs[i]*UL,Rs[i]/Rc);
  fprintf(ffrag,"\t\tCartesian coordinates: ");
  fprintf_vec(ffrag,"%e ",xs[i],3);
  fprintf(ffrag,"\t\tState vector: ");
  fprintf_state(ffrag,"%e ",xs[i]);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //REST OF FRAGMENTS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ilarge=iboulder=idust=0;
  for(i=1;i<nfrag;i++){
    fprintf(ffrag,"\tFragment %d:\n",i);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //PHYSICAL PROPERTIES
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i<nlarge){
      type[i]=LARGE;
      Rs[i]=radiusOveri(i-1,no,q,Rpmax/Rc)*Rc;
      rhos[i]=RHOCOMET;
      Ms[i]=4*PI/3*Rs[i]*Rs[i]*Rs[i]*rhos[i];
      Mlarge+=Ms[i];
      ilarge++;
    }
    else if(i<nlarge+nbould){
      type[i]=BOULDER;
      Rs[i]=radiusGenerate(Rboulder_min,Rboulder_max,q);
      rhos[i]=RHOCOMET;
      Ms[i]=4*PI/3*Rs[i]*Rs[i]*Rs[i]*rhos[i];
      Mboulder+=Ms[i];
      iboulder++;
    }else{
      type[i]=DUST;
      Rs[i]=radiusGenerate(Rdust_min,Rdust_max,1);
      rhos[i]=RHODUST;
      Ms[i]=4*PI/3*Rs[i]*Rs[i]*Rs[i]*rhos[i];
      Mdust+=Ms[i];
      idust++;
    }
    Mtot+=Ms[i];

    fprintf(ffrag,"\t\tType (1:large,2:boulders,3:dust): %d\n",(int)type[i]);
    fprintf(ffrag,"\t\tRadius: %e UL = %e m = %e Rc\n",Rs[i],Rs[i]*UL,Rs[i]/Rc);
    fprintf(ffrag,"\t\tDensity: %e UM/UL^3 = %e kg/m^3\n",rhos[i],rhos[i]*UM/(UL*UL*UL));
    fprintf(ffrag,"\t\tMass: %e UM = %e kg = %e Mc\n",Ms[i],Ms[i]*UM,Ms[i]/Mc);
    fprintf(ffrag,"\t\tCumulative mass: %e UM = %e kg = %e Mc\n",Mtot,Mtot*UM,Mtot/Mc);

    /*
    if(i>=nlarge) exit(0);
    continue;
    */

    //SAVE MASS
    xs[i][6]=Ms[i];
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //INITIAL POSITION
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qfrag=false;
    k=0;
    do{
      //Generate position
      fprintf(ffrag,"\t\t\tGenerating %d random position...\n",k);

      //Region where fragments start
      if(i<nlarge)
	//Inner core
	r=fc*Rc*pow(randReal(),1./3);
      else{
	//Outer crust
	do{
	  r=(1+dth)*fc*Rc*pow(randReal(),1./3);
	}while(r<Rdbmin);
      }
      phi=2*M_PI*randReal();
      theta=acos(1-2*randReal());

      //State vector: cartesian position
      vpack_c(r*sin(theta)*cos(phi),
	      r*sin(theta)*sin(phi),
	      r*cos(theta),
	      xs[i]);
      
      //Check if particle coincide with other particles
      fprintf(ffrag,"\t\t\t\tComparing with %d other particles...\n",ilarge);
      for(j=0;j<ilarge;j++){
	vsub_c(xs[i],xs[j],dx);
	if(vnorm_c(dx)<(Rs[i]+Rs[j])){
	  qfrag=true;
	  break;
	}
	else qfrag=false;
      }
      k++;
    }while(qfrag&&k<100);

    /***DEBUG***
    if(i==1){
      printf("Fragment:%d\n",i);
      vpack_c(1*Rc,0*Rc,0*Rc,xs[i]);
    }
    if(i==2){
      printf("Fragment:%d\n",i);
      vpack_c(0*Rc,1*Rc,0*Rc,xs[i]);
    }
    if(i==3){
      printf("Fragment:%d\n",i);
      vpack_c(-1*Rc,0*Rc,0*Rc,xs[i]);
    }
    //*/

    fprintf(ffrag,"\t\tSpherical coordinates (trials %d): (%e,%e,%e)\n",
	    k,r*UL,theta,phi);
    fprintf(ffrag,"\t\tCartesian coordinates: ");
    fprintf_vec(ffrag,"%e ",xs[i],3);
    rmax=r>rmax?r:rmax;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //INITIAL VELOCITY
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rp=r*sin(theta);
    fprintf(ffrag,"\t\tAxis distance: %e m\n",rp*UL);

    //TANGENTIAL VELOCITY
    v=2*PI*rp/Prot;
    fprintf(ffrag,"\t\tTangential velocity: %e UL/UT = %e m/s\n",v,v*UL/UT);
    vmax=v>vmax?v:vmax;

    //CARTESIAN VELOCITY
    vpack_c(-v*sin(phi),
	    v*cos(phi),
	    0.0,
	    xs[i]+3);

    //RADIAL VELOCITY FOR RADIALLY DRAGGED PARTICLES
    #ifdef RADIALMODE
    if(i>=nlarge){
      //RADIAL VELOCITY
      v=vesc*randReal();
      xf=xs[i];
      vhat_c(xf,ur);
      vxk_c(ur,v,vr);

      vadd_c(xs[i]+3,vr,xs[i]+3);
      /*
      fprintf(ffrag,"Particle velocity before:");fprintf_vec(ffrag,"%e ",xs[i]+3,3);
      fprintf(ffrag,"Radial component:\n");
      fprintf(ffrag,"Escape velocity:%e\n",vesc);
      fprintf(ffrag,"Radial Speed:%e\n",v);
      fprintf(ffrag,"Radial velocity:");fprintf_vec(ffrag,"%e ",vr,3);
      fprintf(ffrag,"Particle velocity after:");fprintf_vec(ffrag,"%e ",xs[i]+3,3);
      exit(0);
      //*/
    }
    #endif
    fprintf(ffrag,"\t\tState vector: ");
    fprintf_state(ffrag,"%e ",xs[i]);

  }
  fclose(ffrag);
  
  //SAVE FRAGMENT DATA
  ffrag=fopen("fragments.dat","w");
  fprintf(ffrag,"%-6s %-6s %-13s %-13s %-13s %-13s\n","#1:i","2:type","3:rhos","4:Ms","5:Rs","6-8:x");
  for(i=0;i<nfrag;i++){
    fprintf(ffrag,"%06d %-6d %-13.5e %-13.5e %-13.5e %-+23.17e %-+23.17e %-+23.17e\n",
	    i,(int)type[i],rhos[i]*UM/(UL*UL*UL),Ms[i]*UM,Rs[i]*UL,
	    (xs[i][0]-xs[0][0])*UL/1E3,(xs[i][1]-xs[0][1])*UL/1E3,(xs[i][2]-xs[0][1])*UL/1E3);
  }
  fclose(ffrag);

  fprintf(stdout,"Maximum distance: %e km\n",rmax*UL/1E3);
  fprintf(stdout,"Maximum velocity: %e m/s\n",vmax*UL/UT);

  #ifdef SAVE_TRAJECTORIES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //TRAJECTORIES PER PARTICLE
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fpfrag=fopen("fragments.gph","w");
  fprintf(fpfrag,
	  "file='comet-fragments.dat'\n"
	  "title='t = %.2f yrs'\n"
	  "facvel=1\n"
	  "Rmax=%e\n"
	  "Rc=%e\n"
	  "rf=%e\n",
	  tini,Rc*UL,Rc*UL,Rp*UL);
  fclose(fpfrag);
  system("rm -rf output/*.dat");
  char fname[100],fpart[100];
  fprintf(stdout,"Creating individual trajectory files for %d fragments...\n",nfrag);
  ftraj=(FILE**)calloc(nfrag,sizeof(FILE*));
  for(i=0;i<nfrag;i++){
    if(i<nlarge) sprintf(fpart,"fragment_large");
    else sprintf(fpart,"fragment_debris");
    sprintf(fname,"output/%s-%05d.dat",fpart,i);
    ftraj[i]=fopen(fname,"w");
  }
  #endif
  #ifdef COMET_PROPERTIES_ONLY
  exit(0);
  #endif

  //////////////////////////////////////////
  //INITIAL CONDITION COMET
  //////////////////////////////////////////
  faxis=fopen("comet-axis.dat","w");
  //STATE
  stateObject(COMETID,tini*UT,xcm);
  xcm[6]=Mc;

  //COMETARY AXIS
  copyState(axis,xcm);

  vpack_c(1,0,0,axis+3);
  mxv_c(IRHA,axis+3,axis+3);
  fprintf_state(faxis,"%-+23.17e ",axis);

  vpack_c(0,1,0,axis+3);
  mxv_c(IRHA,axis+3,axis+3);
  fprintf_state(faxis,"%-+23.17e ",axis);

  vpack_c(0,0,1,axis+3);
  mxv_c(IRHA,axis+3,axis+3);
  fprintf_state(faxis,"%-+23.17e ",axis);
  
  fclose(faxis);

  //////////////////////////////////////////
  //PLACE PARTICLES IN THE HELIOCENTRIC RF
  //////////////////////////////////////////
  ffrag=fopen("fragments.out","a");
  fprintf(ffrag,"\nParticles in Heliocentric RF:\n");
  for(i=0;i<nfrag;i++){
    //HELIOCENTRIC RF
    rotateState(IRHA,xs[i],xs[i]);
    //TRANSLATION
    stateAdd(xs[i],xcm,xs[i]);
    fprintf(ffrag,"Particle %d:",i);
    fprintf_vec(ffrag,"%-+23.17e ",xs[i],NSTATE);
  }
  fclose(ffrag);

  //////////////////////////////////////////
  //PREPARE SYSTEM
  //////////////////////////////////////////
  nsys=NSTATE*(nfrag+1);
  y=realAlloc(nsys);
  xobs=realAlloc(nsys);
  dy=realAlloc(nsys);
  dydt=realAlloc(nsys);

  //FEED THE SYSTEM STATE ARRAY
  memcpy(y,xcm,NSTATE*sizeof(double));
  for(i=0;i<nfrag;i++){
    k=NSTATE*(i+1);
    memcpy(y+k,xs[i],NSTATE*sizeof(double));
  }

  ffrag=fopen("fragments.out","a");
  fprintf(ffrag,"Initial state (nsys = %d):\n",nsys);
  fprintf_vec(ffrag,"%e ",y,nsys);
  fclose(ffrag);

  //////////////////////////////////////////
  //INITIALIZE INTEGRATOR
  //////////////////////////////////////////

  params=realAlloc(2*nfrag+2/*nsys&nlarge*/+2/*Rc&Mc*/);
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

  //////////////////////////////////////////
  //INTEGRATE
  //////////////////////////////////////////

  //POSITION OF FRAGMENTS IN HELIOCENTRIC REFERENCE SYSTEM
  fint=fopen("comet-orbit.dat","w");
  //POSITION OF FRAGMENTS IN A GEOCENTRIC REFERENCE SYSTEM
  fobs=fopen("comet-observations.dat","w");
  //POSITION OF FRAGMENTS IN A COMET-CENTRIC REFERENCE FRAME
  fcom=fopen("comet-comet.dat","w");
  //POSITION OF EARTH IN HELIOCENTRIC REFERENCE SYSTEM
  fearth=fopen("comet-earth.dat","w");
  //DIRECTION OF COMET AXIS, SUN-DIRECTED AND EARTH-DIRECTED VECTOR
  fdir=fopen("comet-dirs.dat","w");
  //DATES
  fdate=fopen("dates.dat","w");
  
  t=0;
  fprintf(stdout,"\nStarting Integration at et = %e UT = %e s\n",tini,tini*UT);
  n=0;
  qfinal=false;
  iout=1;

  do{
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //REPORT
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if((n%nscreen)==0){
      fprintf(stdout,"\tStep %d/%d: t = (%e - %e) UL = (%e - %e) s = (%e - %e) days (nparticles = %d, nlarge = %d, ndebris = %d)\n",n,nsteps,t,t+dtout,t*UT,(t+dtout)*UT,t*UT/DAY,(t+dtout)*UT/DAY,NPARTICLES,NLARGE,NDEBRIS);
      /*
      gravSystem(t,y,dydt,params);
      for(i=0;i<=nfrag;i++){
	k=NSTATE*i;
	xf=y+k;
	fprintf(stdout,"\t\tState %d:",i);
	fprintf_vec(stdout,"%-+17.11e ",y+k,NSTATE);
	fprintf(stdout,"\t\tGradiente %d:",i);
	fprintf_vec(stdout,"%-+17.11e ",dydt+k,NSTATE);
      }
      fprintf(stdout,"\t\tParameters:");
      fprintf_vec(stdout,"%-+17.11e ",params,2*nfrag+2+2);
      //*/
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
    //SAVE TRAJECTORIES
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //fprintf(stdout,"\t\tSAVING SNAPSHOT %d t = %e s = %e days.\n",iout,t,t*UT/DAY);

    #ifdef SAVE_TRAJECTORIES
    for(i=0;i<nfrag;i++){
      k=NSTATE*(i+1);
      xf=y+k;
      //DON'T SAVE IF PARTICLE HAS COLLIDED
      if(xf[6]==0) continue;
      //TRAJECTORY PER PARTICLE
      fprintf(ftraj[i],"%-23.17e ",t);
      fprintf_vec(ftraj[i],"%-+23.17e ",y,NSTATE,false,0,false);
      fprintf_vec(ftraj[i],"%-+23.17e ",y+k,NSTATE,false);
    }
    #endif

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //SAVE HELIOCENTRIC ORBIT
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fprintf(fint,"%-23.17e ",t);
    fprintf_vec(fint,"%-+23.17e ",y,nsys,false);

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //SAVE EARTH
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    stateObject("EARTH",(t+tini)*UT,xE);
    fprintf(fearth,"%-+23.17e ",t);
    fprintf_vec(fearth,"%-+23.17e ",xE,NSTATE,false);

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //SAVE COMET-CENTRIC ORBIT
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fprintf(fcom,"%-+23.17e ",t);

    //DIRECTION OF THE SUN
    memcpy(rcm,y,NSTATE*sizeof(rcom[0]));
    vxk_c(rcm,-1,rcm);
    rotateState(RHA,rcm,rcm);
    vhat_c(rcm,rcm);
    fprintf_vec(fcom,"%-+23.17e ",rcm,3,false,0,false);

    //DIRECTION OF THE EARTH
    xf=y;
    vsubg_c(xE,xf,NSTATE,xCE);
    rotateState(RHA,xCE,xCE);
    vhat_c(xCE,xCE);
    fprintf_vec(fcom,"%-+23.17e ",xCE,3,false,0,false);
      
    //CENTER OF MASS
    xf=y;
    for(i=0;i<nfrag+1;i++){
      k=NSTATE*i;
      memcpy(rcom,y+k,NSTATE*sizeof(rcom[0]));
      vsubg_c(rcom,xf,NSTATE-1,rcom);
      rotateState(RHA,rcom,rcom);
      fprintf_vec(fcom,"%-+23.17e ",rcom,NSTATE,false,0,false);
    }
    fprintf(fcom,"\n");

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //SAVE GEOCENTRIC POSITIONS
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //OBSERVATIONS
    for(i=0;i<=nfrag;i++){
      k=NSTATE*i;
      earthObservations(COMETID,(t+tini)*UT,y+k,xobs+k);
    }
    fprintf(fobs,"%-+23.17e ",t);
    fprintf_vec(fobs,"%-+23.17e ",xobs,nsys,false);

    //DIRECTIONS
    rsun=vnorm_c(y);
    vsub_c(y,xE,rfearth);
    rearth=vnorm_c(rfearth);
    //fprintf(stdout,"Rsun = %e, Rearth = %e\n",rsun,rearth);
    laxis=0.1*rearth;
    memset(dirs,0,2*NSTATE*sizeof(dirs[0]));
    vpack_c(0,0,laxis,dirs);
    //fprintf(stdout,"Dirs initial:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
    mxv_c(IRHA,dirs,dirs);
    //fprintf(stdout,"Dirs rotated:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
    //fprintf(stdout,"State:");fprintf_vec(stdout,"%e ",y,NSTATE);
    vadd_c(dirs,y,dirs);
    //fprintf(stdout,"Dirs translated:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
    earthObservations(COMETID,(t+tini)*UT,dirs,dirs);
    //fprintf(stdout,"Dirs observed:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
    vpack_c(-laxis*y[0]/rsun,-laxis*y[1]/rsun,-laxis*y[2]/rsun,dirs+NSTATE);
    //fprintf(stdout,"Dirs sun initial:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
    vadd_c(dirs+NSTATE,y,dirs+NSTATE);
    //fprintf(stdout,"Dirs sun translated:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
    earthObservations(COMETID,(t+tini)*UT,dirs+NSTATE,dirs+NSTATE);
    //fprintf(stdout,"Dirs sun observed:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
    //exit(0);
    fprintf(fdir,"%-+23.17e ",t);
    fprintf_vec(fdir,"%-+23.17e ",dirs,2*NSTATE,false);

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

  #ifdef SAVE_TRAJECTORIES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //PLOTTING SCRIPT FOR TRAJECT.
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  char fopts[100];
  file fplot=fopen("plot-trajectories.gpl","w");
  fprintf(fplot,""
	  "set view equal xyz\n"
	  "set title 'Following %d fragments'\n"
	  "udir=10\n"
	  "set xlabel 'x(km)';set ylabel 'y(km)';set zlabel 'z(km)'\n"
	  /*"set view 0,0\n"*/
	  "splot \\\n"
	  "'output/fragment_large-00000.dat' u (0*$2):(0*$3):(0*$4):(-$2/sqrt($2*$2+$3*$3+$4*$4)*udir):(-$3/sqrt($2*$2+$3*$3+$4*$4)*udir):(-$4/sqrt($2*$2+$3*$3+$4*$4)*udir) w vec not,\\\n"
	  ,NPARTICLES);
  for(i=0;i<nfrag;i++){
    //Close trajectory files
    fclose(ftraj[i]);
    //OPTIONS
    if(i<nlarge){
      sprintf(fpart,"fragment_large");
      //sprintf(fopts,"w p pt 7 ps 2 lt 1");
      sprintf(fopts,"not w lp pt 7 ps 2 lt 1");
      //sprintf(fopts,"t '%d' w lp pt 7 ps 2",i);
    }
    else{
      sprintf(fpart,"fragment_debris");
      //sprintf(fopts,"w p pt 7 ps 1 lt 3");
      sprintf(fopts,"not w lp pt 7 ps 1 lt 3");
      //sprintf(fopts,"t '%d' w lp pt 7 ps 1",i);
    }
    //Create
    fprintf(fplot,""
	    "'output/%s-%05d.dat' u (($%d-$2)*%e):(($%d-$3)*%e):(($%d-$4)*%e) %s",
	    fpart,i,NSTATE+2,UL/1E3,NSTATE+3,UL/1E3,NSTATE+4,UL/1E3,fopts);
    if(i<nfrag-1)
      fprintf(fplot,",\\\n");
  }
  fprintf(fplot,"\npause -1\n");
  fclose(fplot);
  #endif

  fpfrag=fopen("analysis.cfg","w");
  fprintf(fpfrag,
	  "nfrag=%d\n"
	  "nlarge=%d\n"
	  "ndebris=%d\n"
	  "nfrag_end=%d\n"
	  "nlarge_end=%d\n"
	  "ndebris_end=%d\n"
	  "tini=%e\n"
	  "tint=%e\n"
	  "UM=%e\n"
	  "UL=%e\n"
	  "UT=%e\n",
	  nfrag,nlarge,ndebris,
	  NPARTICLES,NLARGE,NDEBRIS,
	  tini,tint,
	  UM,UL,UT);
  fclose(fpfrag);

  //////////////////////////////////////////
  //CLOSING COMMANDS
  //////////////////////////////////////////
  fclose(fint);
  fclose(fearth);
  fclose(fobs);
  fclose(fcom);
  fclose(fdate);

  #ifdef SAVE_TRAJECTORIES
  //system("gnuplot 'plot-trajectories.gpl'");
  #endif

  return 0;
}
