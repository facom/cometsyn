//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// BEHAVIOR
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//**************************************************
//FORCE
//**************************************************
#define SUN_FORCE /*FORCE FROM THE SUN*/
#define FRAG_FORCE /*FORCE FROM LARGE FRAGMENTS*/
#define RADIATION_FORCE /*RADIATION FORCE*/
#define PR_CORRECTION /*POYINTING-ROBERTSON CORRECTION*/
#define EVAP_FORCE /*EVAPORATION RECOIL*/

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
#include <util.cpp>
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
  double t,tini,tint,dt;
  int iout,nout;
  double xcm[NSTATE],axis[NSTATE];
  file faxis,ffrag,fint,fearth,fobs,fdir,fcom,fpfrag;
  double ul,um,ut,uG;
  double eC[NELEMS],normC[3];
  double RHO[3][3],IRHO[3][3];
  double ROA[3][3],IROA[3][3];
  double RHA[3][3],IRHA[3][3];
  double eE[NELEMS],normE[3];
  double Rc,Prot,rho,Mc;

  int nfrag;
  int nlarge;
  int ndebris;

  double mp,Mp,Rp;
  double dx[3];
  double* type;
  double* Ms;
  double* Rs;
  double** xs;
  double r,rp,phi,theta,v;
  int ilarge;
  bool qfrag;
  int nsys;
  double* y;
  double* xobs;
  double* dy;
  double* dydt;
  double* dirs=realAlloc(2*NSTATE);
  double dth,Rmin;
  double ur[3],vr[3];
  double* xf;
  double* rcom=realAlloc(NSTATE);
  double* rcm=realAlloc(NSTATE);
  double rsun,rfearth[3],rearth,laxis;
  bool lastreport=false;
  int status;
  
  //SET AXIS
  memset(axis,0,NSTATE*sizeof(axis[0]));
  axis[5]=1.0;

  //////////////////////////////////////////
  //SET UNITS
  //////////////////////////////////////////
  ul=1.49597870691E11; //1 AU (m)
  um=1.98892E30; //1 Msun (kg)
  ut=3.1557600E7; //1 year (s)
  uG=0.0; //Compute
  setUnits(ul,um,ut,uG);
  randomInit();

  //////////////////////////////////////////
  //GLOBAL PARAMETERS
  //////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //COMETARY PHYSICAL PROPERTIES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //SIZE DISTRIBUTION OF DEBRIS
  double rmin=100E-6;
  double alpha=2;

  //DENSITY OF ROCKY FRAGMENTS
  RHODUST=2E3 /*kg/m^3*//(UM/(UL*UL*UL));
  
  //MEAN DENSITY OF COMET CORE WHILE HEALTHY
  RHOCOMET=0.6E3 /*kg/m^3*//(UM/(UL*UL*UL));

  //FRACTION OF HEALTHY MASS IN ROCKS+DUST
  FR=0.5;

  //CORE TOTAL RADIUS 
  Rc=10E3 /*m*/ /UL;

  //SIZE OF FRAGMENT AND DEBRIS REGION
  fc=1.5;

  //THICK OF DEBRIS CRUST 
  Rmin=1.0*fc*Rc; /* Minimum radius for debris origin */
  dth=0.1; /* fc Rc */

  //PERIOD OF ROTATION
  Prot=4*HOURS /*secs*/ /UT;

  //HEALTHY CORE MASS
  Mc=4*PI/3*Rc*Rc*Rc*RHOCOMET; 

  //NUMBER OF FRAGMENTS AND DEBRIS 
  nlarge=40; //Number of large particles
  ndebris=400; //Number of debris particles

  //AVERAGE MASS AND RADIUS OF ROCKY+DUST FRAGMENTS
  if(nlarge>0) mp=FR*Mc/nlarge; //Average mass of large fragments
  else mp=0;
  Rp=pow((mp/RHODUST)/(4*PI/3),1./3);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //COMETARY ORBITAL ELEMENTS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  eC[RP]=0.3;
  eC[ECC]=0.8;
  eC[INC]=0.0*D2R;
  eC[LNODE]=0.0*D2R;
  eC[ARGP]=0.0*D2R;
  eC[M0]=0.0*D2R;
  eC[MU]=GPROG*1.0;
  eC[TPER]=0.0;
  eC[PER]=tOrbit(eC);
  eC[AXISZ]=45.0*D2R;
  eC[AXISXY]=90.0*D2R;
  eC[EMASS]=Mc;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //REPORT
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(stdout,"Cometary properties:\n");
  fprintf(stdout,"\tAverage density: %e UM/UL^3 = %e kg/m^3\n",RHODUST,RHODUST*UM/(UL*UL*UL));
  fprintf(stdout,"\tRotation period: %e UT = %e h\n",Prot,Prot*UT/3600);
  fprintf(stdout,"\tCore radius: %e UL = %e km\n",Rc,Rc*UL/1E3);
  fprintf(stdout,"\tMass: %e UM = %e ton\n",Mc,Mc*UM/1E3);
  fprintf(stdout,"\tLarge fragments: %d\n",nlarge);
  fprintf(stdout,"\tAverage mass of fragments: %e UM = %e ton\n",mp,mp*UM/1E3);
  fprintf(stdout,"\tAverage radius of fragments: %e UL = %e km\n",Rp,Rp*UL/1E3);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //DERIVED PHYSICAL PROPERTIES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  vesc=sqrt(2*GPROG*FR*Mc/Rc);
  vtan=2*PI*Rc/Prot;
  Pmin=sqrt(3*PI/(2*GPROG*FR*RHOCOMET));

  fprintf(stdout,"Cometary derived properties:\n");
  fprintf(stdout,"\tMinimum rotation period for breakup: %e UT = %e h\n",Pmin,Pmin*UT/3600);
  fprintf(stdout,"\tComet escape velocity: %e UL/UT = %e km/s\n",vesc,vesc*UL/UT/1E3);
  fprintf(stdout,"\tTangential velocity: %e UL/UT = %e km/s\n",vtan,vtan*UL/UT/1E3);
  //exit(0);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INTEGRATION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tini=-1.0;//yrs before periapse
  tint=1.5; //yrs of integration time
  dt=eC[PER]/1000;//Time step
  /*
  dt=1/365.25/20;
  tint=100*dt;
  //*/
  /*
  tini=-100*DAY/UT;
  tint=20*DAY/UT;
  dt=0.1*DAY/UT;
  */
  nout=1;

  //////////////////////////////////////////
  //COMETARY ORBIT
  //////////////////////////////////////////
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //NORMAL VECTOR
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  orbitNormal(eC,normC);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //ROTATION MATRICES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //HELIO->ORBIT
  twovec_c(xaxis,XAXIS,normC,ZAXIS,RHO);
  invert_c(RHO,IRHO);
  //ORBIT->AXIS
  eul2m_c(eC[AXISZ],eC[AXISXY],0.0,YAXIS,ZAXIS,ZAXIS,ROA);
  invert_c(ROA,IROA);
  //HELIO->AXIS
  mxm_c(RHO,ROA,RHA);
  invert_c(RHA,IRHA);

  //////////////////////////////////////////
  //EARTH ORBITAL ELEMENTS
  //////////////////////////////////////////
  eE[RP]=1.0;
  eE[ECC]=0.016;
  eE[INC]=0.0*D2R;
  eE[LNODE]=0.0*D2R;
  eE[ARGP]=0.0*D2R;
  eE[M0]=120.0*D2R;
  eE[MU]=GPROG*1.0;
  eE[TPER]=0.0;
  eE[PER]=tOrbit(eC);
  eE[AXISZ]=23.5*D2R;
  eE[AXISXY]=0.0*D2R;
  eE[EMASS]=MEARTH/UM;
  orbitNormal(eE,normE);

  //////////////////////////////////////////
  //COMETARY PHYSICAL PROPERTIES
  //////////////////////////////////////////
  fprintf(stdout,"Comet Physical Properties:\n");
  fprintf(stdout,"\tRadius = %e UL = %e m\n",Rc,Rc*UL);
  fprintf(stdout,"\tDensity = %e UM/UL^3 = %e kg/m^3\n",RHODUST,RHODUST*UM/(UL*UL*UL));
  fprintf(stdout,"\tMass = %e UM = %e kg\n",Mc,Mc*UM);
  fprintf(stdout,"\tPeriod of rotation = %e UT = %e secs\n",Prot,Prot*UT);

  //////////////////////////////////////////
  //FRAGMENT PROPERTIES
  //////////////////////////////////////////
  ffrag=fopen("comet-fragments.dat","w");

  //TOTAL FRAGMENTS
  NPARTICLES=nfrag=nlarge+ndebris;
  NLARGE=nlarge;
  NDEBRIS=ndebris;

  //INFORMATION ABOUT FRAGMENTS
  Ms=realAlloc(nfrag);
  Rs=realAlloc(nfrag);
  type=realAlloc(nfrag);//Type of fragment: 1-large, 2-debris
  xs=stateAlloc(nfrag);

  fprintf(stdout,"Fragments:\n");
  fprintf(stdout,"\tNumber: %d\n",nfrag);
  fprintf(stdout,"\tAverage mass of large fragments: %e\n",mp);
  
  fprintf(stdout,"\n");

  vmax=rmax=0;

  //CENTRAL FRAGMENT
  i=0;
  type[i]=1;
  Ms[i]=mp;
  //SAVE MASS
  xs[i][6]=Ms[i];
  Rs[i]=pow(Ms[i]/RHODUST/(4*PI/3),1./3);
  fprintf(stdout,"\tFragment %d:\n",i);
  fprintf(stdout,"\t\tType (1:large,2:debris): %d\n",(int)type[i]);
  fprintf(stdout,"\t\tMass: %e UM = %e kg\n",Ms[i],Ms[i]*UM);
  fprintf(stdout,"\t\tRadius: %e UL = %e m = %e Rc\n",Rs[i],Rs[i]*UL,Rs[i]/Rc);
  fprintf(stdout,"\t\tCartesian coordinates: ");
  fprintf_vec(stdout,"%e ",xs[i],3);
  fprintf(stdout,"\t\tState vector: ");
  fprintf_state(stdout,"%e ",xs[i]);
  fprintf(ffrag,"%e ",type[i]);
  fprintf_state(ffrag,"%-+25.17e ",xs[i]);
  ilarge=0;

  for(i=1;i<nfrag;i++){
    fprintf(stdout,"\tFragment %d:\n",i);
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //PHYSICAL PROPERTIES
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i<nlarge){
      type[i]=LARGE;
      Ms[i]=mp;
      Rs[i]=pow(Ms[i]/RHODUST/(4*PI/3),1./3);
      ilarge++;
    }
    else{
      type[i]=DEBRIS;
      //Rs[i]=1.0E-3 /*m*/ /UL;
      /*GENERATE RANDOM SIZE*/
      Rs[i]=rmin/pow((1-randReal()),1/(alpha-1)) /*m*/ /UL;
      Ms[i]=4*PI/3*RHODUST*Rs[i]*Rs[i]*Rs[i];
    }
    fprintf(stdout,"\t\tType (1:large,2:debris): %d\n",(int)type[i]);
    fprintf(stdout,"\t\tMass: %e UM = %e kg\n",Ms[i],Ms[i]*UM);
    fprintf(stdout,"\t\tRadius: %e UL = %e m = %e Rc\n",Rs[i],Rs[i]*UL,Rs[i]/Rc);

    //SAVE MASS
    xs[i][6]=Ms[i];
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //INITIAL POSITION
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qfrag=false;
    k=0;
    do{
      //Generate position
      fprintf(stdout,"\t\t\tGenerating %d random position...\n",k);

      //Region where fragments start
      if(i<nlarge)
	//Inner core
	r=fc*Rc*pow(randReal(),1./3);
      else{
	//Outer thick
	do{
	  r=(1+dth)*fc*Rc*pow(randReal(),1./3);
	}while(r<Rmin);
      }
	  
	
      phi=2*M_PI*randReal();
      theta=acos(1-2*randReal());

      //State vector: cartesian position
      vpack_c(r*sin(theta)*cos(phi),
	      r*sin(theta)*sin(phi),
	      r*cos(theta),
	      xs[i]);
      
      //Check if particle coincide with other particles
      fprintf(stdout,"\t\t\t\tComparing with %d other particles...\n",ilarge);
      for(j=0;j<ilarge;j++){
	vsub_c(xs[i],xs[j],dx);
	if(vnorm_c(dx)<(Rs[i]+Rs[j])){
	  qfrag=true;
	  break;
	}
	else qfrag=false;
      }
      k++;
    }while(qfrag);

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

    fprintf(stdout,"\t\tSpherical coordinates (trials %d): (%e,%e,%e)\n",
	    k,r*UL,theta,phi);
    fprintf(stdout,"\t\tCartesian coordinates: ");
    fprintf_vec(stdout,"%e ",xs[i],3);
    rmax=r>rmax?r:rmax;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //INITIAL VELOCITY
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rp=r*sin(theta);
    fprintf(stdout,"\t\tAxis distance: %e m\n",rp*UL);

    //TANGENTIAL VELOCITY
    v=2*PI*rp/Prot;
    fprintf(stdout,"\t\tTangential velocity: %e UL/UT = %e m/s\n",v,v*UL/UT);
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
      fprintf(stdout,"Particle velocity before:");fprintf_vec(stdout,"%e ",xs[i]+3,3);
      fprintf(stdout,"Radial component:\n");
      fprintf(stdout,"Escape velocity:%e\n",vesc);
      fprintf(stdout,"Radial Speed:%e\n",v);
      fprintf(stdout,"Radial velocity:");fprintf_vec(stdout,"%e ",vr,3);
      fprintf(stdout,"Particle velocity after:");fprintf_vec(stdout,"%e ",xs[i]+3,3);
      exit(0);
      //*/
    }
    #endif
    fprintf(stdout,"\t\tState vector: ");
    fprintf_state(stdout,"%e ",xs[i]);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //SAVE FRAGMENT
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(ffrag,"%e ",type[i]);
    fprintf_state(ffrag,"%-+25.17e ",xs[i]);
    //scanf("%d",&l);
    //if(ilarge==nlarge-1) break;
  }
  fclose(ffrag);

  fprintf(stdout,"Maximum distance: %e\n",rmax*UL);
  fprintf(stdout,"Maximum velocity: %e\n",vmax*UL/UT/1E3);
  //exit(0);

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
  //exit(0);

  //////////////////////////////////////////
  //INITIAL CONDITION COMET
  //////////////////////////////////////////
  faxis=fopen("comet-axis.dat","w");
  //STATE
  stateOrbit(tini,eC,xcm);
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
  //fprintf(stdout,"Particles in Heliocentric RF:\n");
  for(i=0;i<nfrag;i++){
    //HELIOCENTRIC RF
    rotateState(IRHA,xs[i],xs[i]);
    //TRANSLATION
    stateAdd(xs[i],xcm,xs[i]);
    /*
    fprintf(stdout,"Particle %d:",i);
    fprintf_vec(stdout,"%-+23.17e ",xs[i],NSTATE);
    */
  }

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
  /*
  fprintf(stdout,"Initial state (nsys = %d):\n",nsys);
  fprintf_vec(stdout,"%e ",y,nsys);
  */
  //exit(0);

  //////////////////////////////////////////
  //INTEGRATE ORBITS
  //////////////////////////////////////////
  double* params=realAlloc(2*nfrag+2);
  params[0]=nsys;
  params[1]=nlarge;
  memcpy(params+2,Rs,nfrag*sizeof(double));
  memcpy(params+2+nfrag,Ms,nfrag*sizeof(double));

  gsl_odeiv2_system sys={gravSystem,NULL,nsys,params};
  gsl_odeiv2_driver* driver=
    gsl_odeiv2_driver_alloc_y_new(&sys,
				  gsl_odeiv2_step_rk4,
				  dt,0.0,1E-6);

  fint=fopen("comet-orbit.dat","w");
  fearth=fopen("comet-earth.dat","w");
  fobs=fopen("comet-observations.dat","w");
  fdir=fopen("comet-dirs.dat","w");
  fcom=fopen("comet-comet.dat","w");

  t=0;
  fprintf(stdout,"Integration:\n");
  int n=0;
  int nsteps=(int)(tint/dt);
  int nscreen=nsteps/10;
  bool qfinal=false;
  fprintf(stdout,"Number of steps:%d\n",nsteps);
  iout=1;
  do{
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //REPORT
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if((n%nscreen)==0)
      fprintf(stdout,"\tStep %d: t = %e (nparticles = %d, nlarge = %d, ndebris = %d)\n",n,t,NPARTICLES,NLARGE,NDEBRIS);

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //UPDATE MASSES
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for(i=0;i<nfrag;i++){
      k=NSTATE*(i+1);
      xf=y+k;
      Mp=params[2+nfrag+i];
      if(Mp==0) xf[6]=0;
    }
    
  report:
    if((n%nout)==0 || qfinal){
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //SAVE TRAJECTORIES
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      fprintf(stdout,"\t\tSAVING SNAPSHOT %d t = %e.\n",iout,t);
      iout++;
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
      stateOrbit(t+tini,eE,xE);
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
	earthObservations(t+tini,eE,eC,y+k,xobs+k);
      }
      fprintf(fobs,"%-+23.17e ",t*UT/YEAR);
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
      earthObservations(t+tini,eE,eC,dirs,dirs);
      //fprintf(stdout,"Dirs observed:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
      vpack_c(-laxis*y[0]/rsun,-laxis*y[1]/rsun,-laxis*y[2]/rsun,dirs+NSTATE);
      //fprintf(stdout,"Dirs sun initial:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
      vadd_c(dirs+NSTATE,y,dirs+NSTATE);
      //fprintf(stdout,"Dirs sun translated:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
      earthObservations(t+tini,eE,eC,dirs+NSTATE,dirs+NSTATE);
      //fprintf(stdout,"Dirs sun observed:");fprintf_vec(stdout,"%e ",dirs,2*NSTATE);
      //exit(0);
      fprintf(fdir,"%-+23.17e ",t);
      fprintf_vec(fdir,"%-+23.17e ",dirs,2*NSTATE,false);
    }//END SAVING STATE

    if(qfinal){
      lastreport=true;
      break;
    }
    
    /***DEBUG***
    break;
    //*/

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //INTEGRATE
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    gsl_odeiv2_driver_apply(driver,&t,t+dt,y);
    n++;
  }while(t<tint);
  qfinal=true;
  if(!lastreport) goto report;

  fprintf(stdout,"Integration ended at t = %e yrs = %e days after %d steps\n",
	  t,t*365.25,n);
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

  #ifdef SAVE_TRAJECTORIES
  //system("gnuplot 'plot-trajectories.gpl'");
  #endif

  return 0;
}
