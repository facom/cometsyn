#define RANSEED 1
#include <util.cpp>

int main(int argc,char *argv[])
{
  //////////////////////////////////////////
  //VARIABLES
  //////////////////////////////////////////
  int i,j,k,l;
  double xaxis[]={1,0,0},yaxis[]={0,1,0},zaxis[]={0,0,1};
  double xE[6];
  double t,tini,tint,dt;
  double xcm[6],axis[]={0,0,0,0,0,1};
  file faxis,ffrag,fint,fearth,fobs,fdir;
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

  double mp,Rp;
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
  double dirs[12];

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
  //COMETARY ORBITAL ELEMENTS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  eC[RP]=1.0;
  eC[ECC]=0.8;
  eC[INC]=0.0*D2R;
  eC[LNODE]=0.0*D2R;
  eC[ARGP]=0.0*D2R;
  eC[M0]=0.0;
  eC[MU]=GPROG*1.0;
  eC[TPER]=0.0;
  eC[PER]=tOrbit(eC);
  eC[AXISZ]=45.0*D2R;
  eC[AXISXY]=90.0*D2R;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //COMETARY PHYSICAL PROPERTIES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  RHODUST=2E3 /*kg/m^3*//(UM/(UL*UL*UL));
  Rc=10E3 /*m*/ /UL;
  Prot=17*HOURS /*secs*/ /UT;
  Mc=4*PI/3*Rc*Rc*Rc*RHODUST; 
  nlarge=20; //Number of large particles
  ndebris=80; //Number of debris particles
  if(nlarge>0) mp=Mc/nlarge; //Average mass of large fragments
  else mp=0;
  Rp=pow((mp/RHODUST)/(4*PI/3),1./3);

  fprintf(stdout,"Cometary properties:\n");
  fprintf(stdout,"\tCore radius: %e UL = %e km\n",Rc,Rc*UL/1E3);
  fprintf(stdout,"\tRotation period: %e UT = %e h\n",Prot,Prot*UT/3600);
  fprintf(stdout,"\tAverage density: %e UM/UL^3 = %e kg/m^3\n",RHODUST,RHODUST*UM/(UL*UL*UL));
  fprintf(stdout,"\tMass: %e UM = %e ton\n",Mc,Mc*UM/1E3);
  fprintf(stdout,"\tLarge fragments: %d\n",nlarge);
  fprintf(stdout,"\tAverage mass of fragments: %e UM = %e ton\n",mp,mp*UM/1E3);
  fprintf(stdout,"\tAverage radius of fragments: %e UL = %e km\n",Rp,Rp*UL/1E3);
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //DERIVED PHYSICAL PROPERTIES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  double vesc,vtan;

  vesc=sqrt(2*GPROG*Mc/Rc);
  vtan=2*PI*Rc/Prot;

  fprintf(stdout,"Cometary derived properties:\n");
  fprintf(stdout,"\tComet escape velocity: %e UL/UT = %e km/s\n",vesc,vesc*UL/UT/1E3);
  fprintf(stdout,"\tTangential velocity: %e UL/UT = %e km/s\n",vtan,vtan*UL/UT/1E3);
  exit(0);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //INTEGRATION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tini=-1.0;//yrs before periapse
  tint=1.5; //yrs 
  dt=eC[PER]/500;

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
  nfrag=nlarge+ndebris;

  //INFORMATION ABOUT FRAGMENTS
  Ms=realAlloc(nfrag);
  Rs=realAlloc(nfrag);
  type=realAlloc(nfrag);//Type of fragment: 1-large, 2-debris
  xs=stateAlloc(nfrag);

  fprintf(stdout,"Fragments:\n");
  fprintf(stdout,"\tNumber: %d\n",nfrag);
  fprintf(stdout,"\tAverage mass of large fragments: %e\n",mp);
  
  fprintf(stdout,"\n");

  ilarge=-1;
  for(i=0;i<nfrag;i++){
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
      Rs[i]=1.0E-4 /*m*/ /UL;
      Ms[i]=4*PI/3*RHODUST*Rs[i]*Rs[i]*Rs[i];
    }
    fprintf(stdout,"\t\tType (1:large,2:debris): %d\n",(int)type[i]);
    fprintf(stdout,"\t\tMass: %e UM = %e kg\n",Ms[i],Ms[i]*UM);
    fprintf(stdout,"\t\tRadius: %e UL = %e m = %e Rc\n",Rs[i],Rs[i]*UL,Rs[i]/Rc);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //INITIAL POSITION
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qfrag=false;
    k=0;
    do{
      //Generate position
      fprintf(stdout,"\t\t\tGenerating %d random position...\n",k);
      r=1.5*Rc*pow(randReal(),1./3);
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
    fprintf(stdout,"\t\tSpherical coordinates (trials %d): (%e,%e,%e)\n",
	    k,r/Rc,theta,phi);
    fprintf(stdout,"\t\tCartesian coordinates: ");
    fprintf_vec(stdout,"%e ",xs[i],3);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //INITIAL VELOCITY
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rp=r*sin(theta);
    fprintf(stdout,"\t\tAxis distance: %e m\n",rp*UL);

    //TANGENTIAL VELOCITY
    v=2*PI*rp/Prot;
    fprintf(stdout,"\t\tTangential velocity: %e UL/UT = %e m/s\n",v,v*UL/UT);
    
    //CARTESIAN VELOCITY
    vpack_c(v*cos(phi),
	    v*sin(phi),
	    0.0,
	    xs[i]+3);
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

  file fpfrag=fopen("fragments.gph","w");
  fprintf(fpfrag,
	  "file='comet-fragments.dat'\n"
	  "title='t = %.2f yrs'\n"
	  "nfrag=%d\n"
	  "nlarge=%d\n"
	  "ndebris=%d\n"
	  "tini=%e\n"
	  "UM=%e\n"
	  "UL=%e\n"
	  "UT=%e\n"
	  "Rmax=%e\n"
	  "Rc=%e\n"
	  "rf=%e\n",tini,nfrag,nlarge,ndebris,tini,UM,UL,UT,Rc*UL,Rc*UL,Rp*UL);
  fclose(fpfrag);

  //////////////////////////////////////////
  //INITIAL CONDITION COMET
  //////////////////////////////////////////
  faxis=fopen("comet-axis.dat","w");
  //STATE
  stateOrbit(tini,eC,xcm);
  
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
  for(i=0;i<nfrag;i++){
    //HELIOCENTRIC RF
    rotateState(IRHA,xs[i],xs[i]);
    //TRANSLATION
    stateAdd(xs[i],xcm,xs[i]);
  }

  //////////////////////////////////////////
  //PREPARE SYSTEM
  //////////////////////////////////////////
  nsys=6*(nfrag+1);
  y=realAlloc(nsys);
  xobs=realAlloc(nsys);
  dy=realAlloc(nsys);
  dydt=realAlloc(nsys);

  //FEED THE SYSTEM STATE ARRAY
  memcpy(y,xcm,6*sizeof(double));
  for(i=0;i<nfrag;i++){
    k=6*(i+1);
    memcpy(y+k,xs[i],6*sizeof(double));
  }
  
  //////////////////////////////////////////
  //INTEGRATE ORBITS
  //////////////////////////////////////////
  double* params=realAlloc(nfrag+2);
  params[0]=nsys;
  params[1]=nlarge;
  memcpy(params+2,Rs,nfrag*sizeof(double));

  gsl_odeiv2_system sys={gravSystem,NULL,nsys,params};
  gsl_odeiv2_driver* driver=
    gsl_odeiv2_driver_alloc_y_new(&sys,
				  gsl_odeiv2_step_rk4,
				  dt,0.0,1E-6);

  fint=fopen("comet-orbit.dat","w");
  fearth=fopen("comet-earth.dat","w");
  fobs=fopen("comet-observations.dat","w");
  fdir=fopen("comet-dirs.dat","w");

  t=0;
  fprintf(stdout,"Integration:\n");
  do{
    fprintf(stdout,"\tt = %e\n",t);

    //HELIOCENTRIC ORBIT
    fprintf(fint,"%-23.17e ",t);
    fprintf_vec(fint,"%-+23.17e ",y,nsys,false);

    //EARTH POSITION
    stateOrbit(t+tini,eE,xE);
    fprintf(fearth,"%-+23.17e ",t);
    fprintf_vec(fearth,"%-+23.17e ",xE,6,false);

    //OBSERVATIONS
    for(i=0;i<=nfrag;i++){
      k=6*i;
      earthObservations(t+tini,eE,eC,y+k,xobs+k);
    }
    //CONVERT TO PHYSICAL UNITS
    for(i=0;i<=nfrag;i++){
      k=6*i;
      xobs[0+k]*=UL;
      xobs[1+k]*=UL;
      xobs[2+k]*=UL;
      xobs[3+k]*=UL/UT;
      xobs[4+k]*=UL/UT;
      xobs[5+k]*=UL/UT;
    }
    fprintf(fobs,"%-+23.17e ",t*UT/YEAR);
    fprintf_vec(fobs,"%-+23.17e ",xobs,nsys,false);

    //DIRECTIONS
    vpack_c(0,0,1,dirs);
    mxv_c(IRHA,dirs,dirs);
    vadd_c(dirs,y,dirs);
    earthObservations(t+tini,eE,eC,dirs,dirs);

    vpack_c(-y[0],-y[1],-y[2],dirs+6);
    vadd_c(dirs+6,y,dirs+6);
    earthObservations(t+tini,eE,eC,dirs+6,dirs+6);

    fprintf(fdir,"%-+23.17e ",t);
    fprintf_vec(fdir,"%-+23.17e ",dirs,12,false);
    
    //INTEGRATE
    gsl_odeiv2_driver_apply(driver,&t,t+dt,y);
  }while(t<tint);

  //////////////////////////////////////////
  //CLOSING COMMANDS
  //////////////////////////////////////////
  fclose(fint);
  fclose(fearth);
  fclose(fobs);

  return 0;
}
