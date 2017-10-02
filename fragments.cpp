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
  //////////////////////////////////////////
  //INITIALIZATION AND CONFIGURATION
  //////////////////////////////////////////
  #include<fragments.hpp>

  //////////////////////////////////////////
  //PHYSICAL PROPERTIES
  //////////////////////////////////////////
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

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //DERIVED PROPERTIES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fprintf(stdout,"Cometary derived properties:\n");
  fprintf(stdout,"\tMaximum rotation rate for breakup: %e 1/UT = %e rad/s\n",Omegamin,Omegamin/UT);
  fprintf(stdout,"\tMinimum rotation period for breakup: %e UT = %e h\n",Pmin,Pmin*UT/3600);
  fprintf(stdout,"\tComet escape velocity: %e UL/UT = %e km/s\n",vesc,vesc*UL/UT/1E3);
  fprintf(stdout,"\tTangential velocity: %e UL/UT = %e km/s\n",vtan,vtan*UL/UT/1E3);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //FRAGMENTS GLOBAL PROPERTIES
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //Value of No
  double no,Mlarge;
  noFromfragments(nlarge,q,Rmin/Rc,Rmax/Rc,&no,&Mlarge);
  Mlarge*=Mc;
  
  //Radius of the largest fragment
  double Rpmax=radiusOveri(0,no,q,Rmax/Rc)*Rc;
  
  //Radius of the minimum fragment
  double Rpmin=radiusOveri(nlarge,no,q,Rmax/Rc)*Rc;

  //Total area after disruption
  double Ad=Ac*no/fabs(3+q)*(1/pow(10/*m*//UL/Rc,-(3+q))-1/pow(Rpmax/Rc,-(3+q)));

  //NUMBER OF BOULDERS
  double Rboulder_max=Rpmin;
  int nbould=fboulder*radiusNumber(Rboulder_min/Rc,Rboulder_max/Rc,no,q); 

  //NUMBER OF DEBRIS
  int ndebris=nbould+ndust;

  //COUNTING FRAGMENTS
  int nfrag;
  NPARTICLES=nfrag=nlarge+ndebris;
  NLARGE=nlarge;
  NDEBRIS=ndebris;

  //Report
  fprintf(stdout,"Fragment properties:\n");
  fprintf(stdout,"\tFragment distribution parameters: no = %e, q = %e\n",no,q);
  fprintf(stdout,"\tMass in large fragments: %e UM = %e kg = %e Mc\n",Mlarge,Mlarge*UM,Mlarge/Mc);
  fprintf(stdout,"\tLargest fragment: %e UL = %e m = %e Rc\n",
	  Rpmax,Rpmax*UL,Rpmax/Rc);
  fprintf(stdout,"\tSmallest fragment: %e UL = %e m = %e Rc\n",
	  Rpmin,Rpmin*UL,Rpmin/Rc);
  fprintf(stdout,"\tTotal core area before disruption: %e m^2\n",Ac*UL*UL);
  fprintf(stdout,"\tTotal area after disruption: %e UL^2 = %e m^2 = %e Ac\n",Ad,Ad*UL*UL,Ad/Ac);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //NUMBER OF FRAGMENTS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(stdout,"Fragments:\n");
  fprintf(stdout,"\tLarge Fragments: %d\n",nlarge);
  fprintf(stdout,"\tDebris (boulders): %d (real / %.1f)\n",nbould,1/fboulder);
  fprintf(stdout,"\tDebris (dust): %d (real = %e)\n",ndust,
	  radiusNumber(Rdust_min/Rc,Rdust_max/Rc,no,q));
  fprintf(stdout,"\tTotal: %d\n",nfrag);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //PROPERTIES OF FRAGMENTS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //PROPERTIES OF FRAGMENTS
  double* Ms=realAlloc(nfrag);
  double* Rs=realAlloc(nfrag);
  double* rhos=realAlloc(nfrag);
  double* type=realAlloc(nfrag);//Type of fragment: 1-large, 2-boulders, 3-dust
  double** xs=stateAlloc(nfrag);
  
  //==============================
  //CENTRAL FRAGMENT
  //==============================
  double Afrag,Atot=0;
  double Mboulder,Mdust,Mtot;

  Mlarge=Mboulder=Mdust=Mtot=0;
  int i=0;
  type[i]=1;
  Rs[i]=Rpmax;
  rhos[i]=RHOCOMET;
  Ms[i]=4*PI/3*Rs[i]*Rs[i]*Rs[i]*RHOCOMET;

  //SAVE STATE
  memset(xs[i],0,NSTATE*sizeof(xs[0][0]));
  xs[i][6]=Ms[i];
  Mlarge+=Ms[i];
  Mtot+=Ms[i];
  Afrag=4*PI*Rs[i]*Rs[i]*UL*UL;
  Atot+=Afrag;

  FILE* ffrag=fopen("fragments.txt","w");
  fprintf(ffrag,"\tFragment %d:\n",i);
  fprintf(ffrag,"\t\tType (1:large,2:debris): %d\n",(int)type[i]);
  fprintf(ffrag,"\t\tMass: %e UM = %e kg = %e Mc\n",Ms[i],Ms[i]*UM,Ms[i]/Mc);
  fprintf(ffrag,"\t\tCumulative mass: %e UM = %e kg = %e Mc\n",Mtot,Mtot*UM,Mtot/Mc);
  fprintf(ffrag,"\t\tRadius: %e UL = %e m = %e Rc\n",Rs[i],Rs[i]*UL,Rs[i]/Rc);
  fprintf(ffrag,"\t\tArea : %e m^2\n",Afrag);
  fprintf(ffrag,"\t\tCartesian coordinates: ");
  fprintf_vec(ffrag,"%e ",xs[i],3);
  fprintf(ffrag,"\t\tState vector: ");
  fprintf_state(ffrag,"%e ",xs[i]);

  //==============================
  //REST OF FRAGMENTS
  //==============================
  int ilarge=0,iboulder=0,idust=0;
  double rmax=0,vmax=0;
  for(int i=1;i<nfrag;i++){

    fprintf(ffrag,"\tFragment %d:\n",i);

    //------------------------------
    //PHYSICAL PROPERTIES
    //------------------------------
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
    Afrag=4*PI*Rs[i]*Rs[i]*UL*UL;
    Atot+=Afrag;

    fprintf(ffrag,"\t\tType (1:large,2:boulders,3:dust): %d\n",(int)type[i]);
    fprintf(ffrag,"\t\tRadius: %e UL = %e m = %e Rc\n",Rs[i],Rs[i]*UL,Rs[i]/Rc);
    fprintf(ffrag,"\t\tArea : %e m^2\n",Afrag);
    fprintf(ffrag,"\t\tDensity: %e UM/UL^3 = %e kg/m^3\n",rhos[i],rhos[i]*UM/(UL*UL*UL));
    fprintf(ffrag,"\t\tMass: %e UM = %e kg = %e Mc\n",Ms[i],Ms[i]*UM,Ms[i]/Mc);
    fprintf(ffrag,"\t\tCumulative mass: %e UM = %e kg = %e Mc\n",Mtot,Mtot*UM,Mtot/Mc);

    //SAVE MASS
    xs[i][6]=Ms[i];

    //------------------------------
    //INITIAL POSITION
    //------------------------------
    int qfrag=false;
    int k=0;
    double r,phi,theta,dx[3],vr[3],ur[3];
    double *xf;
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
      for(int j=0;j<ilarge;j++){
	vsub_c(xs[i],xs[j],dx);
	if(vnorm_c(dx)<(Rs[i]+Rs[j])){
	  qfrag=true;
	  break;
	}
	else qfrag=false;
      }
      k++;
    }while(qfrag&&k<100);

    fprintf(ffrag,"\t\tSpherical coordinates (trials %d): (%e,%e,%e)\n",
	    k,r*UL,theta,phi);
    fprintf(ffrag,"\t\tCartesian coordinates: ");
    fprintf_vec(ffrag,"%e ",xs[i],3);
    rmax=r>rmax?r:rmax;

    //------------------------------
    //INITIAL VELOCITY
    //------------------------------
    double rp=r*sin(theta);
    fprintf(ffrag,"\t\tAxis distance: %e m\n",rp*UL);

    //TANGENTIAL VELOCITY
    double v=2*PI*rp/Prot;
    fprintf(ffrag,"\t\tTangential velocity: %e UL/UT = %e m/s\n",v,v*UL/UT);
    vmax=v>vmax?v:vmax;

    //CARTESIAN VELOCITY
    vpack_c(-v*sin(phi),
	    v*cos(phi),
	    0.0,
	    xs[i]+3);

    //------------------------------
    //RADIAL VELOCITY 
    //------------------------------
    #ifdef RADIALMODE
    if(i>=nlarge){
      //RADIAL VELOCITY
      v=vesc*randReal();
      xf=xs[i];
      vhat_c(xf,ur);
      vxk_c(ur,v,vr);

      vadd_c(xs[i]+3,vr,xs[i]+3);
    }
    #endif
    fprintf(ffrag,"\t\tState vector: ");
    fprintf_state(ffrag,"%e ",xs[i]);

  }
  fprintf(ffrag,"\tTotal mass after disruption: %e UM = %e kg\n",Mtot,Mtot*UM);
  fprintf(ffrag,"\tTotal area after disruption: %e m^2\n",Atot);
  fclose(ffrag);
  
  //SAVE FRAGMENT DATA
  ffrag=fopen("fragments.dat","w");
  fprintf(ffrag,"%-6s %-8s%-15s%-13s%-13s%-75s%-75s\n","#1:i","2:type","3:rhos(kg/m3)","4:Ms(kg)","5:Rs(m)","6-8:x(m)","9-11:v(m/s)");
  for(i=0;i<nfrag;i++){
    fprintf(ffrag,"%06d %-8d%-15.5e%-13.5e%-13.5e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e\n",
	    i,(int)type[i],rhos[i]*UM/(UL*UL*UL),Ms[i]*UM,Rs[i]*UL,
	    (xs[i][0]-xs[0][0])*UL,(xs[i][1]-xs[0][1])*UL,(xs[i][2]-xs[0][1])*UL,
	    xs[i][3]*UL/UT,xs[i][4]*UL/UT,xs[i][5]*UL/UT
	    );
  }
  fclose(ffrag);

  fprintf(stdout,"Maximum distance: %e km\n",rmax*UL/1E3);
  fprintf(stdout,"Maximum velocity: %e m/s\n",vmax*UL/UT);

  return 0;
}
