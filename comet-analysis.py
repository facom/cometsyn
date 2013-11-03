"""
COMET DESINTEGRATION
ANALYSIS
"""
from matplotlib import pyplot as plt,patches as pat,cm
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
from sys import *
from os import system

#############################################################
#CONSTANTS
#############################################################
NSTATE=7
NORM=linalg.norm
DAY=86400
DIAG=linalg.eig

#############################################################
#UTIL ROUTINES
#############################################################
def title(text):
    print "*"*80,"\n",text,"\n","*"*80,"\n"

def MoI(xs):
    npart=xs.shape[0]
    I=zeros((3,3))
    for i in xrange(0,npart):
        r=xs[i]
        m=r[6]
        I[0,0]+=m*(r[1]**2+r[2]**2)
        I[1,1]+=m*(r[0]**2+r[2]**2)
        I[2,2]+=m*(r[0]**2+r[1]**2)
        I[0,1]+=m*r[0]*r[1]
        I[1,0]+=m*r[1]*r[0]
        I[0,2]+=m*r[0]*r[2]
        I[2,0]+=m*r[2]*r[0]
        I[1,2]+=m*r[1]*r[2]
        I[2,1]+=m*r[2]*r[1]
    l,v=DIAG(I)
    lmean=l.mean()
    l=l/lmean
    return l,v

#############################################################
#READ DATA
#############################################################
#//////////////////////////////////////////////////
#CONFIGURATION
#//////////////////////////////////////////////////
comconf=dict()
execfile("analysis.cfg",comconf);
nfrag=int(comconf['nfrag'])
nlarge=int(comconf['nlarge'])
ndebris=int(comconf['ndebris'])
nfrag_end=int(comconf['nfrag_end'])
nlarge_end=int(comconf['nlarge_end'])
ndebris_end=int(comconf['ndebris_end'])
tini=comconf['tini']
tint=comconf['tint']
UL=comconf['UL']
UT=comconf['UT']
UM=comconf['UM']

#//////////////////////////////////////////////////
#OBSERVATIONS
#//////////////////////////////////////////////////
try:
    obsdata=loadtxt("comet-observations.dat");
except:
    print "No data available."
    exit(1)

nobs=obsdata.shape[0]
nsys=obsdata.shape[1]

ts=obsdata[:,0]
xcm_obs=obsdata[:,1:1+NSTATE]
xs_obs=[]
for i in xrange(1,1+int(nfrag)):
    k=NSTATE*i+1
    xs_obs+=[obsdata[:,k:k+NSTATE]]
xs_obs=array(xs_obs)
#Indexes array xs_obs: i - num. particle, j - time, k - component

#//////////////////////////////////////////////////
#ORBIT
#//////////////////////////////////////////////////
orbdata=loadtxt("comet-orbit.dat");
norb=orbdata.shape[0]
ts=orbdata[:,0]
xcm_orb=orbdata[:,1:1+NSTATE]
xs_orb=[]
for i in xrange(1,1+int(nfrag)):
    k=NSTATE*i+1
    xs_orb+=[orbdata[:,k:k+NSTATE]]
xs_orb=array(xs_orb)
#Indexes array xs_orb: i - num. particle, j - time, k - component

#//////////////////////////////////////////////////
#COMET
#//////////////////////////////////////////////////
comdata=loadtxt("comet-comet.dat");
ncom=comdata.shape[0]
ts=comdata[:,0]
sun_com=comdata[:,1:4]
ear_com=comdata[:,4:7]
xcm_com=comdata[:,7:7+NSTATE]
xs_com=[]
for i in xrange(1,1+int(nfrag)):
    k=NSTATE*i+7
    xs_com+=[comdata[:,k:k+NSTATE]]
xs_com=array(xs_com)
#Indexes array xs_com: i - num. particle, j - time, k - component

#//////////////////////////////////////////////////
#EARTH DATA
#//////////////////////////////////////////////////
earthdata=loadtxt("comet-earth.dat")
xearth=earthdata[:,1:NSTATE+1]

#//////////////////////////////////////////////////
#DIRECTIONS
#//////////////////////////////////////////////////
dirdata=loadtxt("comet-dirs.dat")
zaxis=dirdata[:,1:NSTATE+1]
dsun=dirdata[:,NSTATE+1:]

#############################################################
#TASKS
#############################################################
#////////////////////////////////////////
#HELP
#////////////////////////////////////////
def Help(**args):
    title("Help")
    print """\
Comet Analysis Script

Usage: 
    
    ./comet-analysis.py <Task> <options>

Where:

<Task>: task to be performed.  
<options>: Options of the task.

"""

    if 'verbose' in args.keys():
        print """\
Available tasks:

   Help: print this help

     Options: (not available)

   Summary: print a summary of the simulation results.

     Options: (not available)

   Snapshot: get snapshot from results.

     iobs: index of snapshot (starting at 1.  0 mean last snapshot)
"""

#////////////////////////////////////////
#SUMMARY
#////////////////////////////////////////
def Summary(**args):
    if not 'notitle' in args.keys():
        title("Summary")
    print "Physical:"
    print "Units: UL = %e m, UT = %e s, UM = %e kg"%(UL,UT,UM)
    print "Initial time of integration = %e yrs = %e days"%(tini,tini*UT/DAY)
    print "Total time of integration = %e yrs = %e days"%(tint,tint*UT/DAY)
    print
    print "Fragments:"
    print "Original number of fragments: %d"%int(nfrag)
    print "\tLarge fragments: %d"%nlarge
    print "\tDebris fragments: %d"%ndebris
    print "Final number of fragments: %d"%int(nfrag_end)
    print "\tLarge fragments: %d"%nlarge_end
    print "\tDebris fragments: %d"%ndebris_end
    print
    print "Integration:"
    print "Number of integration timesteps: %d"%nobs
    print "Fields per timestep: %d"%nsys
    print 

#////////////////////////////////////////
#SNAPSHOT
#////////////////////////////////////////
def Snapshot(iobs=1,facvel=1,**args):
    """
    Plot a snapshot of the fragments positions
    """
    title("Snapshot")
    #Summary(notitle=True)
    if iobs>ncom:
        print "Maximum snapshot %d (t = %e)"%(ncom,ts[-1]+tini)
        exit(1)

    print "Plotting position of Debris Zone at Snapshot %d"%iobs
    iobs-=1

    #**************************************************
    #MOMENT OF INERTIA ANALYSIS
    #**************************************************
    xs=xs_com[:,iobs,:]
    l,v=MoI(xs)

    #**************************************************
    #CREATE FRAGMENT FILE
    #**************************************************
    t=ts[iobs]
    print "Time since integration start: t = %e days"%(t*365.25)
    f=open("comet-fragments-snapshots.dat","w")
    rmax=vmax=0
    rmax_large=zeros(3)
    rmax_debris=zeros(3)

    nlarge_snap=nlarge
    ndebris_snap=ndebris

    #**************************************************
    #SUNWARD DIRECTION
    #**************************************************
    f.write("%e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e\n"%\
                (0,sun_com[iobs,0],sun_com[iobs,1],sun_com[iobs,2],0,0,0))

    #**************************************************
    #EARTH DIRECTION
    #**************************************************
    f.write("%e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e\n"%\
                (0,ear_com[iobs,0],ear_com[iobs,1],ear_com[iobs,2],0,0,0))

    #**************************************************
    #INERTIA AXIS
    #**************************************************
    #X-AXIS
    f.write("%e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e\n"%\
                (0,l[0]*v[0,0],l[0]*v[1,0],l[0]*v[2,0],0,0,0))
    #Y-AXIS
    f.write("%e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e\n"%\
                (0,l[1]*v[0,1],l[1]*v[1,1],l[1]*v[2,1],0,0,0))
    #Z-AXIS
    f.write("%e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e\n"%\
                (0,l[2]*v[0,2],l[2]*v[1,2],l[2]*v[2,2],0,0,0))

    for i in xrange(0,int(nfrag)):
        type=1.0
        if i>=nlarge:
            type=2.0
        rf=xs_com[i,iobs,0:NSTATE]
        Mp=rf[NSTATE-1]
        #print "Fragment %d"%(i+1),rf[0:3]
        if Mp==0:
            if type==1.0:
                nlarge_snap-=1
            else:
                ndebris_snap-=1
            continue
        rfnorm=NORM(rf[0:3])
        vnorm=NORM(rf[3:6])
        rmax=max(rmax,rfnorm)
        vmax=max(vmax,vnorm)

        if type==1:
            for j in xrange(0,3):
                rmax_large[j]=max(rmax_large[j],abs(rf[j])*UL/1E3)
        else:
            for j in xrange(0,3):
                rmax_debris[j]=max(rmax_debris[j],abs(rf[j])*UL/1E3)
        
        f.write("%e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e\n"%\
                    (type,rf[0],rf[1],rf[2],rf[3],rf[4],rf[5]))
    f.close()
    nfrag_snap=nlarge_snap+ndebris_snap
    rmax*=comconf['UL']
    vmax*=comconf['UL']/comconf['UT']
    print "Total number of fragments: %d"%nfrag_snap
    print "\tNumber of large fragments in snapshot: %d"%nlarge_snap
    print "\tNumber of debris in snapshot: %d"%ndebris_snap

    print "Maximum distance: %e km"%(rmax/1E3)
    print "Maximum separation velocity: %e km/s"%(vmax/1E3)

    print "Maximum coordinates large fragments: %e,%e,%e km"%(rmax_large[0],rmax_large[1],rmax_large[2])
    print "Maximum coordinates debris: %e,%e,%e km/s"%(rmax_debris[0],rmax_debris[1],rmax_debris[2])

    #//////////////////////////////////////////
    #//SAVE PLOTTING CONFIGURATION
    #//////////////////////////////////////////
    rcm=xcm_orb[iobs,0:NSTATE]
    rearth=xearth[iobs]
    dre=rcm-rearth
    d=NORM(rcm[0:3])
    D=NORM(dre[0:3])

    f=open("fragments.gph","w")
    f.write("""\
file='comet-fragments-snapshots.dat'
title='t = %.3f days, d = %.2f AU, D = %.2f AU'
Rmax=%e
facvel=%e
"""%(t*365.25,d,D,rmax,facvel))
    f.close()

    #SAVE CONFIGURATION FILE
    comconf['nfrag_end']=nfrag_snap
    comconf['nlarge_end']=nlarge_snap
    comconf['ndebris_end']=ndebris_snap
    f=open("analysis.cfg","w")
    for key in comconf.keys():
        if key=='__builtins__':continue
        f.write("%s=%e\n"%(key,comconf[key]))
    f.close()

    #system("gnuplot plot-fragments.gpl")

#////////////////////////////////////////
#OBSERVATIONS
#////////////////////////////////////////
def Observation(iobs=1,**args):
    """
    Plot a snapshot of the fragments positions
    """
    title("Observation")
    if iobs>ncom:
        print "Maximum snapshot %d (t = %e)"%(ncom,ts[-1]+tini)
        exit(1)

    print "Observation of Debris Zone at Snapshot %d"%iobs
    iobs-=1

    #**************************************************
    #BASIC INFORMATION
    #**************************************************
    t=ts[iobs]
    print "Time since integration start: t = %e days"%(t*365.25)

    #**************************************************
    #PREPARE FIGURE
    #**************************************************
    figfile="animation/plot-snap_%05d.png"%iobs
    fig=plt.figure(figsize=(8,8))
    plt.close("all")
    ax=fig.add_axes([0.1,0.1,0.8,0.8],axisbg='w')

    #**************************************************
    #PLOT
    #**************************************************
    npart=0
    rmax=0
    ax.plot([0],[0],'k+',markersize=20)

    #GET GLOBAL PROPERTIES
    rmaxp=abs(xs_obs[:,iobs,0:2]).max()*UL
    print "Largest distance to particles: %e km"%(rmaxp/1E3)
    rmax=1E8
    expr=floor(log10(rmax))
    facr=10**expr

    D=xs_obs[:,iobs,2].mean()
    zmin=xs_obs[:,iobs,2].min()-D
    zmax=xs_obs[:,iobs,2].max()-D

    #REFERENCE SIZES
    refsize=10
    drefsize=5
    m=abs(drefsize/zmin)
    facdebris=3

    print "Average distance: %f"%D
    print "zmin = %e, zmax = %e"%(zmin,zmax)
    for i in xrange(0,int(nfrag)):
        if i<nlarge:type=1
        if i>=nlarge:type=2
        rf=xs_obs[i,iobs,0:NSTATE]
        M=rf[NSTATE-1]
        if M==0:continue
        npart+=1
        x=rf[0]
        y=rf[1]
        z=rf[2]
        d=z-D
        size=refsize+m*d
        #print "Particle %d: z = %.17e, d = %.17e, size = %f"%(i,rf[2],d,size)

        args=dict(marker='o',markeredgecolor='none')
        if type==1:
            args['color']='r'
            args['markersize']=size
        else:
            args['color']='b'
            args['markersize']=size/facdebris

        ax.plot([x*UL/facr],[y*UL/facr],**args)

    print "Plotting %d fragments..."%npart

    #**************************************************
    #AXIS
    #**************************************************
    zdir=zaxis[iobs,0:3]
    zdir=2*rmax*zdir/facr
    ax.plot([0.0,zdir[0]],[0.0,zdir[1]],'k-',linewidth=1,linestyle=':')
    ax.text(zdir[0],zdir[1],"A",color='k',fontsize=10)

    zsun=dsun[iobs,0:3]
    zsun=2*rmax*zsun/facr
    ax.plot([0.0,zsun[0]],[0.0,zsun[1]],'k-',linewidth=1,linestyle=':')
    ax.text(zsun[0],zsun[1],"S",color='k',fontsize=10)

    #**************************************************
    #DECORATION
    #**************************************************
    ax.set_title(r"Snapshot %d, $t$ = %.3f days, $r_{\rm max}$ = %.2e km, $\Delta$ = %.2f AU"%(iobs+1,t*365.25,rmaxp/1E3,abs(D)),position=(0.5,1.05))
    ax.set_xlabel(r"$\Delta x$ ($\times 10^%d$ m)"%expr)
    ax.set_ylabel(r"$\Delta y$ ($\times 10^%d$ m)"%expr)

    ax.set_xlim((-rmax/facr*1.5,rmax/facr*1.5))
    ax.set_ylim((-rmax/facr*1.5,rmax/facr*1.5))

    #**************************************************
    #SAVEFIG
    #**************************************************
    fig.savefig(figfile);

#////////////////////////////////////////
#SNAPSHOT
#////////////////////////////////////////
def Orbit(iobs=1,fac=1e6,**args):
    """
    Plot orbit of the fragments
    """
    title("Orbit")
    #Summary(notitle=True)
    if iobs>ncom:
        print "Maximum snapshot %d (t = %e)"%(ncom,ts[-1]+tini)
        exit(1)

    print "Plotting position of Debris Zone at Snapshot %d"%iobs
    iobs-=1

    #**************************************************
    #CREATE FRAGMENT FILE
    #**************************************************
    t=ts[iobs]
    print "Time since integration start: t = %e days"%(t*365.25)

    rmax=vmax=0
    nlarge_snap=nlarge
    ndebris_snap=ndebris

    #**************************************************
    #EXTRACT SNAPSHOT INFORMATION
    #**************************************************
    f=open("comet-fragments-orbit.dat","w")
    rcm=xcm_orb[iobs]
    for i in xrange(0,int(nfrag)):
        if i<nlarge:type=1
        if i>=nlarge:type=2
        rf=xs_orb[i,iobs,0:NSTATE]
        M=rf[NSTATE-1]
        if M==0:
            if type==1.0:
                nlarge_snap-=1
            else:
                ndebris_snap-=1
            continue
        rfnorm=NORM(rf[0:3])
        vnorm=NORM(rf[3:6])
        rmax=max(rmax,rfnorm)
        vmax=max(vmax,vnorm)
        #MAGNIFY DIFFERENCE
        dr=rf-rcm
        rf=rcm+dr*fac
        f.write("%e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e\n"%\
                    (type,rf[0],rf[1],rf[2],rf[3],rf[4],rf[5]))
    f.close()
    nfrag_snap=nlarge_snap+ndebris_snap
    print "Number of large: %d"%nlarge_snap
    print "Number of debris: %d"%ndebris_snap

    #**************************************************
    #SAVE CONFIGURATION FILE
    #**************************************************
    comconf['nfrag_end']=nfrag_snap
    comconf['nlarge_end']=nlarge_snap
    comconf['ndebris_end']=ndebris_snap
    f=open("analysis.cfg","w")
    for key in comconf.keys():
        if key=='__builtins__':continue
        f.write("%s=%e\n"%(key,comconf[key]))
    f.close()

    system("gnuplot plot-orbit-fragments.gpl")

#############################################################
#SELECT TASK
#############################################################
"""
Snapshot(iobs=1,facvel=1)
Observation(iobs=1)
exit(0)
Orbit(iobs=1)
exit(0)
#"""

try:
    task=argv[1]
except:
    task="Help"

try:
    options=argv[2]
except:
    options=""

#EXECUTE TASK
cmd="%s(%s)"%(task,options)
print "Task:",cmd
try:
    eval(cmd)
except:
    print "\nError: task '%s' return error."%cmd
    Help(verbose=True)
