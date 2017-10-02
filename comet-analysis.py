"""
COMET DESINTEGRATION
ANALYSIS
"""
from matplotlib import pyplot as plt,patches as pat,cm,image as img
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
from sys import *
from os import system

#############################################################
#CONSTANTS
#############################################################
NSTATE=7
NORM=linalg.norm
D2R=pi/180
DAY=86400
DIAG=linalg.eig
AU=1.496E8

#############################################################
#UTIL ROUTINES
#############################################################
def point(ax,x,y,zorder=10,**args):
    ax.plot(x,y,linestyle='none',marker='o',color='w',alpha=0.3,
            markersize=10,markeredgecolor='none',**args)
    ax.plot(x,y,linestyle='none',marker='o',color='w',alpha=0.4,
            markersize=6,zorder=zorder+1,markeredgecolor='none',**args)
    ax.plot(x,y,linestyle='none',marker='o',color='w',alpha=0.6,
            markersize=4,zorder=zorder+2,markeredgecolor='none',**args)
    ax.plot(x,y,linestyle='none',marker='o',color='w',alpha=1,
            markersize=2,zorder=zorder+3,markeredgecolor='none',**args)

def comet(ax,x,y,zorder=10,**args):
    ax.plot(1.02*x,y,linestyle='none',marker='o',color='w',alpha=0.3,
            markersize=10,markeredgecolor='none',**args)
    ax.plot(1.01*x,y,linestyle='none',marker='o',color='w',alpha=0.4,
            markersize=6,zorder=zorder+1,markeredgecolor='none',**args)
    ax.plot(1.009*x,y,linestyle='none',marker='o',color='w',alpha=0.6,
            markersize=4,zorder=zorder+2,markeredgecolor='none',**args)
    ax.plot(1.00*x,y,linestyle='none',marker='o',color='w',alpha=1,
            markersize=2,zorder=zorder+3,markeredgecolor='none',**args)

def map(x,xmax,n):
    i=n*(x+xmax)/(2*xmax)
    return i

def angsize(r,d):
    t=(r/d)/D2R*3600
    return t

def scinot(value,fmt="%.2f"):
    exponent=floor(log10(value))
    mantisa=value/10**exponent
    strval=fmt%mantisa
    str=r'$%s\times 10^{%d}$'%(strval,exponent)
    return str

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
rhoc=comconf['rhoc']
fr=comconf['fr']
q=comconf['q']
UL=comconf['UL']
UT=comconf['UT']
UM=comconf['UM']

#//////////////////////////////////////////////////
#DATES
#//////////////////////////////////////////////////
fdate=open("dates.dat","r")
dates=[]
for line in fdate:
    tphys,date=line.split("=")
    dates+=[date.strip()]
fdate.close()
#"""
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
    date=dates[iobs]
    print "Time since integration start: t = %e days"%(t*365.25)
    print "Date: %s"%date
    print 

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
    print "Maximum coordinates debris: %e,%e,%e km"%(rmax_debris[0],rmax_debris[1],rmax_debris[2])

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

    system("gnuplot plot-fragments.gpl")

#////////////////////////////////////////
#DISTRIBUTION
#////////////////////////////////////////
def Distribution(iobs=1,rmax=1E8,**args):
    """
    Plot a snapshot of the fragments positions
    """
    title("Distribution")
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
    figfile="animation/distribution-snap_%05d.png"%iobs
    fig=plt.figure(figsize=(8,8))
    plt.close("all")
    ax=fig.add_axes([0.1,0.1,0.8,0.8],axisbg='w')

    #**************************************************
    #PLOT
    #**************************************************
    npart=0
    ax.plot([0],[0],'k+',markersize=20)

    #GET GLOBAL PROPERTIES
    rmaxp=abs(xs_obs[:,iobs,0:2]).max()*UL
    print "Largest distance to particles: %e km"%(rmaxp/1E3)
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
#DISTRIBUTION
#////////////////////////////////////////
def Observation(iobs=1,field=0,**args):
    """
    Plot observation
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
    date=dates[iobs]
    print "Time since integration start: t = %e UL = %e s = %e days"%(t,t*UT,t*UT/DAY)
    print "Date: %s"%date
    print 

    #**************************************************
    #PREPARE FIGURE
    #**************************************************
    figfile="animation/picture-snap_%05d.png"%iobs
    fig=plt.figure(figsize=(8,8))
    plt.close("all")
    ax=fig.add_axes([0.1,0.1,0.8,0.8],axisbg='w')

    #**************************************************
    #PLOT
    #**************************************************
    npart=0
    ax.plot([0],[0],'k+',markersize=20)

    #GET GLOBAL PROPERTIES
    rmaxp=abs(xs_obs[:nlarge,iobs,0:2]).max()*UL
    rmedp=abs(xs_obs[:nlarge,iobs,0:2]).mean()*UL
    print "Largest distance to large particles: %e km"%(rmaxp/1E3)
    print "Average distance to particles: %e km"%(rmedp/1E3)
    print 

    rmaxp=abs(xs_obs[nlarge:,iobs,0:2]).max()*UL
    rmedp=abs(xs_obs[nlarge:,iobs,0:2]).mean()*UL
    print "Largest distance to debris: %e km"%(rmaxp/1E3)
    print "Average distance to debris: %e km"%(rmedp/1E3)
    print 

    D=xs_obs[:,iobs,2].mean()
    print "Average distance : %e UL = %e km"%(abs(D),abs(D)*UL/1E3)
    zmin=xs_obs[:,iobs,2].min()-D
    zmax=xs_obs[:,iobs,2].max()-D
    print "Depth debris zone : (%e : %e) UL = (%e : %e) km"%(zmin,zmax,
                                                             zmin*UL/1E3,zmax*UL/1E3)
    print

    facdebris=2
    dmax=0
    rmax=0
    rsun=NORM(xcm_orb[iobs,0:3])
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

        rmax=max(rmax,NORM(rf[0:2]))
        dalpha=x/abs(D)/D2R*3600
        ddelta=y/abs(D)/D2R*3600

        dmax=max(dmax,abs(dalpha),abs(ddelta))

        size=5
        args=dict(marker='o',markeredgecolor='none')
        if type==1:
            args['color']='r'
            args['markersize']=size
        else:
            args['color']='b'
            args['markersize']=size/facdebris

        ax.plot([dalpha],[ddelta],**args)

    print "Plotting %d fragments..."%npart

    #**************************************************
    #AXIS
    #**************************************************
    zdir=zaxis[iobs,0:3]
    zdir=5*dmax*zdir
    ax.plot([0.0,zdir[0]],[0.0,zdir[1]],'k-',linewidth=2,linestyle='-')
    ax.text(zdir[0],zdir[1],"A",color='k',fontsize=10,
            bbox=dict(boxstyle='square,pad=0.5',fc="1.0",ec="none"))

    zsun=dsun[iobs,0:3]
    zsun=5*dmax*zsun
    ax.plot([0.0,zsun[0]],[0.0,zsun[1]],'k-',linewidth=2,linestyle='-')
    ax.text(zsun[0],zsun[1],"S",color='k',fontsize=10,
            bbox=dict(boxstyle='square,pad=0.5',fc="1.0",ec="none"))

    #**************************************************
    #DECORATION
    #**************************************************
    ax.set_title(r"$t$ = %.2f days, r = %s AU, $\Delta$ = %s AU"%(t*UT/DAY,scinot(rsun),scinot(abs(D))),
                 position=(0.5,1.05))

    ax.text(0.98,0.95,date,
            horizontalalignment='right',
            fontsize=12,
            transform=ax.transAxes)
    
    ax.text(0.02,0.95,r"$r_{\rm max}$ = %s km"%(scinot(rmax*UL/1E3)),
            horizontalalignment='left',
            fontsize=12,
            transform=ax.transAxes)
    
    dscale=10**(floor(log10(dmax)))
    equiv=dscale/3600.0*D2R*abs(D)*UL/1E3
    ax.text(0.98,0.02,"%.0f arcsec = %s km"%(dscale,scinot(equiv)),
            horizontalalignment='right',
            fontsize=10,
            transform=ax.transAxes)

    ax.text(1.02,1.0,"Ferrin & Zuluaga (in prep. 2013)",
            verticalalignment='top',horizontalalignment='left',
            rotation=270,fontsize=10,
            transform=ax.transAxes)

    ax.set_xlabel(r"$\Delta \alpha$ (arcsec)")
    ax.set_ylabel(r"$\Delta \delta$ (arcsec)")

    dmax*=1.5
    ax.set_xlim((-dmax,dmax))
    ax.set_ylim((-dmax,dmax))

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
    date=dates[iobs]
    print "Time since integration start: t = %e days"%(t*365.25)
    print "Date: %s"%date
    print 

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

#////////////////////////////////////////
#MAXIMUM DISTANCE
#////////////////////////////////////////
def MaximumDistance(**args):
    """
    Calculate maximum distance of large fragments to center of mass of debris zone
    """
    title("Maximum Distance")

    fmax=open("maximum-distance.dat","w")
    fmax.write("%7s %-24s %-12s %-12s %-12s %-12s %-12s %-12s\n"%(\
            "#1:iobs","2:t",
            "3:davg","4:ravg","5:pavg",
            "6:dmax","7:rmax","8:pmax"))
            
    for iobs in xrange(0,ncom):    
        #TIME AND DATE
        t=ts[iobs]
        date=dates[iobs]
        print "Time since integration start: t = %e days"%(t*365.25)
        print "Date: %s"%date
        print 
        
        #GET LARGE FRAGMENT SPATIAL POSITIONS
        xs=xs_com[0:nlarge,iobs,:]
        
        #GET LARGE FRAGMENT POSITION OVER PLANE OF SKY
        ps=xs_obs[0:nlarge,iobs,:]

        #GET LARGER DISTANCE
        rmax=0
        pmax=0
        pavg=0
        ravg=0
        davg=0
        for i in xrange(0,nlarge):
            if xs[i,NSTATE-1]==0:continue
            rp=NORM(xs[i,0:3])*UL/1E3
            pp=NORM(ps[i,0:2])*UL/1E3
            davg+=abs(ps[i,2])*UL/1E3
            """
            print "Distance of fragment %d in space: %e"%(i,rp)
            print "Distance of fragment %d on image: %e"%(i,pp)
            print
            """
            rmax=max(rmax,rp)
            pmax=max(pmax,pp)
            ravg+=rp
            pavg+=pp
        davg/=nlarge
        ravg/=nlarge
        pavg/=nlarge
        dmax=(pavg/davg)/D2R*3600
        
        print "Average distance of large fragments: %e km"%davg
        print "Average distance large fragments in image plane: %e km"%pavg
        print "Average distance large fragments in space: %e km"%ravg
        print "Maximum distance large fragments in space: %e km"%(rmax)
        print "Maximum distance large fragments in image: %e km = %e arcsec"%(pmax,dmax)
        
        fmax.write("%-07d %-+24.17e %-+12.5e %-+12.5e %-+12.5e %-+12.5e %-+12.5e %-+12.5e\n"%(iobs,(t+tini),
                                                                                              davg,ravg,pavg,
                                                                                              rmax,pmax,dmax))
    fmax.close()

#////////////////////////////////////////
#DISTANCE DISTRIBUTION
#////////////////////////////////////////
def DistanceDistribution(**args):
    """
    Calculate distribution of distance for large boulders
    """
    title("Distance Distribution")

    #for iobs in xrange(50,ncom):    
    miobs=[1,10,30]
    for iobs in miobs:
        #TIME AND DATE
        t=ts[iobs]
        date=dates[iobs]
        print "Time since integration start: t = %e days"%(t*365.25)
        print "Date time: t = %e yrs"%(t+tini)
        print "Date: %s"%date
        print 
        
        #GET DEBRIS POSITION OVER PLANE OF SKY
        ps=xs_obs[nlarge:,iobs,:]

        #FILTER ACCORDING TO RADII
        Ms=ps[:,6]*UM
        Rs=(Ms/(4*pi/3*rhoc))**(1./3)
        bould=(Rs>1)*(Rs<100)
        #bould=(Rs>30)
        xb=ps[bould]
        nbould=xb.shape[0]
        print "Number of boulders: %d"%nbould
        
        #CALCULATE DISTANCE OVER SKY PLANE OF BOULDERS
        davg=abs(xb[:,2]).mean()
        rb=array([NORM(xb[i,0:2]) for i in xrange(0,nbould)])*UL/1E3
        rb=rb/AU/D2R*3600

        #CALCULATE HISTOGRAM
        H,B=histogram(rb,bins=10,normed=False)
        
        fhis=open("boulder-distribution-%05d.dat"%iobs,"w")
        F=0
        for i in xrange(0,H.shape[0]):
            F+=H[i]
            fhis.write("%e %e\n"%(B[i],nbould-F))

        print F
        fhis.close()

#////////////////////////////////////////
#DISTANCE DISTRIBUTION
#////////////////////////////////////////
def DustDensityMap(**args):
    """
    Calculate the density of dust in the debris-zone
    """
    title("Dust density map")

    miobs=[16]
    for iobs in miobs:

        #TIME AND DATE
        t=ts[iobs]
        date=dates[iobs]
        print "Time since integration start: t = %e days"%(t*365.25)
        print "Date time: t = %e yrs"%(t+tini)
        print "Date: %s"%date
        print 

        figfile="animation/densitymap-snap_%05d.png"%iobs
        fig=plt.figure(figsize=(8,8))
        plt.close("all")
        ax=fig.add_axes([0.1,0.1,0.8,0.8],axisbg='w')

        #DISTANCE OF DEBRIS ZONE CENTER
        dist=abs(xcm_obs[iobs,2])
        print "Distance of debris zone: %e AU"%(dist)
        
        #GET POSITION OF LARGE FRAGMENTS
        xlarge=xs_obs[:nlarge,iobs,:]
        surviving=(xlarge[:,6]>0)
        xlarge=xlarge[surviving]
        nl=xlarge.shape[0]
        print "Surviving large fragments: %d"%nl

        #GET POSITION OF DEBRIS
        xdust=xs_obs[nlarge:,iobs,:]
        surviving=(xdust[:,6]>0)
        xdust=xdust[surviving]
        nd=xdust.shape[0]
        Ms=xdust[:,6]*UM
        Rs=(Ms/(4*pi/3*rhoc))**(1./3)
        print "Surviving dust particles: %d"%nd

        #GET POSITION OF BOULDERS
        bould=(Rs>1)*(Rs<100)
        xbould=xdust[bould]
        nb=xbould.shape[0]
        print "Number of boulders: %d"%nb

        #GET POSITION AND RADIUS OF PURE DUST
        rhod=2E3
        pdust=(Rs<10)
        xpdust=xdust[pdust]
        npd=xpdust.shape[0]
        md=xpdust[:,6]*UM
        rd=(md/(4*pi/3*rhod))**(1./3)
        logrd=log10(rd)
        print "Number of pure dust particles: %d"%npd

        print 

        #GET MAXIMUM FRAGMENT DISTANCE
        rlarge=array([NORM(xlarge[i,0:2]) for i in xrange(0,nl)])
        rmax=4*rlarge.max()
        print "Maximum large fragment distance: %e km\n"%(rmax*UL/1E3)

        #SCAN REGION
        Ncell=20
        dx=dy=2*rmax/Ncell
        ND=zeros((Ncell+1,Ncell+1))
        AD=zeros((Ncell+1,Ncell+1))
        sbins=arange(-6,3,1)
        nbins=sbins.shape[0]
        weights=[]
        for i in xrange(0,nbins-1):
            R1=10**sbins[i]
            R2=10**sbins[i+1]
            weights+=[1/(3+q)*(R2**(3+q)-R1**(3+q))]
        weights=array(weights)
        #print sbins,weights
        #exit(0)
            
        for i in xrange(0,Ncell+1):
            xl=-rmax+i*dx
            #print "i = %d, X = %e"%(i,xl)
            for j in xrange(0,Ncell+1):
                yl=-rmax+j*dx
                #print "\tj = %d, Y = %e"%(j,yl)
                inside=((xpdust[:,0]-xl)>0)*((xpdust[:,0]-xl)<dx)*((xpdust[:,1]-yl)>0)*((xpdust[:,1]-yl)<dy)
                ms=xpdust[inside,6]
                logrds=logrd[inside]
                try:
                    ND[i,j]+=ms.shape[0]
                    h,b=histogram(logrds,sbins)
                    AD[i,j]=(h*weights).sum()
                    #print AD[i,j]
                except:
                    pass

        AD/=AD.sum()
        """
        for i in xrange(0,Ncell+1):
            xl=-rmax+i*dx
            print "i = %d, X = %e"%(i,xl*UL/1E3)
            for j in xrange(0,Ncell+1):
                yl=-rmax+j*dy
                print "\tj = %d, Y = %e"%(j,yl*UL/1E3),
                print "%e"%AD[i,j]
        """

        AD=AD.transpose()
        imgplot=ax.imshow(AD,interpolation='bicubic',origin='lower',cmap=cm.gray)
        #imgplot=ax.imshow(AD,interpolation='nearest',origin='lower')
        #fig.colorbar(imgplot)
        ninside=ND.sum()
        print "Fraction of dust particles inside region: %e"%(ninside/npd) 
        
        comet(ax,
              map(xlarge[:,0],rmax,Ncell),
              map(xlarge[:,1],rmax,Ncell),
              zorder=10)
        """
        ax.plot(map(xlarge[:,0],rmax,Ncell),
                map(xlarge[:,1],rmax,Ncell),
                'ow',zorder=10,markersize=3,markeredgecolor='none')
        """

        ax.plot(map(xbould[:,0],rmax,Ncell),
                map(xbould[:,1],rmax,Ncell),
                'o',color='w',zorder=10,markersize=1,markeredgecolor='none')

        comet(ax,
              map(rmax/2,rmax,Ncell),
              map(rmax/2,rmax,Ncell),
              zorder=10)

        #**************************************************
        #AXIS
        #**************************************************
        zdir=zaxis[iobs,0:3]
        zdir=5*rmax*zdir
        ax.plot(map([0.0,zdir[0]],rmax,Ncell),
                map([0.0,zdir[1]],rmax,Ncell),
                'w-',linewidth=2,linestyle='-')
        ax.text(map(zdir[0],rmax,Ncell),
                map(zdir[1],rmax,Ncell),
                "A",color='w',fontsize=10,
                bbox=dict(boxstyle='square,pad=0.5',fc="0.0",ec="none"))
        
        zsun=dsun[iobs,0:3]
        zsun=5*rmax*zsun
        ax.plot(map([0.0,zsun[0]],rmax,Ncell),
                map([0.0,zsun[1]],rmax,Ncell),
                'w-',linewidth=2,linestyle='-')
        ax.text(map(zsun[0],rmax,Ncell),
                map(zsun[1],rmax,Ncell),
                "S",color='w',fontsize=10,
                bbox=dict(boxstyle='square,pad=0.5',fc="0.0",ec="none"))


        #**************************************************
        #DECORATION
        #**************************************************
        ax.set_xlim(-0.5,Ncell+0.5)
        ax.set_ylim(-0.5,Ncell+0.5)

        ax.set_xticks(concatenate((linspace(0,Ncell/2,5),linspace(Ncell/2,Ncell,5))))
        ax.set_yticks(concatenate((linspace(0,Ncell/2,5),linspace(Ncell/2,Ncell,5))))
        
        xt=ax.get_xticks()
        xl=[]
        for x in xt:
            xp=-rmax+x*dx
            xl+=["%.1f"%(xp*UL/1E3)]
        ax.set_xticklabels(xl)
        ax.set_yticklabels(xl)

        ax.axhline(map(0,rmax,Ncell),color='w',linewidth=0.5,linestyle='--')
        ax.axvline(map(0,rmax,Ncell),color='w',linewidth=0.5,linestyle='--')
        #**************************************************
        #SAVE FIG
        #**************************************************
        fig.savefig(figfile);
        

        
#############################################################
#SELECT TASK
#############################################################
DustDensityMap()
exit(0)

#DistanceDistribution()
#exit(0)

#MaximumDistance()
#exit(0)

"""
Snapshot(iobs=1,facvel=1)
Observation(iobs=1)
exit(0)
Orbit(iobs=1)
exit(0)
Observation(iobs=1)
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
