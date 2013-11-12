from matplotlib import pyplot as plt,patches as pat,cm,image as img
from numpy import *
from scipy import interpolate,optimize
from sys import *
from os import system

#############################################################
#CONSTANTS
#############################################################
INTERP=interpolate.interp1d
BISECT=optimize.bisect
NORM=linalg.norm
NSTATE=7
D2R=pi/180

#############################################################
#UTIL ROUTINES
#############################################################
def title(text):
    print "*"*80,"\n",text,"\n","*"*80,"\n"

def angsize(r,d):
    t=(r/d)/D2R*3600
    return t

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
FAC=UL/1E3 #km
AU=1.496E8 #km

#//////////////////////////////////////////////////
#DATES
#//////////////////////////////////////////////////
fdate=open("dates.dat","r")
dates=[]
for line in fdate:
    tphys,date=line.split("=")
    dates+=[date.strip()]
fdate.close()

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
#SPECIAL ROUTINES
#############################################################

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
#DISTANCE DISTRIBUTION
#////////////////////////////////////////
RBINS=arange(-6,6,1)
NRBINS=RBINS.shape[0]
WEIGHTS=[]
for i in xrange(0,NRBINS-1):
    R1=10**RBINS[i]
    R2=10**RBINS[i+1]
    WEIGHTS+=[1/(3+q)*(R2**(3+q)-R1**(3+q))]
WEIGHTS+=[0]
WEIGHTS=array(WEIGHTS)

def radialDistribution(rdust,Rd,Rmin=1E-3,Rmax=2E3,nbin=20,q=0.2,
                       weights=WEIGHTS):

    #//////////////////////////////////////////////////
    #FILTER
    #//////////////////////////////////////////////////
    select=(Rd>Rmin)*(Rd<Rmax)
    Rd=Rd[select]
    rdust=rdust[select]
    logrdust=log10(rdust)

    #//////////////////////////////////////////////////
    #BINNING
    #//////////////////////////////////////////////////
    #RADIAL BINNING
    logrmin=logrdust.min()
    logrmax=logrdust.max()
    dlogr=(logrmax-logrmin)/nbin

    logrd=[]
    A=[]
    Ac=[]
    F=0

    for i in xrange(0,nbin):
        logr1=logrmin+i*dlogr
        logr2=logr1+dlogr

        inbinr=(logrdust>logr1)*(logrdust<logr2)
        logrbin=log10(Rd[inbinr])
        
        h,b=histogram(logrbin,RBINS)

        logrd+=[logr1]
        Abin=(h*weights[:-1]).sum()
        A+=[1.0*Abin]
        F+=Abin
        Ac+=[F]

    #//////////////////////////////////////////////////
    #NORMALIZATION
    #//////////////////////////////////////////////////
    logrd=array(logrd)
    A=array(A)/(1.0*F)

    #//////////////////////////////////////////////////
    #REPORT
    #//////////////////////////////////////////////////
    rd=10**logrd
    Ac=1-array(Ac)/(1.0*F)

    #//////////////////////////////////////////////////
    #LIMIT
    #//////////////////////////////////////////////////
    Acint=INTERP(rd,Ac)
    cross=lambda r:Acint(r)-q
    rq=BISECT(cross,rd[0],rd[-1])

    return rd,Ac,rq

def AreaDistribution(iobs=1,**args):
    """
    Calculate the distribution of area as a function of angular
    distance to debris zone center
    """
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #TIME AND DATE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t=ts[iobs]
    date=dates[iobs]
    print "Snapshot: %d"%iobs
    print "Time since integration start: t = %e days"%(t*365.25)
    print "Date time: t = %e yrs"%(t+tini)
    print "Date calendar: %s"%date

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #PREPARE FIGURE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figfile="animation/areadistrib-snap%05d.png"%iobs
    fig=plt.figure(figsize=(8,8))
    plt.close("all")
    ax=fig.add_axes([0.1,0.1,0.8,0.8],axisbg='w')

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #PREPARE DATA
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #DISTANCE OF DEBRIS ZONE CENTER
    dist=abs(xcm_obs[iobs,2])
    rdist=NORM(xcm_orb[iobs,0:3])
    print "Distance of debris zone: %e AU"%(dist)
    print "Distance heliocentric: %e AU"%(rdist)
    
    #//////////////////////////////////////////////////
    #GET POSITION OF LARGE FRAGMENTS
    #//////////////////////////////////////////////////
    xlarge=xs_obs[:nlarge,iobs,:]
    surviving=(xlarge[:,6]>0)
    xlarge=xlarge[surviving]
    nl=xlarge.shape[0]
    rlarge=array([NORM(xlarge[i,0:2]) for i in xrange(0,nl)])
    rmax=rlarge.max()
    print "Surviving large fragments: %d"%nl
    print "Maximum distance of large particles: %e"%(rmax*FAC)

    #//////////////////////////////////////////////////
    #GET POSITION OF BOULDERS
    #//////////////////////////////////////////////////
    xbould=xs_obs[nlarge:,iobs,:]
    Mb=xbould[:,6]*UM
    Rb=(Mb/(4*pi/3*rhoc))**(1./3)
    surviving=(Rb>10)
    xbould=xbould[surviving]
    Rb=Rb[surviving]
    nb=xbould.shape[0]
    rbould=array([NORM(xbould[i,0:2]) for i in xrange(0,nb)])
    rbmax=rbould.max()
    print "Surviving bould particles: %d"%nb
    print "Maximum distance of bould particles: %e"%(rbmax*FAC)

    #//////////////////////////////////////////////////
    #GET POSITION OF DEBRIS
    #//////////////////////////////////////////////////
    xdust=xs_obs[nlarge:,iobs,:]
    surviving=(xdust[:,6]>0)
    xdust=xdust[surviving]
    nd=xdust.shape[0]
    rdust=array([NORM(xdust[i,0:2]) for i in xrange(0,nd)])
    rdmax=rdust.max()
    print "Surviving dust particles: %d"%nd
    print "Maximum distance of dust particles: %e"%(rdmax*FAC)
    Md=xdust[:,6]*UM
    Rd=(Md/(4*pi/3*rhoc))**(1./3)

    rd,Ac,rq=radialDistribution(rdust,Rd,nbin=50,q=0.1,Rmin=1E-3,Rmax=10)

    dmax=angsize(rmax,dist)
    print "dmax = %e"%(dmax)

    dq=angsize(rq,dist)
    print "rq = %e"%(dq)

    rb,Ab,rqb=radialDistribution(rbould,Rb,nbin=30,q=0.1,Rmin=10,Rmax=500,
                                 weights=ones_like(RBINS))

    dbmax=angsize(rbmax,dist)
    print "dbmax = %e"%(dbmax)

    dqb=angsize(rqb,dist)
    print "rqb = %e"%(dqb)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #PLOT
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax.plot(angsize(rd,dist),100*Ac,color='b',label=r'Dust ($1\,{\rm mm}<R_i<10\,{\rm m}$')
    ax.plot(angsize(rb,dist),100*Ab,color='r',label=r'Boulders ($R_i>10\,{\rm m}$)')

    random.seed(1)

    height=15+5*random.rand(rlarge.shape[0])
    ax.plot(angsize(rlarge,dist),height,'ok')
    ax.text(dmax/20,22,"Large fragments",fontsize=10,
            horizontalalignment='left',verticalalignment='bottom')

    height=5+5*random.rand(rbould.shape[0])
    ax.plot(angsize(rbould,dist),height,'ok',markersize=2)
    ax.text(dmax/20,12,"Boulders",fontsize=10,
            horizontalalignment='left',verticalalignment='bottom')

    ax.axvspan(0,dmax,color=cm.gray(0.1),alpha=0.2)
    ax.axvline(dmax,color='k',linewidth=1)
    ax.text(dmax,102,r'$d_{\rm max}$',horizontalalignment='left')

    ax.axvspan(0,dq,color=cm.gray(0.5),alpha=0.2)
    ax.axvline(dq,color='b',linewidth=1)
    ax.text(dq,102,r'$d_{90,D}$',horizontalalignment='left')

    ax.axvspan(0,dqb,color=cm.gray(0.8),alpha=0.2)
    ax.axvline(dqb,color='r',linewidth=1)
    ax.text(dqb,102,r'$d_{90,B}$',horizontalalignment='left')
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #DECORATION
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax.legend(loc='upper right')
    ax.set_xlim((0,1.5*angsize(rqb,dist)))

    ax.set_xlabel('$d$ (arcsec)')
    ax.set_ylabel('$F_A(d)$')

    date=" ".join(date.split()[0:3])
    ax.set_title("%s: $r$ = %.2f AU, $\Delta$ = %.2f AU"%(date,rdist,dist),
                 position=(0.5,1.05))

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #SAVE FIGURE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print "Saving %s..."%figfile
    fig.savefig(figfile);

def ExpansionEvolution(did='0000',imax=ncom,**args):
    """
    Calculate the distribution of area as a function of angular
    distance to debris zone center
    """

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #PREPARE FIGURE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dateini="_".join(dates[0].split()[0:3])
    figfile="animation/expansion-%s.png"%dateini
    print "Creating: %s"%figfile
    fig=plt.figure(figsize=(8,8))
    plt.close("all")
    ax=fig.add_axes([0.1,0.1,0.8,0.8],axisbg='w')
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #TIME AND DATE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #CALCULATE PERIHELION AND PERIGEO
    dperh=1E100
    dperg=1E100
    tperg=0
    tperh=0
    for iobs in xrange(0,ncom):
        t=ts[iobs]
        dist=abs(xcm_obs[iobs,2])
        rdist=NORM(xcm_orb[iobs,0:3])
        if dist<=dperg:
            tperg=t
            dperg=dist
        if rdist<=dperh:
            tperh=t
            dperh=rdist

    dmaxs=[]
    dqs=[]
    dqbs=[]
    tes=[]
    for iobs in xrange(0,imax):
        t=ts[iobs]
        tes+=[t+tini]
        date=dates[iobs]
        print "Snapshot: %d"%iobs

        dist=abs(xcm_obs[iobs,2])
        rdist=NORM(xcm_orb[iobs,0:3])
    
        xlarge=xs_obs[:nlarge,iobs,:]
        surviving=(xlarge[:,6]>0)
        xlarge=xlarge[surviving]
        nl=xlarge.shape[0]
        rlarge=array([NORM(xlarge[i,0:2]) for i in xrange(0,nl)])
        rmax=rlarge.max()
        dmax=angsize(rmax,dist)
        dmaxs+=[dmax]

        xbould=xs_obs[nlarge:,iobs,:]
        Mb=xbould[:,6]*UM
        Rb=(Mb/(4*pi/3*rhoc))**(1./3)
        surviving=(Rb>10)
        xbould=xbould[surviving]
        Rb=Rb[surviving]
        nb=xbould.shape[0]
        rbould=array([NORM(xbould[i,0:2]) for i in xrange(0,nb)])
        rbmax=rbould.max()
        dbmax=angsize(rbmax,dist)

        xdust=xs_obs[nlarge:,iobs,:]
        surviving=(xdust[:,6]>0)
        xdust=xdust[surviving]
        nd=xdust.shape[0]
        rdust=array([NORM(xdust[i,0:2]) for i in xrange(0,nd)])
        rdmax=rdust.max()
        Md=xdust[:,6]*UM
        Rd=(Md/(4*pi/3*rhoc))**(1./3)

        rd,Ac,rq=radialDistribution(rdust,Rd,nbin=50,q=0.1,Rmin=1E-3,Rmax=10)
        dq=angsize(rq,dist)
        dqs+=[dq]
        
        rb,Ab,rqb=radialDistribution(rbould,Rb,nbin=30,q=0.1,Rmin=10,Rmax=500,
                                     weights=ones_like(RBINS))
        dqb=angsize(rqb,dist)
        dqbs+=[dqb]

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #SAVE FILES
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data=zeros((len(tes),4))
    data[:,0]=array(tes)
    data[:,1]=array(dmaxs)
    data[:,2]=array(dqs)
    data[:,3]=array(dqbs)
    savetxt("output/expansion-%s-%s.dat"%(dateini,did),data)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #PLOT
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax.plot(tes,dmaxs,label='Large')
    ax.plot(tes,dqs,label='Dust')
    ax.plot(tes,dqbs,label='Boulders')
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #DECORATION
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ndates=len(dates)
    pdates=ndates/10
    xt=[]
    xtl=[]
    for i in xrange(0,ndates,pdates):
        xt+=[ts[i]+tini]
        d=dates[i]
        xtl+=["-".join(d.split()[1:3])]
    ax.set_xticks(xt)
    ax.set_xticklabels(xtl,rotation=-20,
                       rotation_mode='anchor',
                       horizontalalignment='left',
                       fontsize=12)

    ax.set_yscale("log")
    ax.legend(loc='lower right')

    ax.set_ylabel("$d$ (arcsec)")

    dMin,dMax=ax.get_ylim()
    ax.axvline(tini,linestyle='--',color='k')
    ax.text(tini,1.1*dMax,"Disruption",horizontalalignment='left')

    ax.axvline(tperh+tini,linestyle='--',color='k')
    ax.text(tperh+tini,1.1*dMax,"Perihelion",horizontalalignment='left')

    ax.axvline(tperg+tini,linestyle='--',color='k')
    ax.text(tperg+tini,1.1*dMax,"Perigeo",horizontalalignment='right')

    ax.set_ylim((10E-3,dMax))

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #SAVE FIGURE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print "Saving %s..."%figfile
    fig.savefig(figfile);

#############################################################
#SELECT TASK
#############################################################
#ExpansionEvolution()
#AreaDistribution()
#exit(0)

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
