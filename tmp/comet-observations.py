"""
COMET DESINTEGRATION
OBSERVATIONS SCRIPTS
"""
from matplotlib import pyplot as plt,patches as pat,cm
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
from sys import *
from os import system

#############################################################
#CONSTANTS
#############################################################
NSTATE=6
R2S=180/pi*3600
NORM=linalg.norm
system("rm -rf animation/*.png")

#############################################################
#OBSERVATIONAL DATA
#############################################################
obsdata=loadtxt("comet-observations.dat");
nobs=obsdata.shape[0]
nsys=obsdata.shape[1]
nfrag=nsys/NSTATE-1

print "Observations properties:"
print "\tNumber of observations:%d"%nobs
print "\tNumber of fragments:%d"%nfrag

ts=obsdata[:,0]
xcm=obsdata[:,1:7]
xs=[]
for i in xrange(0,nfrag):
    k=NSTATE*(i+1)+1
    xs+=[obsdata[:,k:k+NSTATE]]
xs=array(xs)

#############################################################
#ORBITAL DATA
#############################################################
orbdata=loadtxt("comet-orbit.dat")
xorbs=orbdata[:,1:7]

earthdata=loadtxt("comet-earth.dat")
xearth=earthdata[:,1:7]

#############################################################
#SELECT SNAPSHOT
#############################################################
config={}
execfile("fragments.gph",config);

nfrag=config['nfrag']
nlarge=config['nlarge']
ndebris=config['ndebris']
tini=config['tini']

print "Orbital properties:"
print "\tNumber of points:%d"%nobs
print "\tNumber of fragments:%d"%nfrag
print "\tLarge fragments:%d"%nlarge
print "\tDebris:%d"%ndebris

try:
    iobs=int(argv[1])
except:
    iobs=0
if iobs>=nobs:
    print "Maximum snapshot %d (t = %e)"%(norb,ts[nobs-1]+tini)
    exit(1)

t=ts[iobs]
print "Snaptshot: %d, t - tper = %.2e yrs"%(iobs,t+tini)
print "Integration time: %.2e yrs = %.2e days"%(t,t*365.25)
exit(0)

#############################################################
#DIRECTIONS OBSERVED
#############################################################
dirdata=loadtxt("comet-dirs.dat")
zaxis=dirdata[:,1:7]
dsun=dirdata[:,7:13]

#############################################################
#PLOT OBSERVATIONS
#############################################################
obsmax=20
for iobs in xrange(0,nobs):
#for iobs in xrange(0,1):
    print "Plotting for iobs = %d..."%iobs
    fig=plt.figure(figsize=(8,8))
    plt.close("all")
    #ax=fig.add_axes([0.1,0.1,0.8,0.8],axisbg='k')
    ax=fig.add_axes([0,0,1,1],axisbg='k')
    
    xo=xcm[iobs,0]
    yo=xcm[iobs,1]
    zo=xcm[iobs,2]

    ax.plot([xo/abs(zo)*R2S],[yo/abs(zo)*R2S],'cx',markersize=20)
    alphamax=0
    deltamax=0
    for i in xrange(0,nfrag):
        x=xs[i,iobs,0]
        y=xs[i,iobs,1]
        z=xs[i,iobs,2]
        alpha=(x/abs(z))*R2S
        delta=(y/abs(z))*R2S
        ax.plot([alpha],[delta],'o',markersize=2,
                markeredgecolor='none',color='w')

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #AXIS
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    zdir=zaxis[iobs,0:3]
    zdir=obsmax*zdir/NORM(zdir)/2
    ax.plot([0.0,zdir[0]],[0.0,zdir[1]],'w-',linewidth=1,linestyle=':')
    ax.text(zdir[0],zdir[1],"Axis",color='w',fontsize=10)
    
    zsun=dsun[iobs,0:3]
    zsun=obsmax*zsun/NORM(zsun)/2
    ax.plot([0.0,zsun[0]],[0.0,zsun[1]],'y-',linewidth=1,linestyle=':')
    ax.text(zsun[0],zsun[1],"Sunward",color='w',fontsize=10)

    Delta=NORM(xearth[iobs,0:3]-xcm[iobs,0:3])

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #DECORATION
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    angle=2*obsmax*0.1
    ax.axhline(-0.9*obsmax,xmin=0.1,xmax=0.2,color='w')
    ax.text(-(1-2*0.1)*obsmax,-0.88*obsmax,"%.2f arcsec"%angle,color='w')
    """
    ax.plot([0],[0],'k-',label='Rotational Axis')
    ax.plot([0],[0],'r-',label='Sun direction')
    ax.legend(loc='best')
    """
    ax.set_xlabel(r"$\Delta\alpha$ (arcsec)")
    ax.set_ylabel(r"$\Delta\delta$ (arcsec)")
    ax.set_xlim((-obsmax,obsmax))
    ax.set_ylim((-obsmax,obsmax))
    ax.set_title("t = %.2f yrs, $\Delta$ = %.2f AU"%(ts[iobs],Delta/1.496E11),position=(0.5,0.95),color='w')
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #SAVE
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig.savefig("animation/plot-t_%.2f.png"%ts[iobs]);

print "Creating animation..."
system("convert -delay 20 animation/*.png -loop 0 animation/disintegration.gif")
