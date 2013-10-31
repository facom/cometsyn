"""
COMET DESINTEGRATION
SNAPSHOT
"""
from matplotlib import pyplot as plt,patches as pat,cm
from mpl_toolkits.mplot3d import Axes3D
from numpy import *
from sys import *
from os import system

#############################################################
#CONSTANTS
#############################################################
NORM=linalg.norm

#############################################################
#OBSERVATIONAL DATA
#############################################################
orbdata=loadtxt("comet-orbit.dat");
config={}
execfile("fragments.gph",config);

norb=orbdata.shape[0]
nfrag=config['nfrag']
nlarge=config['nlarge']
ndebris=config['ndebris']
tini=config['tini']

print "Orbital properties:"
print "\tNumber of points:%d"%norb
print "\tNumber of fragments:%d"%nfrag
print "\tLarge fragments:%d"%nlarge
print "\tDebris:%d"%ndebris

ts=orbdata[:,0]
xcm=orbdata[:,1:7]
xs=[]
print orbdata.shape
for i in xrange(0,nfrag):
    k=6*(i+1)+1
    xs+=[orbdata[:,k:k+6]]
xs=array(xs)

#############################################################
#SELECT SNAPSHOT
#############################################################
try:
    iobs=int(argv[1])
except:
    iobs=0
if iobs>=norb:
    print "Maximum snapshot %d (t = %e)"%(norb,ts[norb-1]+tini)
    exit(1)

t=ts[iobs]
print "Snaptshot: %d, t - tper = %.2f yrs"%(iobs,t+tini)

rcm=xcm[iobs,0:6]
d=NORM(rcm[0:3])
D=d
f=open("comet-fragments-snapshots.dat","w")
rmax=0

for i in xrange(0,nfrag):
    rs=xs[i,iobs,0:6]
    rf=rs-rcm
    rfnorm=NORM(rf[0:3])
    rmax=max(rmax,rfnorm)
    type=1.0
    if i>=nlarge:type=2.0
    f.write("%e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e %-+23.17e\n"%(type,
                                                                                rf[0],
                                                                                rf[1],
                                                                                rf[2],
                                                                                rf[3],
                                                                                rf[4],
                                                                                rf[5]))
f.close()
rmax=rmax*config['UL']
print "Maximum distance: %e"%rmax

#############################################################
#SAVE fragments.gph
#############################################################
f=open("fragments.gph","w")
f.write("""\
file='comet-fragments-snapshots.dat'
title='t = %.2f yrs, d = %.2f AU, D = %.2f AU'
nfrag=%d
nlarge=%d
ndebris=%d
tini=%e
UM=%e
UL=%e
UT=%e
Rmax=%e
Rc=%e
rf=%e
"""%(t+tini,d,D,nfrag,nlarge,ndebris,tini,
     config['UM'],config['UL'],config['UT'],
     rmax,config['Rc']*10,config['rf']))
f.close()

#############################################################
#PLOT
#############################################################
system("gnuplot plot-fragments.gpl")
