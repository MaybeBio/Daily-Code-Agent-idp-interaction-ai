import matplotlib
import matplotlib.pyplot as plt
from matplotlib import pyplot;
import os,sys,string

savename=sys.argv[1]
filename=sys.argv[2]
sf=open(filename,'r')
dist=[]
for i,ln in enumerate(sf.readlines()):
    line =ln.strip().split()
    dist.append(float(line[0]))
y,bine=np.histogram(dist,normed=True,bins=list(range(0,100,1)))
binc=0.5*(bine[1:] +bine[:-1])
x_cutoff=[35]*100
y_cutoff=np.linspace(-1,0.2,100)
fig=plt.figure(1)
plt.bar(binc,y) # ,color='b',linewidth=3.0)
plt.hold(True)
plt.plot(x_cutoff,y_cutoff,'r--',linewidth=3.0)
plt.ylim([0,0.16])
plt.xlim([0.,40])
ax = fig.add_subplot(111)
plt.setp(ax.spines.values(), linewidth=5)
plt.savefig(savename +"XLdisthist" + ".png",format='png',dpi=600)
plt.hold(True)

