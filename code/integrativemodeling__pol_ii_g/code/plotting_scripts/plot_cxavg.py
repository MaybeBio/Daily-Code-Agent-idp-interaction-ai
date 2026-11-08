import matplotlib
import matplotlib.pyplot as plt
from matplotlib import pyplot;
import os,sys,string

savename=sys.argv[1]
filename=sys.argv[2]
sf=open(filename,'r')
mu=[];sig=[];name=[]
for i,ln in enumerate(sf.readlines()):
    line =ln.strip().split()
    mu.append(float(line[5]))
    sig.append(float(line[6]))
    name.append(str(line[0])+"-"+str(line[1])+" "+str(line[2])+"-"+str(line[3]))
XL=range(1,len(name)+1)
tot=[35]*len(name)

fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
fig.set_size_inches(13, 9.5)
plt.plot(XL,tot,'--',color='r',linewidth=3.0)
plt.ylim([0,36])
plt.xlim([0,42])
plt.errorbar(XL,mu,sig,fmt='o')
XLc=range(1,len(name)+1)
plt.xticks(XLc, name, rotation=90)
ax.xaxis.labelpad=20
ax.xaxis.set_tick_params(pad=10)
ax.yaxis.set_tick_params(pad=10)
plt.subplots_adjust(bottom=0.25)
# Set the tick labels font
for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(11)
plt.ylabel('Distance')
plt.savefig(savename +"XL" + ".png", format='png',dpi=100)
plt.hold(True)

