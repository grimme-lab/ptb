import os
import math
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator) 

thres=["0.0001","0.00001","0.000001"]
modes=["full_matrix","kmeans_r_heuristic"]
systems=[]
charges=[]
dims=[]

sygvd={}

for s in os.listdir("systems/"):
    if s.endswith(".xyz"):
        s2=s.replace(".xyz","")
        systems.append(s2)
        f = open("systems/"+s2+".charge")
        charge=int(f.read())
        f.close()
        charges.append(charge)
print(systems)
print(charges)

data={}

for m in modes:
    for f in thres:
        for i in range(len(systems)):
            tsolve2=0
            tsm=0
            dim=0
            nsm=0
            de=0
            mse=0
            try:
                fi=open("runs/"+systems[i]+"_"+str(charges[i])+"/filter_"+f+"/"+m+"/out")
                for l in fi.read().splitlines():
                    if l.startswith("Total solve2"):
                        tsolve2=float(l.split()[2])
                    if l.startswith("Submatrix "):
                        tsm=float(l.split()[1])
                    if l.startswith(" basis setup done. Ndim"):
                        dim=int(l.split()[4])
                    if l.startswith(" number of submatrices="):
                        nsm=int(l.split()[3])
                    if l.startswith("* sygvd solver"):
                        sygvd[dim]=float(l.split()[3])
                    if l.startswith("PTB H matrix iteration                      2"):
                        break
                    if l.startswith(" energy deviation per atom"):
                        de=abs(float(l.split()[7].replace("E","e")))
                    if l.startswith(" MSE of density matrix"):
                        mse=abs(float(l.split()[4].replace("E","e")))
                fi.close()
                data[m+"_"+str(f)+"_"+str(systems[i])+"_nsm"]=nsm
                data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"]=dim
                data[m+"_"+str(f)+"_"+str(systems[i])+"_tsolve2"]=tsolve2
                data[m+"_"+str(f)+"_"+str(systems[i])+"_tsm"]=tsm
                data[m+"_"+str(f)+"_"+str(systems[i])+"_de"]=de
                data[m+"_"+str(f)+"_"+str(systems[i])+"_mse"]=mse
            except:
                pass
print(data)
#from scipy.optimize import curve_fit
#import numpy as np
#def func(x, a, b):
#    return a * x**b
#
#x=np.zeros(len(sygvd))
#y=np.zeros(len(sygvd))
#
#x=sorted(sygvd)
#i=0
#for xi in x:
#    y[i]=sygvd[xi]
#    i=i+1
#print(x,y)
#
#popt, pcov = curve_fit(func, x[8:], y[8:])
#print(popt)

#plot
col=[]
for i in mcolors.TABLEAU_COLORS:
    col.append(i)
ph=4
pw=ph*0.6666
bb=0.02

fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph, pw))
ax.set_xlabel(r'Number of Basis Functions')
ax.set_ylabel(r'Speedup')
ax.set_xscale("log")
#ax.set_yscale("log")

c=0
f=0.0001
m="full_matrix"
x0=[]
y0=[]
for i in range(len(systems)):
    x0.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
    y0.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_tsolve2"])
y0 = sorted(y0, key=lambda z: x0[y0.index(z)])
x0=sorted(x0)
#ax.plot(x0,y0,linewidth=1,label="truncation "+str(f),color=col[c])
#c=c+1

for f in thres:
    m="kmeans_r_heuristic"
    x=[]
    y=[]
    for i in range(len(systems)):
        x.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
        y.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_tsm"])
    y = sorted(y, key=lambda z: x[y.index(z)])
    x=sorted(x)
    for i in range(len(systems)):
        y[i]=y0[i]/y[i]
    print(x)
    print(y)
    ax.plot(x,y,linewidth=1,label=r'truncation $10^{'+str(int(math.log10(float(f))))+'}$',color=col[c],marker=".",markersize=4)
    c=c+1


plt.legend(loc='upper left',ncol=1,fontsize="7")
plt.tight_layout()
fig.savefig("submatrix_speedup.pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_speedup.png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)

#######################################################################################
fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph, pw))
ax.set_xlabel(r'Number of Basis Functions')
ax.set_ylabel(r'Number of Submatrices')
ax.set_xscale("log")
ax.set_yscale("log")
c=0

f="0.0001"
m="full_matrix"

for f in thres:
    m="kmeans_r_heuristic"
    x=[]
    y=[]
    for i in range(len(systems)):
        x.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
        y.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_nsm"])
    y = sorted(y, key=lambda z: x[y.index(z)])
    x=sorted(x)
    print(x)
    print(y)
    ax.plot(x,y,linewidth=1,label=r'truncation $10^{'+str(int(math.log10(float(f))))+'}$',color=col[c],marker=".",markersize=4)
    c=c+1


plt.legend(loc='upper left',ncol=1,fontsize="7")
plt.tight_layout()
fig.savefig("submatrix_nsm.pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_nsm.png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)


#######################################################################################
fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph,pw))
ax.set_xlabel(r'Number of Basis Functions')
ax.set_ylabel(r'$|\Delta E_\mathrm{tr}$| / a.u.')
ax.set_xscale("log")
ax.set_yscale("log")
c=0

f="0.0001"
m="full_matrix"
x0=[]
y0=[]
for i in range(len(systems)):
    x0.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
    y0.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_tsolve2"])
y0 = sorted(y0, key=lambda z: x0[y0.index(z)])
x0=sorted(x0)

for f in thres:
    m="kmeans_r_heuristic"
    x=[]
    y=[]
    for i in range(len(systems)):
        x.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
        y.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_de"])
    y = sorted(y, key=lambda z: x[y.index(z)])
    x=sorted(x)
    print(x)
    print(y)
    ax.plot(x,y,linewidth=1,label=r'truncation $10^{'+str(int(math.log10(float(f))))+'}$',color=col[c],marker=".",markersize=4)
    c=c+1


plt.legend(loc='upper left',ncol=1,fontsize="7")
plt.tight_layout()
fig.savefig("submatrix_error.pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_error.png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)

#######################################################################################
fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph, pw))
ax.set_xlabel(r'Number of Basis Functions')
ax.set_ylabel(r'MSE of Density Matrix')
ax.set_xscale("log")
ax.set_yscale("log")
c=0

f="0.0001"
m="full_matrix"
x0=[]
y0=[]
for i in range(len(systems)):
    x0.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
    y0.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_tsolve2"])
y0 = sorted(y0, key=lambda z: x0[y0.index(z)])
x0=sorted(x0)

for f in thres:
    m="kmeans_r_heuristic"
    x=[]
    y=[]
    for i in range(len(systems)):
        x.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
        y.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_mse"])
    y = sorted(y, key=lambda z: x[y.index(z)])
    x=sorted(x)
    print(x)
    print(y)
    ax.plot(x,y,linewidth=1,label=r'truncation $10^{'+str(int(math.log10(float(f))))+'}$',color=col[c],marker=".",markersize=4)
    c=c+1


plt.legend(loc='upper left',ncol=1,fontsize="7")
plt.tight_layout()
fig.savefig("submatrix_errorP.pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_errorP.png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)

