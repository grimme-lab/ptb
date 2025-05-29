import os
import sys

thres=["0.0001","0.00001","0.000001"]
modes=["full_matrix","kmeans_r_heuristic"]
systems=[]
charges=[]
dims=[]

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

for m in modes:
    for i in range(len(systems)):
        for f in thres:
            tsolve2=0
            tsm=0
            dim=0
            nsm=0
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
                    if l.startswith("PTB H matrix iteration                      2"):
                        break
                fi.close()
                print(f,systems[i],m,dim,nsm,tsolve2,tsm)
            except:
                pass



