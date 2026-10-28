Path="/home/ywitzky/HPS/HOOMD/Test/Calvados2_CPU"

from Submit_HOOMD import run
import os
import numpy as np
from multiprocessing import Pool
import time

Paths= []
for name in os.listdir(Path):
    if os.path.isdir(os.path.join(Path, name)+"/RUN_001/"):
        Paths.append(os.path.join(Path, name)+"/RUN_001/")

#Paths = Paths[0:2]
NJobs = len(Paths)
print(NJobs)

#with Pool(NJobs) as p:
#    p.map(run, Paths)
#    print("asdfghj")
    
if True:
    p = Pool(NJobs)

    results = []
    for x in range(NJobs):
        results.append( p.apply_async(run, [Paths[x]]))

    p.close()

    for r in results:
        r.get()