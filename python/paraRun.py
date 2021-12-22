#%%
import os
import shutil
import time

import numpy as np
import sobol_seq

#%%
t = [25e-6, 50e-6]
R = [0.02, 0.04]
# c = [0, 1]

data = np.vstack([t, R])

size = 400
vec = sobol_seq.i4_sobol_generate(data.shape[0], size)

lb = data[:, 0]
ub = data[:, 1]

out = vec * (ub - lb) + lb
# print(out)

#%%
# source_dir = r"D:\projects\CYLINDER_Du"
root = r"C:\Users\edison\workspace\matlab"
source_dir = root + "\CYLINDER_Du"

nproc = 10

# %%
for i, geo in enumerate(np.array_split(out, nproc)):
    # print(i, geo)
    destination_dir = root + rf"\case_{i}"
    shutil.copytree(source_dir, destination_dir, ignore=shutil.ignore_patterns(".git"))

#%%
for i, geo in enumerate(np.array_split(out, nproc)):
    destination_dir = root + rf".\case_{i}"
    np.savetxt(destination_dir + "/geo.csv", geo, delimiter=",")

# %%
for i in range(nproc):
    # f"cd d:\projects\case_{i}\ && matlab -nosplash -nodesktop -r GENERATE_CDB"
    os.system(rf"cd {root}\case_{i}\ && matlab -nosplash -nodesktop -r GENERATE_CDB")
    time.sleep(45)

# %%
for i in range(nproc):
    # f"cd d:\projects\case_{i}\ && matlab -nosplash -nodesktop -r NONLINEAR_SOLVE"
    os.system(rf"cd {root}\case_{i}\ && matlab -nosplash -nodesktop -r NONLINEAR_SOLVE")
    time.sleep(1)

#%%
for i in range(nproc):
    shutil.rmtree(root + rf"\case_{i}")
