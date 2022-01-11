#%%
import datetime
import os
import shutil
import time

import numpy as np
import sobol_seq

#%%
size = 32
dict = {
    "t": [75e-6, 75e-6],
    "R": [0.02, 0.04],
    "p": [1e4, 5e4],
    "L": [0.5, 0.7],
}

params = ["t", "R", "p", "L"]
data = np.vstack([dict[i] for i in params])

vec = sobol_seq.i4_sobol_generate(data.shape[0], size)

lb = data[:, 0]
ub = data[:, 1]

out = vec * (ub - lb) + lb

np.random.shuffle(out)
print(out[:5, :])

#%%
root = r"C:\Users\edison\workspace\matlab"
source_dir = root + "\CYLINDER_Du"

nproc = 16

# %%
for i, geo in enumerate(np.array_split(out, nproc)):
    # print(i, geo)
    destination_dir = root + rf"\case_{i}"
    shutil.copytree(
        source_dir, destination_dir, ignore=shutil.ignore_patterns(".git", "*.h5")
    )

#%%
for i, geo in enumerate(np.array_split(out, nproc)):
    destination_dir = root + rf".\case_{i}"
    np.savetxt(destination_dir + "/geo.csv", geo, delimiter=",")

# %%
for i in range(nproc):
    print(f"{i+1}/{nproc}")
    # f"cd d:\projects\case_{i}\ && matlab -nosplash -nodesktop -r GENERATE_CDB"
    os.system(rf"cd {root}\case_{i}\ && matlab -nosplash -nodesktop -r GENERATE_CDB")
    while not (os.path.isfile(root + rf"/case_{i}/cdb_creation_finished.csv")):
        time.sleep(1)

# %%
for i in range(nproc):
    # f"cd d:\projects\case_{i}\ && matlab -nosplash -nodesktop -r NONLINEAR_SOLVE"
    os.system(rf"cd {root}\case_{i}\ && matlab -nosplash -nodesktop -r NONLINEAR_SOLVE")
    time.sleep(2)


#%%
import post

# %%
case = "".join(params) + str(size)
date = datetime.datetime.now().strftime("%y-%m-%d")
data_repo = root + rf"/LPRES_3000_{case}_{date}"

for i in range(nproc):
    case_dir = root + rf"\case_{i}"
    shutil.copytree(
        case_dir,
        data_repo + rf"\case{i}",
        ignore=shutil.ignore_patterns("*.cdb", "*.h5", ".git"),
    )

#%%
for i in range(nproc):
    shutil.rmtree(root + rf"\case_{i}")

# %%
