#%%
import os
import shutil

import numpy as np
import sobol_seq

#%%

t = [25e-6, 75e-6]
R = [0.01, 0.04]
c = [0, 1]

data = np.vstack([t, R, c])

size = 10
vec = sobol_seq.i4_sobol_generate(data.shape[0], size)

lb = data[:, 0]
ub = data[:, 1]

out = vec * (ub - lb) + lb
# print(out)
#%%
# source_dir = r"D:\projects\CYLINDER_Du"
# source_dir = r"C:\Users\admin\workspace\matlab\CYLINDER_Du"
root = "C:\\work\\workspace\\matlab\\"
source_dir = root + "CYLINDER_Du"

nproc = 10

# for i in range(4):
for i, geo in enumerate(np.array_split(out, nproc)):
    # print(i, geo)

    # destination_dir = r"D:\projects\case_" + str(i)
    # destination_dir = r"C:\Users\admin\workspace\matlab\case_" + str(i)
    destination_dir = root + r".\case_" + str(i)
    shutil.copytree(source_dir, destination_dir)

    np.savetxt(destination_dir + "/geo.csv", geo, delimiter=",")

# %%
for i in range(nproc):
    # f"cd d:\projects\case_{i}\ && matlab -nosplash -nodesktop -r GENERATE_CDB"
    os.system(
        rf"cd C:\Users\admin\workspace\matlab\case_{i}\ && matlab -nosplash -nodesktop -r GENERATE_CDB"
    )

# %%
for i in range(nproc):
    # f"cd d:\projects\case_{i}\ && matlab -nosplash -nodesktop -r NONLINEAR_SOLVE"
    os.system(
        rf"cd C:\Users\admin\workspace\matlab\case_{i}\ && matlab -nosplash -nodesktop -r NONLINEAR_SOLVE"
    )

#%%
for i in range(nproc):
    shutil.rmtree(root + rf"\case_{i}")
# %%
shutil.rmtree(root + r"\case_0")
# %%
