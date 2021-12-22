#%%
import scipy.io as sio
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#%%
geo = pd.read_csv("../../case_0/geo.csv", names=["t", "R"], header=None)
geo

# %%
doe_res = sio.loadmat("../../case_0/doe_res.mat")

# %%
doe_res

# %%
# plt.plot(-doe_res["doe_res"]["u"][0, 0], -doe_res["doe_res"]["load"][0, 0])
i=39
print(i,geo.iloc[i])
u = -doe_res["doe_res"]["u"][0, i]
load = -doe_res["doe_res"]["load"][0, i]
plt.plot(u, load)
# %%
doe_res["doe_res"]["u"][0, :].shape
# %%

# %%
geo.iloc[39]
# %%
u = doe_res["doe_res"]["u"][0,:]
load = doe_res["doe_res"]["load"][0,:]
# %%
u[0].reshape(-1)
load[0].reshape(-1)
# %%
load[0]
# %%
a = np.vstack([-u[0].reshape(-1),-load[0].reshape(-1)]).T
a.shape
# %%

# %%
data = np.hstack([a,np.ones([a.shape[0],1])*geo.values[0]])
data.shape
# %%
plt.plot(-data[:,0],-data[:,1])
# %%

# %%
