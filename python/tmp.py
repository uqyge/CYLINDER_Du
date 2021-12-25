#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%
df = pd.read_hdf("./LPRES1000.h5")

# %%
for R in df["R"].unique():
    df_curve = df[df["R"] == R]
    # print(R,df[df['R']==R].u.shape)
    plt.plot(df_curve.u, df_curve.load)

# plt.show()

# %%
plt.plot(np.gradient(df_curve.u))

#%%
plt.plot(df_curve.u)
# %%
grd = np.gradient(df_curve.load,df_curve.u)
grd.size

#%%
for i in range(grd.size):
    if grd[i]/grd[i+1] >1.005:
        break
    # print(i,grd[i])
plt.figure()
plt.plot(grd[:i],'rd')
plt.plot(grd)
plt.figure()
plt.plot(df_curve.u,df_curve.load)
plt.plot(df_curve.u[:i],df_curve.load[:i],'rd')


#%%

plt.plot(df_curve.u, 
grd,
'd-')
# %%
plt.plot(df_curve.u,df_curve.load,'d-')

#%%
