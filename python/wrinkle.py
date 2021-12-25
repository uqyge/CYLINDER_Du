#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio

pio.renderers.default = "notebook_connected"
#%%
df = pd.read_hdf("./LPRES3000.h5")

# %%
for R in df["R"].unique():
    df_curve = df[df["R"] == R]
    # print(R,df[df['R']==R].u.shape)
    plt.plot(df_curve.u, df_curve.load)


# %%
grd = np.gradient(df_curve.load, df_curve.u)

for i in range(grd.size):
    if grd[i] / grd[i + 1] > 1.005:
        break
    # print(i,grd[i])
plt.figure()
plt.plot(grd[:i], "rd")
plt.plot(grd)
plt.figure()
plt.plot(df_curve.u, df_curve.load)
plt.plot(df_curve.u[:i], df_curve.load[:i], "rd")

#%%
idx = []
thres = 1.02
for R in df["R"].unique():
    df_curve = df[df["R"] == R]
    grd = np.gradient(df_curve.load, df_curve.u)

    for i in range(grd.size):
        if grd[i] / grd[i + 1] > thres:
            break
        id = df.index[df.u == df_curve.u.values[i]].values[0]
        # print(id)
    idx.append(id)
len(idx)
# %%
px.scatter_3d(df.iloc[idx], x="R", y="t", z="load")
# %%
case_incr_max = [
    df["incr"][df["R"] == unique_R].idxmax() for unique_R in df["R"].unique()
]
#%%
px.scatter_3d(df.iloc[case_incr_max], x="R", y="t", z="load")

# %%
idx = []
thres = 1.01
for R in df["R"].unique():
    df_curve = df[df["R"] == R]

    t = df_curve.t.unique()[0]
    x = np.linspace(df_curve.u.min(), df_curve.u.max(), 500)
    f = np.interp(x, df_curve.u, df_curve.load)
    grd = np.gradient(f, x)

    for i in range(grd.size):
        if grd[i] / grd[i + 1] > thres:
            break

    idx.append([f[i], t, R])

print(len(idx))
df_wrinkle = pd.DataFrame(np.asarray(idx), columns=["load", "t", "R"])
# %%
px.scatter_3d(df_wrinkle, x="t", y="R", z="load")
# %%
