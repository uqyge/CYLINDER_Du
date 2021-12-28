#%%
from chaospy.descriptives.expected import E
import chaospy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px

import plotly.io as pio

pio.renderers.default = "notebook_connected"
#%%
# df = pd.read_hdf("./LPRES3000.h5")
df = pd.read_hdf("./LPRES3000_32_new.h5")

# %%
for R in df["R"].unique():
    df_curve = df[df["R"] == R]
    plt.plot(df_curve.u, df_curve.load)


# %%
grd = np.gradient(df_curve.load, df_curve.u)

for i in range(grd.size):
    if grd[i] / grd[i + 1] > 1.05:
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
thres = 1.01
for R in df["R"].unique():
    df_curve = df[df["R"] == R]
    grd = np.gradient(df_curve.load, df_curve.u)

    for i in range(grd.size):
        if grd[i] / grd[i + 1] > thres:
            break
    id = df.index[df.u == df_curve.u.values[i]].values[0]
    # print(id)
    idx.append(id)
print(len(idx))
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
    x = np.linspace(df_curve.u.min(), df_curve.u.max(), 1000)
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


# %%
t = [25e-6, 50e-6]
R = [0.02, 0.04]

dist_t = chaospy.Uniform(*t)
dist_r = chaospy.Uniform(*R)
joint = chaospy.J(dist_t, dist_r)

expansion = chaospy.generate_expansion(8, joint)
# %%
doe = df.iloc[case_incr_max]

samples = doe[["t", "R"]].values
evaluations = doe["load"]

pce_model = chaospy.fit_regression(expansion, samples.T, evaluations)

# %%
doe["pce"] = pce_model(*samples.T)


# %%
px.scatter_3d(doe, x="t", y="R", z="pce")
# %%
px.scatter_3d(doe, x="t", y="R", z="load")
# %%
plt.plot(doe.pce, doe.load, "d")
# %%

# %%
from sklearn import linear_model as lm


#%%
model = lm.Lars(fit_intercept=False)
# %%
lars_pce = chaospy.fit_regression(expansion, samples.T, evaluations, model=model)

# %%
doe["lars"] = lars_pce(*samples.T)
# %%
plt.plot(doe.lars, doe.pce, "d")
# %%
plt.plot(doe.load, doe.pce, "d")
# %%
plt.plot(doe.load, doe.lars, "d")
# %%
px.scatter_3d(doe, x="t", y="R", z="lars")
# %%
px.scatter_3d(doe, x="t", y="R", z="pce")
# %%
t_test = np.linspace(*t, 20)
r_test = 0.02 * np.ones(t_test.shape)

test = np.vstack([t_test, r_test])
plt.plot(t_test, pce_model(*test))

# %%
r_test = np.linspace(*R, 20)
t_test = 20e-6 * np.ones(r_test.shape)

test = np.vstack([t_test, r_test])
plt.plot(r_test, pce_model(*test))

# %%
t_mean = 37.5e-6
R_mean = 0.03
L = 0.6

V = 2 * np.pi * R_mean * t_mean * L
print(f"{V=}")


def t_constraint(r):
    return V / (r * (2 * np.pi * L))


# %%
r = np.linspace(0.02, 0.04, 100)
t = t_constraint(r)
# %%
test = np.vstack([t, r])
pred = pce_model(*test)
# %%
opt = pd.DataFrame(np.vstack([t, r, pred]).T, columns=["t", "R", "pce"])

# %%
px.scatter_3d(opt, x="t", y="R", z="pce")
# %%
plt.plot(r, t)

# %%
opt.iloc[opt["pce"].idxmin()]

# %%
opt.iloc[opt["pce"].idxmax()]
# %%
px.scatter(df, x="u", y="load", color="R")
# %%
import matplotlib.tri as tri
import matplotlib.pyplot as plt

x = doe.t
y = doe.R
z = doe.load
plt.tricontour(x, y, z, 15, linewidths=0.5, colors="k")
plt.tricontourf(x, y, z, 15)
plt.colorbar()
# %%
