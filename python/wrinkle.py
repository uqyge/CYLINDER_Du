#%%
import chaospy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
from chaospy.descriptives.expected import E
from sklearn import linear_model as lm
from sklearn.metrics import r2_score

pio.renderers.default = "notebook_connected"
#%%
df = pd.read_hdf("./LPRES3000_tRpL32_22-01-04.h5")

# %%
case_incr_max = [
    df["incr"][df["R"] == unique_R].idxmax() for unique_R in df["R"].unique()
]
doe = df.iloc[case_incr_max]
px.scatter_3d(doe, x="R", y="L", z="load")

# %%
for R in df["R"].unique():
    df_curve = df[df["R"] == R]
    plt.plot(df_curve.u, df_curve.load)

grd = np.gradient(df_curve.load, df_curve.u)

for i in range(grd.size):
    if grd[i] / grd[i + 1] > 1.05:
        break

plt.subplot(2, 1, 1)
plt.plot(df_curve.u[:i], grd[:i], "rd")
plt.plot(df_curve.u, grd)
plt.subplot(2, 1, 2)
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
print(f"{len(idx)=}")
px.scatter_3d(df.iloc[idx], x="R", y="p", z="load")

# %%
idx = []
thres = 1.01
for R in df["R"].unique():
    df_curve = df[df["R"] == R]

    t = df_curve.t.unique()[0]
    p = df_curve.p.unique()[0]
    x = np.linspace(df_curve.u.min(), df_curve.u.max(), 1000)
    f = np.interp(x, df_curve.u, df_curve.load)
    grd = np.gradient(f, x)

    for i in range(grd.size):
        if grd[i] / grd[i + 1] > thres:
            break

    idx.append([f[i], t, R, p])

print(f"{len(idx)=}")

df_wrinkle = pd.DataFrame(np.asarray(idx), columns=["load", "t", "R", "p"])
px.scatter_3d(df_wrinkle, x="t", y="R", z="load")

# %%
t = [75e-6, 75e-6]
R = [0.02, 0.04]
p = [1e4, 8e4]

dist_t = chaospy.Uniform(*t)
dist_r = chaospy.Uniform(*R)
dist_p = chaospy.Uniform(*p)

# joint = chaospy.J(dist_t, dist_r, dist_p)
joint = chaospy.J(dist_r, dist_p)
expansion = chaospy.generate_expansion(3, joint)
# %%
# samples = doe[["t", "R", "p"]].values
samples = doe[["R", "p"]].values
evaluations = doe["load"]

# %%
pce_model = chaospy.fit_regression(expansion, samples.T, evaluations)
doe["pce"] = pce_model(*samples.T)
plt.plot(doe.pce, doe.load, "d")
plt.title(f"{r2_score(doe.load,doe.pce)=}")
# %%
model = lm.Lars(fit_intercept=False)
lars_pce = chaospy.fit_regression(expansion, samples.T, evaluations, model=model)
doe["lars"] = lars_pce(*samples.T)

plt.plot(doe.load, doe.lars, "d")
plt.title(f"{r2_score(doe.load,doe.lars)=}")
# %%
px.scatter_3d(doe, x="p", y="R", z="lars")
# %%
px.scatter_3d(doe, x="p", y="R", z="pce")
#%%
px.scatter_3d(doe, x="p", y="R", z="load")

# %%
t_test = np.linspace(*t, 20)
r_test = 0.02 * np.ones(t_test.shape)
p_test = 1e4 * np.ones(t_test.shape)

test = np.vstack([t_test, r_test, p_test])
plt.plot(p_test, pce_model(*test))

# %%
r_test = np.linspace(*R, 20)
t_test = 25e-6 * np.ones(r_test.shape)
p_test = 3e4 * np.ones(t_test.shape)

test = np.vstack([t_test, r_test, p_test])
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
p = np.ones(r.shape) * 1e4

# %%
test = np.vstack([t, r, p])
pred = pce_model(*test)
# pred = lars_pce(*test)
# %%
opt = pd.DataFrame(np.vstack([t, r, pred]).T, columns=["t", "R", "pce"])
px.scatter_3d(opt, x="t", y="R", z="pce")

# %%
opt.iloc[opt["pce"].idxmin()]

# %%
opt.iloc[opt["pce"].idxmax()]

# %%
x = doe.p
y = doe.R
z = doe.load
plt.tricontour(x, y, z, 15, linewidths=0.5, colors="k")
plt.tricontourf(x, y, z, 15)
plt.colorbar()
# %%
# p*r**2 = m
p = np.linspace(2e4, 8e4, 100)
m = 0.04 ** 2 * p.min()
r = np.sqrt(m / p)
t = np.ones(r.shape) * 75e-6
#%%
# test = np.vstack([t, r, p])
test = np.vstack([r, p])
# pred = pce_model(*test)
pred = lars_pce(*test)
# %%
opt = pd.DataFrame(np.vstack([p, r, pred]).T, columns=["p", "R", "pce"])

# px.scatter_3d(opt[(opt.R > 0.02) & (opt.R < 0.035)], x="p", y="R", z="pce")
px.scatter_3d(opt, x="p", y="R", z="pce")

# %%
opt.pce.idxmax()

# %%
opt.p
# %%
opt.describe()
# %%
df.describe()
# %%
df.min()
# %%
df.max()
# %%
df[["R", "p"]]
# %%
df[df["R"] < 2.5e-2]


# %%
r2_score(doe.load, doe.pce)
# %%
r2_score(doe.load, doe.lars)
# %%

px.scatter_3d(doe, x="p", y="R", z="load")
# %%

px.scatter_3d(
    doe[(doe.R > 0.02) & (doe.R < 0.035) & (doe.p < 4e4)], x="p", y="R", z="lars"
)


# %%
doe['m']=doe.R**2*doe.p
# %%
px.scatter_3d(doe,x='R',y='m',z='load')

# %%
