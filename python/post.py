#%%
import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import scipy.io as sio

pio.renderers.default = "notebook_connected"

# %%
df_case = pd.DataFrame()

# for case in [0, 1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15]:
# for case in [8]:
for case in range(16):
    print(case)
    geo = pd.read_csv(
        f"../../case_{case}/geo.csv",
        # names=["t", "R"],
        names=["t", "R", "p", "L"],
        header=None,
    )
    doe_res = sio.loadmat(f"../../case_{case}/doe_res.mat")

    for i in range(len(geo)):
        # print(i, geo.iloc[i])
        if i < doe_res["doe_res"]["u"].shape[1]:
            geo_val = geo.values[i]
            u = -doe_res["doe_res"]["u"][0, i][1::].reshape(-1)
            load = -doe_res["doe_res"]["load"][0, i][1::].reshape(-1)
            incr = np.asarray(range(len(u)))
            incr_u_load = np.vstack([incr, u, load]).T

            FAILDISP = (doe_res["doe_res"]["res"][0, i]["FAILDISP"],)
            FAILLOAD = (doe_res["doe_res"]["res"][0, i]["FAILLOAD"],)

            # print(FAILDISP, FAILLOAD)
            # print(u[-1], load[-1])

            assert (FAILLOAD == load[-1]) & (FAILDISP == u[-1])
            # plt.plot(u, load)
            # plt.show()

            data = np.hstack([incr_u_load, np.ones(len(u)).reshape(-1, 1) * geo_val])
            df_incr = pd.DataFrame(
                data, columns=["incr", "u", "load"] + geo.columns.tolist()
            )
            df_case = pd.concat([df_case, df_incr])

df_case.reset_index(inplace=True)
# df_case.tail()

# %%
for R in df_case["R"].unique():
    df_curve = df_case[df_case["R"] == R]
    # print(R,df_case[df_case['R']==R].u.shape)
    plt.plot(df_curve.u, df_curve.load)
    # px.scatter(df_curve, x="u", y="load")
plt.show()

# %%
case_incr_max = [
    df_case["incr"][df_case["R"] == unique_R].idxmax()
    for unique_R in df_case["R"].unique()
]
print(f"{len(case_incr_max)=}")
# %%
px.scatter_3d(df_case.iloc[case_incr_max], x="R", y="p", z="load")

# %%
case = "".join(geo.columns) + str(len(case_incr_max))
date = datetime.datetime.now().strftime("%y-%m-%d")
filename = f"LPRES3000_{case}_{date}.h5"
print(f"{filename=}")
df_case.to_hdf(filename, key="df_case")

# %%
