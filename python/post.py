#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
import scipy.io as sio

pio.renderers.default = "notebook_connected"

#%%


# %%
df_case = pd.DataFrame()
case = 0
for case in [0, 1, 2, 3, 4, 6, 7, 8, 9]:
    print(case)
    geo = pd.read_csv(f"../../case_{case}/geo.csv", names=["t", "R"], header=None)
    doe_res = sio.loadmat(f"../../case_{case}/doe_res.mat")

    for i in range(len(geo)):
        # i = 0
        # print(i, geo.iloc[i])

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

# plt.plot(df_incr.u, df_incr.load)
df_case.reset_index(inplace=True)
# df_case.tail()
# %%
for R in df_case["R"].unique():
    df_curve = df_case[df_case["R"] == R]
    # print(R,df_case[df_case['R']==R].u.shape)
    plt.plot(df_curve.u, df_curve.load)

# %%
case_incr_max = [
    df_case["incr"][df_case["R"] == unique_R].idxmax()
    for unique_R in df_case["R"].unique()
]

# %%
px.scatter_3d(df_case.iloc[case_incr_max], x="R", y="t", z="load")

#%%

# %%

# %%
len(df_case['R'].unique())
# %%
