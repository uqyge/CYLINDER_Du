#%%
import numpy as np
import matplotlib.pyplot as plt

# %%
x = np.linspace(0, 20, 1000)
# %%
interval0 = [1 if (i < 5) else 0 for i in x]
interval1 = [1 if (i >= 5 and i < 10) else 0 for i in x]
interval2 = [1 if (i >= 10) else 0 for i in x]


# %%
y = x ** 2 * interval0 + x * interval1 + np.sin(x) * interval2
# %%
plt.plot(x, y)
# %%
def interval0(x):
    return 1 if (x < 5) else 0


def interval1(x):
    return 1 if (x >= 5 and x < 10) else 0


def interval2(x):
    return 1 if (x >= 10) else 0


def fun(x):
    return x ** 2 * interval0(x) + x * interval1(x) + np.sin(x) * interval2(x)


# %%
plt.plot(x, [fun(i) for i in x])
# %%
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(19680801)

mu = 200
sigma = 25
n_bins = 50
x = np.random.normal(mu, sigma, size=100)

fig, ax = plt.subplots(figsize=(8, 4))

# plot the cumulative histogram
n, bins, patches = ax.hist(
    x, n_bins, density=True, histtype="step", cumulative=True, label="Empirical"
)

# Add a line showing the expected distribution.
y = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu)) ** 2)
y = y.cumsum()
y /= y[-1]

ax.plot(bins, y, "k--", linewidth=1.5, label="Theoretical")

# Overlay a reversed cumulative histogram.
ax.hist(
    x, bins=bins, density=True, histtype="step", cumulative=-1, label="Reversed emp."
)

# tidy up the figure
ax.grid(True)
ax.legend(loc="right")
ax.set_title("Cumulative step histograms")
ax.set_xlabel("Annual rainfall (mm)")
ax.set_ylabel("Likelihood of occurrence")

plt.show()
# %%
y.cumsum()
# %%
import matplotlib.pyplot as plt
import numpy

data = numpy.random.randn(5)
print("The data is-", data)
sorted_random_data = numpy.sort(data)
p = 1.0 * numpy.arange(len(sorted_random_data)) / float(len(sorted_random_data) - 1)
print("The CDF result is-", p)

fig = plt.figure()
fig.suptitle("CDF of data points")
ax2 = fig.add_subplot(111)
ax2.plot(sorted_random_data, p)
ax2.set_xlabel("sorted_random_data")
ax2.set_ylabel("p")
# %%
plt.hist(data, cumulative=True, histtype="step", density=True)
# %%
p
# %%
