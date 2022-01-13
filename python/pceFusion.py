#%%
import chaospy
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from sklearn.metrics import r2_score
from scipy.special import kl_div

# %%
coordinates = np.linspace(0, 10, 1000)


def model_solver(parameters):
    alpha, beta, coordinates = parameters
    return alpha * np.e ** (-coordinates * beta)


for params in [(1.3, 0.13), (1.7, 0.17), (1.1, 0.19), (1.9, 0.11)]:
    samples = [(*params, coordinate) for coordinate in coordinates]
    plt.plot(coordinates, [model_solver(i) for i in samples])

# %%
alpha = chaospy.Normal(1.5, 0.2)
beta = chaospy.Uniform(0.1, 0.2)
gamma = chaospy.Uniform(0, 10)
joint = chaospy.J(alpha, beta, gamma)


order = 2
expansion = chaospy.generate_expansion(order, joint)
expansion.shape

# %%
samples = joint.sample(1000, rule="sobol")
evaluations = np.array([model_solver(sample) for sample in samples.T])
evaluations.shape

# %%
approx_solver = chaospy.fit_regression(expansion, samples, evaluations)

plt.scatter(samples[2], approx_solver(*samples))
#%%
length = coordinates.shape[0]
alpha_1, beta_1 = 1.5, 0.15
test_1 = np.vstack([[alpha_1] * length, [beta_1] * length, coordinates])

ax1 = plt.subplot(1, 2, 1)
ax1.plot(coordinates, approx_solver(*test_1), label="pce")
ax1.plot(coordinates, [model_solver(i) for i in test_1.T], label="true")
ax1.legend()
ax1.set_title(f"{alpha_1,beta_1=}")
ax2 = plt.subplot(1, 2, 2)
ax2.scatter([model_solver(i) for i in test_1.T], approx_solver(*test_1))
r2 = r2_score([model_solver(i) for i in test_1.T], approx_solver(*test_1)).round(5)
ax2.set_title(f"{order=},{r2=}")

# %%
n_samples, location = 1_000_000, 5
test = np.vstack([joint.sample(n_samples)[:2], np.ones(n_samples) * location])
pce_sampes = approx_solver(*test)
# plt.hist(pce_sampes, bins=20, density=True)

dist = chaospy.GaussianKDE(pce_sampes)
# %%
t = np.linspace(0, 1.5, 20)
plt.plot(t, dist.pdf(t))
plt.hist(pce_sampes, bins=20, density=True)

# %%
coeff = approx_solver.coefficients
indet = approx_solver.indeterminants
expon = approx_solver.exponents


def f0(coeff):
    pce_model = chaospy.sum(coeff * chaospy.prod(indet ** expon, axis=1))
    n_samples, location = 1_000, 5
    test = np.vstack([joint.sample(n_samples)[:2], np.ones(n_samples) * location])
    pce_sampes = pce_model(*test)

    return chaospy.GaussianKDE(pce_sampes)


# %%
x0 = [i * 0.86 for i in coeff]

plt.subplot(1, 2, 1)
plt.plot(t, f0(x0).pdf(t))
plt.plot(t, f0(coeff).pdf(t))
plt.subplot(1, 2, 2)
plt.plot(t, f0(x0).cdf(t))
plt.plot(t, f0(coeff).cdf(t))


# %%
exp = np.asarray([0.3, 0.5, 0.55, 0.6, 0.7, 0.95])
exp.sort()

y_cum = np.arange(1, len(exp) + 1) / len(exp)

# plt.hist(exp, cumulative=True, density=True, histtype="step")
# plt.plot(t, dist.cdf(t))
# plt.scatter(exp, y_cum)


def exp_cdf(x, epsilon):
    # epsilon = 0.1
    for i in range(len(exp)):
        if i == len(exp) - 1:
            return 1
        if x < exp[0]:
            return 0
        if (x >= exp[i]) & (x < exp[i + 1]):
            return max(min(y_cum[i] + epsilon, 1), 0)


# %%
t = np.linspace(0, 1.5, 200)
plt.plot(t, [exp_cdf(i, 0.3) for i in t], label="upper")
plt.plot(t, [exp_cdf(i, 0) for i in t])
plt.plot(t, [exp_cdf(i, -0.3) for i in t], label="lower")
# plt.hist(exp, cumulative=True, density=True, histtype="step")
plt.plot(t, f0(coeff).cdf(t), label="F0")
plt.plot(t, f0(x0).cdf(t), label="x0")
plt.legend()

# %%
def obj_func(x, exp=0):
    t = np.linspace(0, 1.5, 20)
    SMALL = 1e-10
    # return kl_divergence(f0(x).pdf(t) + SMALL, f0(coeff).pdf(t) + SMALL)
    return sum(kl_div(f0(x).pdf(t) + SMALL, f0(coeff).pdf(t) + SMALL)) / len(t)


print(f"{obj_func(x0)=},{obj_func(coeff)=}")

# %%
cons = (
    {
        "type": "ineq",
        "fun": lambda x, exp: exp_cdf(exp[0], epsilon=0.3) - f0(x).cdf(exp[0]),
        "args": (exp,),
    },
    {
        "type": "ineq",
        "fun": lambda x, exp: f0(x).cdf(exp[0]) - exp_cdf(exp[0], epsilon=-0.3),
        "args": (exp,),
    },
    {
        "type": "ineq",
        "fun": lambda x, exp: exp_cdf(exp[1], epsilon=0.3) - f0(x).cdf(exp[1]),
        "args": (exp,),
    },
    {
        "type": "ineq",
        "fun": lambda x, exp: f0(x).cdf(exp[1]) - exp_cdf(exp[1], epsilon=-0.3),
        "args": (exp,),
    },
    {
        "type": "ineq",
        "fun": lambda x, exp: exp_cdf(exp[2], epsilon=0.3) - f0(x).cdf(exp[2]),
        "args": (exp,),
    },
    {
        "type": "ineq",
        "fun": lambda x, exp: f0(x).cdf(exp[2]) - exp_cdf(exp[2], epsilon=-0.3),
        "args": (exp,),
    },
)

#%%
# res = minimize(obj_func, x0, args=(exp), method="trust-constr", constraints=cons,)
res = minimize(obj_func, x0, args=(exp), method="COBYLA", constraints=cons,)
res

# %%
plt.plot(t, [exp_cdf(i, 0.3) for i in t], "k--", label="upper")
plt.plot(t, [exp_cdf(i, 0) for i in t], "g")
plt.plot(t, [exp_cdf(i, -0.3) for i in t], "k--", label="lower")
plt.fill_between(
    t,
    [exp_cdf(i, 0.3) for i in t],
    [exp_cdf(i, -0.3) for i in t],
    alpha=0.3,
    color="grey",
)

plt.plot(t, f0(coeff).cdf(t), "b", label="f0")
plt.plot(t, f0(res.x).cdf(t), "r", label="opt")
# plt.plot(t, f0(x0).cdf(t), label="init")
plt.legend()

# %%
plt.plot(t, f0(coeff).pdf(t), "b", label="f0")
plt.plot(t, f0(res.x).pdf(t), "r", label="opt")
# plt.plot(t, f0(x0).pdf(t), label="init")
plt.legend()

# %%
