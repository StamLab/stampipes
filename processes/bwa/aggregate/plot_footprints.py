#!/usr/bin/env python2
# coding: utf-8

# In[2]:


import sys

import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.stats
from footprint_tools.modeling import dispersion
from matplotlib.ticker import MaxNLocator
from pylab import rcParams

# In[3]:


# get_ipython().magic(u'matplotlib inline')
matplotlib.use("agg")


rcParams["pdf.fonttype"] = 42

# plt.switch_backend('agg')

# In[4]:


# dm = dispersion.read_dispersion_model(sys.argv[1])

dm = dispersion.read_dispersion_model(sys.argv[1])


# In[11]:


xx = np.arange(0, 100)
yy = np.array([dm.fit_mu(x) for x in xx])

r = np.array([dm.r[x] for x in xx])
p = np.array([dm.p[x] for x in xx])
mu = p * r / (1.0 - p)

fit_mu = np.array([dm.fit_mu(x) for x in xx])
fit_r = np.array([dm.fit_r(x) for x in xx])

fig = plt.figure()
gs = gridspec.GridSpec(1, 2, wspace=0.5)


ax = fig.add_subplot(gs[0, 0])

ax.plot(xx, mu, label="mle fit")
# ax.plot(xx, fit_mu, label = "smoothed parameters fit")
ax.plot([0, 100], [0, 100], color="grey", ls="--", zorder=-10)

ax.set_xlabel("Expected cleavages")
ax.set_ylabel("Mean (mu) observed cleavages")

[ax.spines[loc].set_color("none") for loc in ["top", "right"]]
ax.xaxis.set_ticks_position("bottom")
ax.xaxis.set_tick_params(direction="out")
ax.xaxis.set(major_locator=MaxNLocator(4))

ax.yaxis.set_ticks_position("left")
ax.yaxis.set_tick_params(direction="out")
ax.yaxis.set(major_locator=MaxNLocator(4))

ax.legend()


ax = fig.add_subplot(gs[0, 1])

ax.plot(xx[1:], 1 / r[1:])
ax.plot(xx[1:], 1 / fit_r[1:])

ax.set_xlabel("Expected cleavages")
ax.set_ylabel("1/r")

[ax.spines[loc].set_color("none") for loc in ["top", "right"]]
ax.xaxis.set_ticks_position("bottom")
ax.xaxis.set_tick_params(direction="out")
ax.xaxis.set(major_locator=MaxNLocator(4))

ax.yaxis.set_ticks_position("left")
ax.yaxis.set_tick_params(direction="out")
ax.yaxis.set(major_locator=MaxNLocator(4))

fig.set_size_inches(6, 2.5)

plt.savefig("dispersion.model.pdf", transparent=True)


# In[6]:


def step(arr, xaxis=False, interval=0):
    if xaxis and interval == 0:
        interval = abs(arr[1] - arr[0]) / 2.0
    newarr = np.array(zip(arr - interval, arr + interval)).ravel()
    return newarr


def fill_between(arr, ax, **kwargs):
    ax.fill_between(
        step(np.arange(arr.shape[0]), xaxis=True),
        step(np.zeros(arr.shape[0])),
        step(arr),
        **kwargs,
    )


def make_density_fit_plot(i, dm, ax, lo=0, hi=125, include_poisson=False):
    xx = np.arange(lo, hi)

    mu = dm.fit_mu(i)
    r = dm.fit_r(i)
    p = r / (r + mu)

    # ax.step(xx, dm.h[i, xx]/np.sum(dm.h[i,:]), label = "Observed cleavages")
    fill_between(
        dm.h[i, xx] / np.sum(dm.h[i, :]), ax, facecolor="lightgrey", edgecolor="none"
    )

    ax.plot(xx, scipy.stats.nbinom.pmf(xx, r, p), label="Negative binomial fit")

    if include_poisson:
        ax.plot(xx, scipy.stats.poisson.pmf(xx, i), label="Poisson (lambda = %d)" % i)

    ax.set_xlabel("Cleavages")
    ax.set_ylabel("Density")

    [ax.spines[loc].set_color("none") for loc in ["top", "right"]]
    ax.xaxis.set_ticks_position("bottom")
    ax.xaxis.set_tick_params(direction="out")
    ax.xaxis.set(major_locator=MaxNLocator(4))

    ax.yaxis.set_ticks_position("left")
    ax.yaxis.set_tick_params(direction="out")
    ax.yaxis.set(major_locator=MaxNLocator(4))


# In[9]:


fig = plt.figure()
gs = gridspec.GridSpec(1, 4, wspace=0.5)

for i, j in enumerate([5, 15, 25, 65]):
    ax = fig.add_subplot(gs[0, i])
    make_density_fit_plot(j, dm, ax, include_poisson=True)
    ax.set_xlim(right=125)


ax.legend()

fig.set_size_inches(12, 2)

plt.savefig("dispersion.hist.pdf", transparent=True)


# In[8]:


fig, ax = plt.subplots()

deltas = np.arange(0.5, 0, step=-0.1)

xx = np.arange(250, dtype=np.float, step=10)

for delta in deltas:
    y = dm.p_values(xx, xx * delta)
    plt.plot(xx, np.log10(y), label=delta)

ax.set_xlabel("Expected cleavage")
ax.set_ylabel("Significance (-log10 p-value) of depletion")

[ax.spines[loc].set_color("none") for loc in ["top", "right"]]
ax.xaxis.set_ticks_position("bottom")
ax.xaxis.set_tick_params(direction="out")
ax.xaxis.set(major_locator=MaxNLocator(6))

ax.yaxis.set_ticks_position("left")
ax.yaxis.set_tick_params(direction="out")
ax.yaxis.set(major_locator=MaxNLocator(4))

ax.grid(axis="y")

fig.set_size_inches(6, 3)

# ax.legend(bbox_to_anchor=(0., 1.05, 1., .105), loc=3,
#           ncol=5, mode="expand", borderaxespad=0., title = "Cleavage depletion ratio (obs/exp) at nucleotide")

plt.savefig("dispersion.power.analysis.pdf")


# In[ ]:
