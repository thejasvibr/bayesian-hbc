#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Trying to analyse the horseshoe bat data using PyMC 3
Created on Wed Feb 24 05:53:20 2021

@author: Thejasvi Beleyur
"""

import arviz as az
import matplotlib.pyplot as plt
import natsort
import numpy as np 
import pandas as pd
import patsy
import scipy.stats as stats
import pymc3 as pm

# %% 
# Load and order rows by time of observation (indicated by video annotation ID)
d = pd.read_csv('one_call_per_row_2020-12-17.csv')
# fix the 21502300 --> 2123 to ease sorting 

d['timestamp'] = d['video_annot_id'].apply(lambda X: X[7:])
d['timestamp'] = d['timestamp'].str.replace('21502300','2123')
# annotation sorted
sorted_annots = d['timestamp'].sort_values(key=natsort.natsort_keygen())
d_timesorted = d.loc[sorted_annots.index,:].reset_index(drop=True)

# %% 
# Generate the dummy variable for multi-single bats
d_timesorted['group_state'] = d_timesorted['num_bats'].apply(lambda X: 0 if X==1 else 1)
d_timesorted = d_timesorted[~pd.isna(d_timesorted['tfm_duration'])].reset_index()
# %% 
# Build the linear model for tFM duration

lm = pm.Model()
with lm:
    baseline = pm.Uniform('baseline', lower=10**-4, upper=5*10**-3)
    alpha = pm.Uniform('alpha', lower=-5*10**-3,upper=5*10**-3)
    mu = pm.Deterministic('mu', baseline + alpha*d_timesorted['group_state'])
    sigma = pm.Uniform('sigma', lower=10**-4, upper=10**-3)
    pred_tfm = pm.Normal('pred_tfm', mu=mu, sd=sigma, observed=d_timesorted['tfm_duration'])
    trace_lm = pm.sample(2000, tune=1000)
    lm_map = pm.find_MAP()
    

# %% 
# Now generate predictions using the MAP estimates and inspect the residuals
d_timesorted['map_preds'] = lm_map['baseline'] + lm_map['alpha']*d_timesorted['group_state']
d_timesorted['resids'] = d_timesorted['map_preds']-d_timesorted['tfm_duration']

plt.figure()
plt.subplot(311)
plt.plot(d_timesorted['resids']*10**3,'*')
plt.subplot(312)
plt.plot(np.correlate(d_timesorted['resids'].dropna(),
                      d_timesorted['resids'].dropna(),'same'))
plt.subplot(313)
plt.hist(d['resids'])

# %%
plt.figure()
a0 = plt.subplot(111)
stats.probplot(d_timesorted['resids'], dist=stats.norm,
               sparams=(0,lm_map['sigma']),
               plot=a0)
plt.plot(np.linspace(-0.002,0.002,10),
         np.linspace(-0.002,0.002,10))
import seaborn as sns
from seaborn_qqplot import pplot

pplot(d_timesorted,x='resids', y = stats.norm)

