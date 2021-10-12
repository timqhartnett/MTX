#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 17:33:06 2021

composition CP plots

@author: timothy
"""
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np

files = glob.glob('/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/CP_Plots/Fe-*.csv')
mpl.rcParams['font.family'] = 'Arial'
dfs = []
fig = plt.figure(figsize = (10,6))
ax = plt.axes((0.1,0.2,0.8,0.7))
for i in range(len(files)):
    dfs.append(pd.read_csv(files[i]))
    upper = dfs[i]['mean']+dfs[i]['sd']
    lower = dfs[i]['mean']-dfs[i]['sd']
    ax.plot(dfs[i]['X'],dfs[i]['mean'],label=str(i)+'mean prediction')
    ax.fill_between(dfs[i]['X'],lower,upper,alpha = 0.2)
    ax.set_xlim(np.min(dfs[i]['X']),np.max(dfs[i]['X']))
    ax.tick_params(axis='y',labelsize=12)
    ax.tick_params(axis='x',labelsize=12)
    ax.set_xlabel('X substituion (%)',fontsize = 20)
    ax.set_ylabel('Transition Temperature (K)',fontsize = 20)
plt.legend(['M-site','T-site','Both sites'],fontsize=12,loc='upper right')

files = glob.glob('/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/CP_Plots/*-distribution.csv')
mpl.rcParams['font.family'] = 'Arial'
dfs = []
fig = plt.figure(figsize = (10,6))
ax = plt.axes((0.1,0.2,0.8,0.7))
for i in range(len(files)):
    dfs.append(pd.read_csv(files[i]))
    upper = dfs[i]['mean']+dfs[i]['sd']
    lower = dfs[i]['mean']-dfs[i]['sd']
    ax.plot(dfs[i]['X'],dfs[i]['mean'],label=str(i)+'mean prediction')
    ax.fill_between(dfs[i]['X'],lower,upper,alpha = 0.2)
    ax.set_xlim(np.min(dfs[i]['X']),np.max(dfs[i]['X']))
    ax.tick_params(axis='y',labelsize=12)
    ax.tick_params(axis='x',labelsize=12)
    ax.set_xlabel('X-substituion (%)',fontsize = 20)
    ax.set_ylabel('Transition Temperature (K)',fontsize = 20)
plt.legend(['Ge','Ga','Al','Sn'],fontsize=12,loc='upper right')