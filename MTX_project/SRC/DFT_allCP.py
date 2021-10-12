#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
all CP plots DFT 9-23-21
"""
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np


for i in range(1,10):
    files = glob.glob('/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/CP_Plots/tx'+str(i)+'.cp*.csv')
    mpl.rcParams['font.family'] = 'Arial'
    dfs = []
    fig = plt.figure(figsize = (10,6))
    ax = plt.axes((0.1,0.2,0.8,0.7))
    prediction = pd.read_csv('/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/CP_Plots/tx-'+str(i)+'-predictions.csv')
    compositions = []
    for j in range(len(files)):
        composition = files[j].split('-')[3].split('.csv')[0]
        compositions.append(composition)
        dfs.append(pd.read_csv(files[j]))
        upper = dfs[j]['mean']+dfs[j]['sd']
        lower = dfs[j]['mean']-dfs[j]['sd']
        ax.plot(dfs[j]['Pnma.T.X'],dfs[j]['mean'],label=str(j)+'mean prediction')
        ax.fill_between(dfs[j]['Pnma.T.X'],lower,upper,alpha = 0.2)
        ax.set_xlim(np.min(dfs[j]['Pnma.T.X']),np.max(dfs[j]['Pnma.T.X']))
        ax.tick_params(axis='y',labelsize=12)
        ax.tick_params(axis='x',labelsize=12)
        ax.set_xlabel('Pnma T-X bond length ($\AA$)',fontsize = 20)
        ax.set_ylabel('Transition Temperature (K)',fontsize = 20)
    plt.axhline(297.232,color='r',linestyle='--')
    plt.scatter(x=prediction['V1'],y=prediction['V2'],color='r')
    compositions.append('BD intercept')
    plt.legend(compositions,fontsize=12)
    plt.savefig('/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/CP_Plots/DFT_cluster-'+str(i)+'-cp.png')
