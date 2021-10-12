#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 19:19:19 2021

@author: timothy
"""
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np
CP_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/CP_Plots/homebrew_8-17-21/'

mpl.rcParams['font.family'] = 'Arial'
plt.ioff()

dataset = pd.read_csv("/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/CP_Plots/homebrew_8-17-21/predicted_composition_8-18-21.csv",
                      index_col=0)
dataset.index = range(dataset.shape[0])
for i in range(77):
    fig,ax = plt.subplots(nrows=3,ncols=3,figsize=(10,6),gridspec_kw={'wspace':0.4,'hspace':0.4})
    for j in range(3):
        for k in range(3):
            file = glob.glob(CP_DIR+'observation-'+str(i+1)+'-'+dataset.columns[j*3+k]+'.csv')
            fig_temp = plt.figure(figsize=(10,6))
            ax_temp = plt.axes((0.1,0.2,0.8,0.7))
            df = pd.read_csv(file[0],index_col=0)
            upper = df['mean']+df['sd']
            lower = df['mean']-df['sd']
            x = dataset.loc[i,df.columns[0]]
            y = dataset.loc[i,'Theating']
            ax[j][k].plot(df.iloc[:,0],df['mean'],label='mean prediction')
            ax[j][k].plot(x,y,'-ro',linewidth=0.5)
            ax[j][k].fill_between(df.iloc[:,0],lower,upper,alpha = 0.2)
            ax[j][k].set_xlim(np.min(df.iloc[:,0]),np.max(df.iloc[:,0]))
            ax[j][k].tick_params(axis='y',labelsize=8)
            ax[j][k].tick_params(axis='x',labelsize=8)
            ax[j][k].set_xlabel(df.columns[0]+' substitution (%)',fontsize = 12)
            ax[j][k].set_ylabel('T$_t$ (K)',fontsize = 12)
            ax_temp.plot(df.iloc[:,0],df['mean'],label='mean prediction')
            ax_temp.plot(x,y,'-ro',linewidth=0.5)
            ax_temp.fill_between(df.iloc[:,0],lower,upper,alpha = 0.2)
            ax_temp.set_xlim(np.min(df.iloc[:,0]),np.max(df.iloc[:,0]))
            ax_temp.tick_params(axis='y',labelsize=12)
            ax_temp.tick_params(axis='x',labelsize=12)
            ax_temp.set_xlabel(df.columns[0]+' substitution (%)',fontsize = 20)
            ax_temp.set_ylabel('T$_t$ (K)',fontsize = 20)
            fig_temp.savefig(CP_DIR+'individual/observation-'+str(i+1)+'-'+df.columns[0]+'.png')
    fig.savefig(CP_DIR+'total/observation-'+str(i+1)+'-total.png')
            
            
'''
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
'''