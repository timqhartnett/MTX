#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 19:58:34 2021

@author: timothy
"""
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import numpy as np

BD_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/breakdown/'
plt.ioff()
mpl.rcParams['font.family'] = 'Arial'
breakdown = pd.read_csv(BD_DIR+'breakdown_composition_8-17-21.csv',index_col=0)
breakdown.index = range(breakdown.shape[0])
cols = list(breakdown.columns.values)
breakdown = breakdown[['intercept','Mn','Fe','Co','Ni','Si','Ge','Sn','Al','Ga','prediction']]
breakdown = breakdown.rename({'intercept': 'Intercept','prediction': 'Prediction'},axis=1)

for i in range(breakdown.shape[0]):
    fig = plt.figure(figsize=(10,6))
    ax = plt.axes((0.1,0.2,0.8,0.7))
    contribution = np.array(breakdown.loc[i,:])
    red_value = contribution[contribution<0]
    blue_value = contribution[contribution>=0]
    red_y = np.where(contribution<0)[0]
    blue_y = np.where(contribution>=0)[0]
    ypos = np.arange(len(contribution))
    ax.barh(red_y,red_value,color='red')
    ax.barh(blue_y,blue_value,color='blue')
    ax.set_yticks(ypos)
    ax.set_yticklabels(breakdown.columns)
    ax.invert_yaxis()
    ax.set_xlabel('Predicted T$_t$ (K)',fontsize = 20)
    ax.tick_params(axis='y',labelsize=12)
    ax.axvline(x=0,color='black',linestyle='--',linewidth=2)
    for j in range(len(blue_value)):
        ax.text(blue_value[j]+3,blue_y[j]+0.2,str(np.round(blue_value[j],decimals=2)),color='blue')
    for j in range(len(red_value)):
        ax.text(3,red_y[j],str(np.round(red_value[j],decimals=2)),color='red')      
    plt.xlim([np.min(contribution)-50,np.max(contribution)+50])
    fig.savefig(BD_DIR+'homebrew/observation-'+str(i+1)+'.png')