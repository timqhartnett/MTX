#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 12:03:09 2021

@author: timothy
"""

import pandas as pd
import sys
sys.path.insert(1,'/home/timothy/magpie-python')
from magpie import MagpieServer
import numpy as np

import os
os.chdir('/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/SRC')

import Virtual_set_creation as vsc

'''directory definitions'''
DATA_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Data/'

'''grab datasets'''
elementwise = pd.read_excel(DATA_DIR+'MnNiSi_VCU_Dataset_for_ML.xlsx',engine='openpyxl',index_col=0)
elementwise.to_csv(DATA_DIR+'MnNISi_VCU_ML_3-16-21.csv')
elementwise.drop_duplicates()
Compositions = elementwise['1A'] 
for i in range(11):
    Compositions = Compositions + elementwise.iloc[:,i+1].astype(str) 
Compositions = Compositions.to_list()
Compositions.append('Mn22.3Ni22.2Si16.65Fe22.2Ge16.65')
Compositions.append('Mn20Ni26.7Si11Fe26,7Ga15.6')
Compositions.append('MnNiSi0.2Ge0.2Sn0.2Al0.2Ga0.2')
Compositions.extend(vsc.create_end_points())
Compositions.extend(vsc.create_high_entropic())


m = MagpieServer()
x = m.generate_attributes('oqmd-dH', Compositions)
x.index = Compositions
Theating = elementwise['Theating'].to_list()
Tcooling = elementwise['Tcooling'].to_list()

#reduced['Tcooling'] = elementwise['Tcooling'].copy().to_list().append('NaN')
#reduced.to_csv('/home/timothy/Desktop/MnNiSi/transition_temp_modeling/MagpieFeatures.csv')

comps_sheet  = pd.read_excel(DATA_DIR+'composition_3-20-21.xlsx',engine='openpyxl',index_col=0)
comps_sheet = comps_sheet.round(3)
new_comps = 'Mn'+comps_sheet['Mn'].astype(str)
for element in comps_sheet.columns[1:9]:
    new_comps = new_comps+element+comps_sheet[element].astype(str)
new_comps = new_comps.tolist()
new_comps.extend(vsc.create_high_entropic())
new_comps.extend(vsc.create_Fe_Ni())
m = MagpieServer()
m2 = MagpieServer()
y = m.generate_attributes('oqmd-dH',new_comps)
y.index = new_comps
z = y.corr()
numeric = y.columns
var = y.var()
variable = []
for i in range(0,len(var)):
    if var[i] >= 1:
        variable.append(numeric[i])
        
Theating = comps_sheet['Theating'].to_list()
Tcooling = comps_sheet['Tcoolilng'].to_list()
reduced = y[['mean_CovalentRadius','dev_MeltingT','mean_AtomicWeight','mean_MendeleevNumber','mean_NdValence']]
#virtual = virt[['mean_CovalentRadius','dev_MeltingT','mean_AtomicWeight','mean_MendeleevNumber','mean_NdValence']]
Theating.extend([np.nan]*(reduced.shape[0]-len(Theating)))
Tcooling.extend([np.nan]*(reduced.shape[0]-len(Tcooling)))
reduced['Theating'] = Theating
reduced['Tcooling'] = Tcooling
reduced['composition'] = new_comps

reduced.iloc[:82,:].to_csv(DATA_DIR+'MagpyFeatures_3-20-21.csv')
reduced.iloc[82:,:].to_csv(DATA_DIR+'Magpie_virtual_3-25-21.csv')
''' homewbrewed '''
m = MagpieServer()
comps_sheet  = pd.read_excel(DATA_DIR+'composition_3-20-21.xlsx',engine='openpyxl',index_col=0)
comps_sheet = comps_sheet.round(3)
element_descriptors = m.generate_attributes('oqmd-dH', ['Mn','Ni','Fe','Co','Si','Ge','Sn','Al','Ga','In',
                                                        ])
element_descriptors = element_descriptors[['mean_CovalentRadius','mean_MeltingT','mean_AtomicWeight','mean_MendeleevNumber','mean_NdValence']]
element_descriptors.index = ['Mn','Ni','Fe','Co','Si','Ge','Sn','Al','Ga','In']
x = vsc.create_compositions_DF()
In = vsc.create_compositions_In_DF()
Theating = elementwise['Theating'].dropna().to_list()
Tcooling = elementwise['Tcooling'].dropna().to_list()
Theating.extend([np.nan]*(comps_sheet.shape[0]-len(Theating)))
Tcooling.extend([np.nan]*(comps_sheet.shape[0]-len(Tcooling)))

comps = pd.DataFrame(np.zeros([comps_sheet.shape[0],element_descriptors.shape[1]]),columns=element_descriptors.columns)
for feature in element_descriptors.columns:
    tempDF = pd.DataFrame([comps_sheet[el]*element_descriptors.loc[el,feature] for el in comps_sheet.columns[:9]]).T
    comps[feature] = np.array(tempDF.sum(axis=1))
comps['Theating'] = Theating
comps['Tcooling'] = Tcooling

comps.to_excel(DATA_DIR+'MagpyFeatures_3-21-21.xlsx')
norm = vsc.create_compositions_DF()
norm.to_csv(DATA_DIR+'virtual_norm.csv')
high_entropy = vsc.create_composition_High_Entropy()
high_entropy.to_csv(DATA_DIR+'virtual_entropy.csv')

virtual_high_entropy =  pd.DataFrame(np.zeros([high_entropy.shape[0],element_descriptors.shape[1]]),columns=element_descriptors.columns)
for feature in element_descriptors.columns:
    tempDF = pd.DataFrame([high_entropy[el]*element_descriptors.loc[el,feature] for el in high_entropy.columns[:10]]).T
    virtual_high_entropy[feature] = np.array(tempDF.sum(axis=1))

virtual_In =  pd.DataFrame(np.zeros([In.shape[0],element_descriptors.shape[1]]),columns=element_descriptors.columns)
for feature in element_descriptors.columns:
    tempDF = pd.DataFrame([In[el]*element_descriptors.loc[el,feature] for el in In.columns[:10]]).T
    virtual_In[feature] = np.array(tempDF.sum(axis=1))

virtual = pd.concat([virtual_In,virtual_high_entropy],ignore_index=True)
virtual.to_excel(DATA_DIR+'virtual.xlsx')
