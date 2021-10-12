#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 07:41:23 2021

@author: timothy
"""
import pandas as pd
import sys
sys.path.insert(1,'/home/timothy/magpie-python')
from magpie import MagpieServer
import numpy as np

DATA_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Data/'

ele
m = MagpieServer()
x = m.generate_attributes('oqmd-dH', elements)
x.index = elements

reduced = x[['mean_CovalentRadius','mean_MeltingT','mean_AtomicWeight','mean_MendeleevNumber','mean_NdValence']]
reduced.to_csv(DATA_DIR+'element_magpie.csv')
