#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 13:36:38 2021

@author: timothy
"""
import pandas as pd
import numpy as np
from sqlalchemy import create_engine
import pymysql
import mysql.connector as sql
import re
import math

MnNiSi_master= pd.read_csv('/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Data/MASTER.csv')
references = MnNiSi_master['Reference']
families = MnNiSi_master['Family']
def remove_parenthesis(array):
    updated_list = []
    for value in array:
        if type(value) == str:
            updated_list.append(re.sub('[()]','',value))
        else:
            updated_list.append(value)
    return updated_list

references_new = remove_parenthesis(references)
families_new = remove_parenthesis(families)
MnNiSi_master['Reference'] = references_new
MnNiSi_master['Family'] = families_new

table_name = 'VCU'
user = 'tqh8pn'
pwd = 'Hart7355'
host = 'localhost'
db = 'mtx'

engine = create_engine("mysql+pymysql://{user}:{pw}@{host}/{db}?charset=utf8mb4".format(host=host, db=db, user=user, pw=pwd))
MnNiSi = MnNiSi_master.drop(['Reference','Family'],axis=1) ### need to find a way to get the references inserted
# currently says invalid string likely due to wrong character set?
MnNiSi.to_sql(table_name,engine,index=False, if_exists = 'replace')
#features_2 = pd.read_sql("select * from Magpy", engine)