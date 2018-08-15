#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 11:53:36 2018

@author: fourier
"""

import pandas as pd
import numpy as np
import os

path = '/home/fourier/Documentos/Chicago_Crimes'

os.chdir(path)


C1 = pd.read_csv('Chicago_Crimes_2001_to_2004.csv', error_bad_lines=False)



C1 = C1.fillna(0)

C1["Date"] = pd.to_datetime(C1["Date"])

C1.drop("Case Number", axis=1, inplace=True)
C1.drop("Block", axis=1, inplace=True)
C1.drop("Description", axis=1, inplace=True)
C1.drop("IUCR", axis=1, inplace=True)
C1.drop("Unnamed: 0", axis=1, inplace=True)
C1.drop("Location", axis=1, inplace=True)
C1.drop("Updated On", axis=1, inplace=True)
C1.drop("ID", axis=1, inplace=True)
C1.drop("Year", axis=1, inplace=True)


Types = {}

for i in range(len(C1)):
  temp_fbi = C1["FBI Code"][i]
  temp_type = C1["Primary Type"][i]
  if not temp_fbi in Types.keys():
    Types[temp_fbi] = []
  if not temp_type in Types[temp_fbi]:
    Types[temp_fbi].append(temp_type)
    
  

C1.drop("Primary Type", axis=1, inplace=True)
C1.drop("X Coordinate", axis=1, inplace=True)
C1.drop("Y Coordinate", axis=1, inplace=True)
C1.drop("Location Description", axis=1, inplace=True)
C1.drop("Arrest", axis=1, inplace=True)
C1.drop("Domestic", axis=1, inplace=True)


C_test = C1.copy()
C_test["Quantity"] = 1
C_test = C_test.groupby([pd.Grouper(key='Date', freq='D'), 'District','Ward', 'FBI Code'])['Quantity'].sum().reset_index().sort_values('Date')
.0
C_test = C_test.groupby('FBI Code').resample('W-Mon', on='Date').sum().reset_index().sort_values(by='Date')
