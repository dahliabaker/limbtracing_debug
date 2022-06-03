#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 13:52:50 2022

@author: dahliabaker
"""

from matplotlib import pyplot as plt



plt.style.use('ggplot')

plt.figure(figsize=(12,8))

plt.subplot(projection="aitoff")

plt.grid(True)

plt.xlabel('Long. in deg')
plt.ylabel('Lat. in deg')

plt.savefig('empty_aitoff.png',dpi=300)

#first, get data that's associated with cartesian coordinates

#convert cartesian to lat long

#create a data grid in lat long

#plot as a heat map surface on this aitoff projection

