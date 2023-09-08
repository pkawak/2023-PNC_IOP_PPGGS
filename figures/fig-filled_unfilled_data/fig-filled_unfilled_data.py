#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on Fri March 3 05:48:00 2023

@author: pierrekawak
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('tkagg')
import pandas as pd

test = 0
fontsize  = 12
labelsize = 10
figsize   = (2.75, 2.5)

neat_file_name   = "Unfilled.csv"
filled_file_name = "CB.csv"

neat_df          = pd.read_csv(neat_file_name, index_col=False, delimiter=',', header=None) 
filled_df        = pd.read_csv(filled_file_name, index_col=False, delimiter=',', header=None)

neat_df   = neat_df.sort_values(by=neat_df.columns[0])
filled_df = filled_df.sort_values(by=filled_df.columns[0])

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)

ax.plot(neat_df.iloc[:, 0]/100, neat_df.iloc[:, 1], color='r', label="Neat Rubber")
ax.plot(filled_df.iloc[:, 0]/100, filled_df.iloc[:, 1], color='#808080', label="Filled with CB")

ax.tick_params(axis='both', labelsize=labelsize, pad=0)
#ax.set_yscale('log')
ax.set_ylabel(r"Stress, $\sigma$ (MPa)", fontsize=fontsize)#, labelpad=10)
ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)#, labelpad=10)
ax.legend(loc='best', fontsize=labelsize, handlelength=1.3)

#ax.set_xlim(0, 450)
#ax.set_ylim(0, 23)

ax.scatter(neat_df.iloc[-1, 0]  /100, neat_df.iloc[-1, 1], color='r', marker='x')
ax.scatter(filled_df.iloc[-1, 0]/100, filled_df.iloc[-1, 1], color='#808080', marker='x')

if test == 0:
  plt.subplots_adjust(left=0.167, bottom=0.162, right=0.998, top=0.998)
import os
basename = os.path.splitext(__file__)[0]
plt.savefig(basename+".pdf")
plt.savefig(basename+".png")

ax.set_xlim(ax.get_xlim()[0], ax.get_xlim()[1])
ax.plot([ax.get_xlim()[0], ax.get_xlim()[1]], [2, 2], color='k', ls='--', label="Loading at 2 MPa")
ax.legend(loc='best', fontsize=labelsize, handlelength=1.3)
plt.savefig(basename+".2.pdf")
