#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on Fri Sept 8 05:48:00 2023

@author: pierrekawak
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('tkagg')
import numpy as np
import sys

def common_format(ax):
  ax.tick_params(axis='both', labelsize=labelsize, pad=0)
#  ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)
  if test == 0:
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)

# figure settings
test = 1
fontsize         = 9
labelsize        = 8
figsize          = (3.25, 2.5)
pt_to_in         = 72
neat_marker      = "^"
filled_marker    = "o"
neat_color       = "r"
filler_color     = "b"
filled_color     = "#808080"
polymer_color    = "r"
other_color      = "c"
alpha = 1.0

def simple_derivative(x, y): 
  der = []
  for i in range(1,len(x)-1):
    der.append((y[i+1]-y[i-1])/(x[i+1]-x[i-1]))
    if i == 1:
      der.append((y[i+1]-y[i-1])/(x[i+1]-x[i-1]))
  der.append((y[i+1]-y[i-1])/(x[i+1]-x[i-1]))
  return(der)

df = pd.read_csv("data.csv", header=0, delimiter=' ', skiprows=1)
df["modulus"] = simple_derivative(np.array(df["strain"]), np.array(df["stress"]))

fig, ax = plt.subplots(1, figsize=figsize, sharex=True, constrained_layout=True)

ax.plot(df["strain"].iloc[1:-1], df["modulus"].iloc[1:-1], marker="None", color=neat_color  , alpha=alpha)#, label="Neat")
ax.vlines(25, np.amin(df["modulus"]), np.amax(df["modulus"]), color='k', ls='--')
ax.vlines(55, np.amin(df["modulus"]), np.amax(df["modulus"]), color='k', ls='--')
ax.vlines(150, np.amin(df["modulus"]), np.amax(df["modulus"]), color='k', ls='--')
ax.vlines(230, np.amin(df["modulus"]), np.amax(df["modulus"]), color='k', ls='--')

common_format(ax)
ax.set_ylabel(r"Modulus, $E=\mathrm{d}\sigma/\mathrm{d}\varepsilon$", fontsize=fontsize)
ax.set_xlabel(r"Strain (%), $\varepsilon$", fontsize=fontsize)

out_file_name = os.path.basename(sys.argv[0]).split('.')[0]+".pdf"
print("Writing", out_file_name)
fig.savefig(out_file_name)
