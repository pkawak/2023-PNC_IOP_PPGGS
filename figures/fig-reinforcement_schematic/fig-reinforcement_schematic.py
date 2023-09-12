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

df = pd.read_csv("data.csv", header=0, delimiter=' ', skiprows=1)

fig, ax = plt.subplots(1, figsize=figsize, sharex=True, constrained_layout=True)

ax.plot(df["strain"], df["stress"], marker="None", color=neat_color  , alpha=alpha)#, label="Neat")
num_scatter = 20
#idx = np.round(np.linspace(0, len(df["Neat-Avg-poisy-"+kind]) - 1, num_scatter)).astype(int)
#ax.errorbar(df["Neat-Avg-shear"][idx], df["Neat-Avg-poisy-"+kind][idx], yerr=df["Neat-StdErr-poisy-"+kind][idx], ls='', marker=neat_marker  , color=neat_color  , alpha=alpha, label="Neat")
#ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)

common_format(ax)
ax.set_ylabel(r"Stress (MPa), $\sigma$", fontsize=fontsize)
ax.set_xlabel(r"Strain (%), $\varepsilon$", fontsize=fontsize)

# inset axes....
axins = ax.inset_axes([0.1, 0.7, 0.3, 0.25])
axins.plot(df["strain"], df["stress"], marker="None", color=neat_color  , alpha=alpha)#, label="Neat")
# subregion of the original image
x1, x2, y1, y2 = 0.0, 20.0, 0.0, 2.0
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.tick_params(axis='both', labelsize=labelsize-1, pad=0)
axins.set_ylabel(r"$\sigma$", fontsize=fontsize-1, labelpad=-6)
axins.set_xlabel(r"$\varepsilon$", fontsize=fontsize-1, labelpad=-6)
axins.set_yticks([2])
axins.set_xticks([0, 20])
#axins.set_xticklabels([])
#axins.set_yticklabels([])

#ax.indicate_inset_zoom(axins, edgecolor="black")

out_file_name = os.path.basename(sys.argv[0]).split('.')[0]+".pdf"
print("Writing", out_file_name)
fig.savefig(out_file_name)
