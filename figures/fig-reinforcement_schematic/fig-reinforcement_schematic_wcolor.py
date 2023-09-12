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

ax.plot(df["strain"], df["stress"], marker="None", color="k", alpha=alpha)#, label="Neat")
is_make_vlines=0
if is_make_vlines:
  ax.vlines(25, np.amin(df["stress"]), np.amax(df["stress"]), color='k', ls='--')
  ax.vlines(55, np.amin(df["stress"]), np.amax(df["stress"]), color='k', ls='--')
  ax.vlines(150, np.amin(df["stress"]), np.amax(df["stress"]), color='k', ls='--')
  ax.vlines(230, np.amin(df["stress"]), np.amax(df["stress"]), color='k', ls='--')
ax.fill_betweenx([np.amin(df["stress"]), np.amax(df["stress"])], [0  , 0  ], [25 , 25 ], color='#E1CE7A')
ax.fill_betweenx([np.amin(df["stress"]), np.amax(df["stress"])], [25 , 25 ], [55 , 55 ], color='#FBFFB9')
ax.fill_betweenx([np.amin(df["stress"]), np.amax(df["stress"])], [55 , 55 ], [150, 150], color='#FDD692')
ax.fill_betweenx([np.amin(df["stress"]), np.amax(df["stress"])], [150, 150], [230, 230], color='#EC7357')
ax.fill_betweenx([np.amin(df["stress"]), np.amax(df["stress"])], [230, 230], [np.amax(df["strain"]), np.amax(df["strain"])], color='#754F44')

common_format(ax)
ax.set_ylabel(r"Stress (MPa), $\sigma$", fontsize=fontsize)
ax.set_xlabel(r"Strain (%), $\varepsilon$", fontsize=fontsize)
ax.set_xlim([np.amin(df["strain"]), np.amax(df["strain"])])
ax.set_ylim([np.amin(df["stress"]), np.amax(df["stress"])])
ax.set_yticks([0, 10, 20])
ax.set_xticks([0, 100, 200, 300])
ax.text((0  + 25 )/2                  , 0.92*np.amax(df["stress"]), "I"  , horizontalalignment="center")
ax.text((25 + 55 )/2                  , 0.92*np.amax(df["stress"]), "II" , horizontalalignment="center")
ax.text((55 + 150)/2                  , 0.92*np.amax(df["stress"]), "III", horizontalalignment="center")
ax.text((150+ 230)/2                  , 0.92*np.amax(df["stress"]), "IV" , horizontalalignment="center")
ax.text((230+ np.amax(df["strain"]))/2, 0.92*np.amax(df["stress"]), "V"  , horizontalalignment="center")

is_make_inset=0
if is_make_inset:
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
