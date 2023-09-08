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

# figure settings
test = 1
fontsize  = 9
labelsize = 8
figsize   = (2.65, 2.6)
bin       = 1
num_bins  = 20
pt_to_in  = 72
filter_shear = 0.2
#figsize   = (180/pt_to_in, 190/pt_to_in)
neat_marker = "^"
filled_marker = "o"

# constants
Npolymer = 102500
Nfiller  = 51450
Ntotal   = Npolymer+Nfiller

neat_file_name   = "../data_files/neat_deform_dataframe.out"
filled_file_name = "../data_files/filled_deform_dataframe.out"
neat_base_file_name   = "../data_files/neat_npt_dataframe.out"
filled_base_file_name = "../data_files/filled_npt_dataframe.out"

neat_df          = pd.read_csv(neat_file_name)
filled_df        = pd.read_csv(filled_file_name)
neat_base_df     = pd.read_csv(neat_base_file_name)
filled_base_df   = pd.read_csv(filled_base_file_name)

neat_df          = neat_df[neat_df["v_shear_shift"] < filter_shear]
filled_df        = filled_df[filled_df["v_shear_shift"] < filter_shear]
neat_base_df     = neat_base_df[neat_base_df["v_shear_shift"] < filter_shear]
filled_base_df   = filled_base_df[filled_base_df["v_shear_shift"] < filter_shear]

neat_base_df     = neat_base_df.mean()
filled_base_df   = filled_base_df.mean()

if bin:
  bins = np.linspace(np.amin(neat_df["Time"]), np.amax(neat_df["Time"])+1, int(len(neat_df["Time"])/num_bins))
  neat_df["Time_Binned"] = pd.cut(neat_df["Time"], bins, right=False)
  neat_df = neat_df.groupby("Time_Binned", as_index=False).mean()

  bins = np.linspace(np.amin(filled_df["Time"]), np.amax(filled_df["Time"])+1, int(len(filled_df["Time"])/num_bins))
  filled_df["Time_Binned"] = pd.cut(filled_df["Time"], bins, right=False)
  filled_df = filled_df.groupby("Time_Binned", as_index=False).mean()

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)

neat_shear    = np.append(0, neat_df["v_shear_shift"])
filled_shear  = np.append(0, filled_df["v_shear_shift"])
neat_stress   = np.append(neat_base_df["v_pxx"], neat_df["v_pxx"])
filled_stress = np.append(filled_base_df["v_pxx"], filled_df["v_pxx"])
ax.plot(neat_shear  , neat_stress  , marker=neat_marker, label="Neat")
ax.plot(filled_shear, filled_stress, marker=filled_marker, label="Filled")

#ax.set_xticks([0, 0.5, 1.0, 1.5, 2])
ax.tick_params(axis='both', labelsize=labelsize, pad=0)
ax.set_ylabel(r"Stress, $\sigma$", fontsize=fontsize)#, labelpad=10)
ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)#, labelpad=10)
ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)

#ax.set_xlim(0, 450)
#ax.set_ylim(0, 23)
if test == 0:
  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)
import os
basename = os.path.splitext(__file__)[0]
plt.savefig(basename+".pdf")

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)

# gotta add the quiescent at 0
neat_shear   = np.append(0, neat_df["v_shear_shift"])
filled_shear = np.append(0, filled_df["v_shear_shift"])
neat_stress   = np.append(neat_base_df["v_pxx"], neat_df["v_pxx"])
filled_stress = np.append(filled_base_df["v_pxx"], filled_df["v_pxx"])

filled_polymer_stress = np.append(filled_base_df["c_pxx_polymer[1]"]/Npolymer, filled_df["c_pxx_polymer[1]"]/Npolymer)
filled_filler_stress  = np.append(filled_base_df["c_pxx_filler[1]"]/Nfiller, filled_df["c_pxx_filler[1]"]/Nfiller)

ax.plot(neat_shear  , neat_stress  , marker=neat_marker, alpha=0.3, label="Neat")
ax.plot(filled_shear, filled_stress, marker=filled_marker, label="Filled")

ax.plot(filled_shear, filled_filler_stress , marker=filled_marker, label="Filled; Filler")
ax.plot(filled_shear, filled_polymer_stress, marker=filled_marker, label="Filled; Polymer")

#ax.set_xticks([0, 0.5, 1.0, 1.5, 2])
ax.tick_params(axis='both', labelsize=labelsize, pad=0)
ax.set_ylabel(r"Stress, $\sigma$", fontsize=fontsize)#, labelpad=10)
ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)#, labelpad=10)
ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)

#ax.set_xlim(0, 450)
#ax.set_ylim(0, 23)
if test == 0:
  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)
plt.savefig(basename+".groups.pdf")

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)

ax.plot(filled_shear, filled_stress-neat_stress, marker=filled_marker, label="Total")

ax.plot(filled_shear, filled_filler_stress , marker=filled_marker, label="Filler")
ax.plot(filled_shear, filled_polymer_stress-neat_stress, marker=filled_marker, label="Polymer")

#ax.set_xticks([0, 0.5, 1.0, 1.5, 2])
ax.tick_params(axis='both', labelsize=labelsize, pad=0)
ax.set_ylabel(r"$\sigma_{\mathrm{Filled}} - \sigma_{\mathrm{Neat}}$", fontsize=fontsize)#, labelpad=10)
ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)#, labelpad=10)
ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)

#ax.set_xlim(0, 450)
#ax.set_ylim(0, 23)
if test == 0:
  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)
plt.savefig(basename+".diff.pdf")

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)

ax.plot(filled_shear[1:], filled_stress[1:]/neat_stress[1:], marker=filled_marker, label="Total")
#ax.plot(filled_shear, filled_filler_stress , marker=filled_marker, label="Filler")
ax.plot(filled_shear[1:], filled_polymer_stress[1:]/neat_stress[1:], marker=filled_marker, label="Polymer")

#ax.set_xticks([0, 0.5, 1.0, 1.5, 2])
ax.tick_params(axis='both', labelsize=labelsize, pad=0)
ax.set_ylabel(r"$\sigma_{\mathrm{Filled}}/\sigma_{\mathrm{Neat}}$", fontsize=fontsize)#, labelpad=10)
ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)#, labelpad=10)
ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)

#ax.set_xlim(0, 450)
#ax.set_ylim(0, 23)
if test == 0:
  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)
plt.savefig(basename+".ratio.pdf")
