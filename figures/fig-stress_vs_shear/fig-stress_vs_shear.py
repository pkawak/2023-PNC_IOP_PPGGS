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
num_bins  = 200
pt_to_in  = 72
filter_min_shear = 0.0
filter_shear = 500
neat_marker = "^"
filled_marker = "o"
neat_color = "r"
filler_color = "b"
filled_color = "#808080"
polymer_color = "r"

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
neat_df       ["pnorm"] = neat_df["v_pyy"] + neat_df["v_pzz"]
filled_df     ["pnorm"] = filled_df["v_pyy"] + filled_df["v_pzz"]
neat_base_df  ["pnorm"] = neat_base_df["v_pyy"] + neat_base_df["v_pzz"]
filled_base_df["pnorm"] = filled_base_df["v_pyy"] + filled_base_df["v_pzz"]

neat_df          = neat_df[neat_df["v_shear_shift"] < filter_shear]
filled_df        = filled_df[filled_df["v_shear_shift"] < filter_shear]
neat_base_df     = neat_base_df[neat_base_df["v_shear_shift"] < filter_shear]
filled_base_df   = filled_base_df[filled_base_df["v_shear_shift"] < filter_shear]
neat_df          = neat_df[neat_df["v_shear_shift"] > filter_min_shear]
filled_df        = filled_df[filled_df["v_shear_shift"] > filter_min_shear]
neat_base_df     = neat_base_df[neat_base_df["v_shear_shift"] > filter_min_shear]
filled_base_df   = filled_base_df[filled_base_df["v_shear_shift"] > filter_min_shear]

neat_base_df     = neat_base_df.mean()
filled_base_df   = filled_base_df.mean()

if bin:
  bins = np.linspace(np.amin(neat_df["Time"]), np.amax(neat_df["Time"])+1, int(len(neat_df["Time"])/num_bins))
  neat_df["Time_Binned"] = pd.cut(neat_df["Time"], bins, right=False)
  neat_df = neat_df.groupby("Time_Binned", as_index=False).mean()

  bins = np.linspace(np.amin(filled_df["Time"]), np.amax(filled_df["Time"])+1, int(len(filled_df["Time"])/num_bins))
  filled_df["Time_Binned"] = pd.cut(filled_df["Time"], bins, right=False)
  filled_df = filled_df.groupby("Time_Binned", as_index=False).mean()

# gotta add the quiescent at 0
neat_shear    = np.append(0, neat_df["v_shear_shift"])
filled_shear  = np.append(0, filled_df["v_shear_shift"])
neat_stress   = np.append(neat_base_df["v_pxx"], neat_df["v_pxx"])
filled_stress = np.append(filled_base_df["v_pxx"], filled_df["v_pxx"])
neat_stress_norm   = -np.append(neat_base_df["pnorm"], neat_df["pnorm"])
filled_stress_norm = -np.append(filled_base_df["pnorm"], filled_df["pnorm"])
neat_stress_norm   = -neat_df["pnorm"]
filled_stress_norm = -filled_df["pnorm"]

filled_polymer_pi = np.append(filled_base_df["c_pxx_polymer[1]"]/Npolymer, filled_df["c_pxx_polymer[1]"]/Npolymer)
filled_filler_pi  = np.append(filled_base_df["c_pxx_filler[1]"]/Nfiller, filled_df["c_pxx_filler[1]"]/Nfiller)

filled_polymer_pi_norm = -np.append(filled_base_df["c_pxx_polymer[2]"]/Npolymer + filled_base_df["c_pxx_polymer[3]"]/Npolymer, filled_df["c_pxx_polymer[2]"]/Npolymer + filled_df["c_pxx_polymer[3]"]/Npolymer)
filled_filler_pi_norm  = -np.append(filled_base_df["c_pxx_filler[2]"]/Nfiller   + filled_base_df["c_pxx_filler[3]"]/Nfiller  , filled_df["c_pxx_filler[2]"]/Nfiller   + filled_df["c_pxx_filler[3]"]/Nfiller  )
filled_polymer_pi_norm = -filled_df["c_pxx_polymer[2]"]/Npolymer - filled_df["c_pxx_polymer[3]"]/Npolymer
filled_filler_pi_norm  = -filled_df["c_pxx_filler[2]"]/Nfiller   - filled_df["c_pxx_filler[3]"]/Nfiller

filled_polymer_stress = np.append(filled_base_df["c_pxx_polymer[1]"]/Ntotal, filled_df["c_pxx_polymer[1]"]/Ntotal)
filled_filler_stress  = np.append(filled_base_df["c_pxx_filler[1]"]/Ntotal, filled_df["c_pxx_filler[1]"]/Ntotal)

filled_polymer_stress_norm = -np.append(filled_base_df["c_pxx_polymer[2]"]/Ntotal + filled_base_df["c_pxx_polymer[3]"]/Ntotal, filled_df["c_pxx_polymer[2]"]/Ntotal + filled_df["c_pxx_polymer[3]"]/Ntotal)
filled_filler_stress_norm  = -np.append(filled_base_df["c_pxx_filler[2]"]/Ntotal   + filled_base_df["c_pxx_filler[3]"]/Ntotal  , filled_df["c_pxx_filler[2]"]/Ntotal   + filled_df["c_pxx_filler[3]"]/Ntotal  )
filled_polymer_stress_norm = -filled_df["c_pxx_polymer[2]"]/Ntotal - filled_df["c_pxx_polymer[3]"]/Ntotal
filled_filler_stress_norm  = -filled_df["c_pxx_filler[2]"]/Ntotal  - filled_df["c_pxx_filler[3]"]/Ntotal

import os
basename = os.path.splitext(__file__)[0]

def common_format(ax):
  ax.set_ylabel(r"Stress, $\sigma$", fontsize=fontsize)
  #ax.set_ylabel(r"Stress, $\sigma_{\mathrm{extensional}}$", fontsize=fontsize)
  ax.tick_params(axis='both', labelsize=labelsize, pad=0)
  ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)
  ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)
  if test == 0:
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
ax.plot(neat_shear  , neat_stress  , marker=neat_marker  , color=neat_color, label="Neat")
ax.plot(filled_shear, filled_stress, marker=filled_marker, color=filled_color, label="Filled")
ax.set_xticks([0, 0.5, 1.0, 1.5, 2])
common_format(ax)
plt.savefig(basename+".pdf")

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
#ax.plot(neat_shear  , neat_stress  , marker=neat_marker  , color=neat_color  , alpha=0.3, label="Neat")
ax.plot(filled_shear, filled_stress, marker=filled_marker, color=filled_color, label="Total")
ax.plot(filled_shear, filled_filler_stress , marker=filled_marker, color=filler_color, label="Filler")
ax.plot(filled_shear, filled_polymer_stress, marker=filled_marker, color=polymer_color, label="Polymer")
common_format(ax)
plt.savefig(basename+".groups.pdf")

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
ax.plot(filled_shear, filled_stress-neat_stress, marker=filled_marker, label="Total")
#ax.plot(filled_shear, filled_filler_stress , marker=filled_marker, label="Filler")
#ax.plot(filled_shear, filled_polymer_stress-neat_stress, marker=filled_marker, label="Polymer")
#common_format(ax)
ax.set_ylabel(r"$\sigma_{\mathrm{Filled}} - \sigma_{\mathrm{Neat}}$", fontsize=fontsize)#, labelpad=10)
ax.tick_params(axis='both', labelsize=labelsize, pad=0)
ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)
plt.savefig(basename+".diff.pdf")

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
ax.plot(neat_shear[1:]  , neat_stress_norm  , marker=neat_marker  , color=neat_color, label="Neat")
ax.plot(filled_shear[1:], filled_stress_norm, marker=filled_marker, color=filled_color, label="Filled")
common_format(ax)
ax.set_ylabel(r"$P_{\mathrm{normal}}$", fontsize=fontsize)
plt.savefig(basename+".norm.pdf")

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
#ax.plot(neat_shear  , neat_stress_norm  , marker=neat_marker  , color=neat_color  , alpha=0.3, label="Neat")
ax.plot(filled_shear[1:], filled_stress_norm, marker=filled_marker, color=filled_color, label="Total")
ax.plot(filled_shear[1:], filled_filler_stress_norm , marker=filled_marker, color=filler_color, label="Filler")
ax.plot(filled_shear[1:], filled_polymer_stress_norm, marker=filled_marker, color=polymer_color, label="Polymer")
common_format(ax)
ax.set_ylabel(r"$P_{\mathrm{normal}}$", fontsize=fontsize)
plt.savefig(basename+".norm.groups.pdf")

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
ax.plot(filled_shear[1:], filled_stress_norm-neat_stress_norm, marker=filled_marker, label="Total")
ax.plot(filled_shear[1:], filled_filler_stress_norm , marker=filled_marker, label="Filler")
ax.plot(filled_shear[1:], filled_polymer_stress_norm-neat_stress_norm, marker=filled_marker, label="Polymer")
common_format(ax)
ax.set_ylabel(r"$\sigma_{\mathrm{Filled}} - \sigma_{\mathrm{Neat}}$", fontsize=fontsize)#, labelpad=10)
plt.savefig(basename+".norm.diff.pdf")

figsize   = (150/pt_to_in, 110/pt_to_in)
fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
ax.plot(filled_shear[1:], filled_stress[1:]/neat_stress[1:], marker=filled_marker, label="Total")
#ax.plot(filled_shear, filled_filler_stress , marker=filled_marker, label="Filler")
#ax.plot(filled_shear[1:], filled_polymer_stress[1:]/neat_stress[1:], marker=filled_marker, label="Polymer")
ax.plot([0, 2], [1, 1], ls='--', color='k')

#ax.set_xticks([0, 0.5, 1.0, 1.5, 2])
ax.tick_params(axis='both', labelsize=labelsize, pad=0)
ax.set_ylabel(r"$\sigma_{\mathrm{Filled}}/\sigma_{\mathrm{Neat}}$", fontsize=fontsize)#, labelpad=10)
ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)#, labelpad=10)
#ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)

#ax.set_xlim(0, 450)
#ax.set_ylim(0, 23)
if test == 0:
  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)
plt.savefig(basename+".ratio.pdf")

figsize   = (180/pt_to_in, 190/pt_to_in)
fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
ax.plot(neat_shear  , neat_stress  , marker=neat_marker  , color=neat_color, label="Neat")
ax.plot(filled_shear, filled_stress, marker=filled_marker, color=filled_color, label="Filled")
ax.set_xticks([0, 0.5, 1.0, 1.5, 2])
common_format(ax)
plt.savefig(basename+".small.pdf")

figsize   = (150/pt_to_in, 110/pt_to_in)
fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
ax.plot(neat_shear  , neat_stress  , marker=neat_marker  , color=neat_color, label="Neat")
ax.plot(filled_shear, filled_stress, marker=filled_marker, color=filled_color, label="Filled")
ax.set_xticks([0, 0.5, 1.0, 1.5, 2])
common_format(ax)
plt.savefig(basename+".smaller.pdf")

