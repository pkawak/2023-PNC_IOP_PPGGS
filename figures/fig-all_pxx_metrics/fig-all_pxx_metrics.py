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
#matplotlib.use('tkagg')
import pandas as pd

# figure settings
test = 1
fontsize  = 9
labelsize = 8
figsize   = (2.65, 2.6)
bin       = 1
num_bins  = 200
pt_to_in  = 72
figsize   = (180/pt_to_in, 190/pt_to_in)
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

neat_base_df     = neat_base_df.mean()
filled_base_df   = filled_base_df.mean()

if bin:
  bins = np.linspace(np.amin(neat_df["Time"]), np.amax(neat_df["Time"])+1, int(len(neat_df["Time"])/num_bins))
  neat_df["Time_Binned"] = pd.cut(neat_df["Time"], bins, right=False)
  neat_df = neat_df.groupby("Time_Binned", as_index=False).mean()

  bins = np.linspace(np.amin(filled_df["Time"]), np.amax(filled_df["Time"])+1, int(len(filled_df["Time"])/num_bins))
  filled_df["Time_Binned"] = pd.cut(filled_df["Time"], bins, right=False)
  filled_df = filled_df.groupby("Time_Binned", as_index=False).mean()

neat_shear    = np.append(0, neat_df["v_shear_shift"])
filled_shear  = np.append(0, filled_df["v_shear_shift"])
neat_stress   = np.append(neat_base_df["v_pxx"], neat_df["v_pxx"])
filled_stress = np.append(filled_base_df["v_pxx"], filled_df["v_pxx"])

filled_polymer_pi = np.append(filled_base_df["c_pxx_polymer[1]"]/Npolymer, filled_df["c_pxx_polymer[1]"]/Npolymer)
filled_filler_pi  = np.append(filled_base_df["c_pxx_filler[1]"]/Nfiller, filled_df["c_pxx_filler[1]"]/Nfiller)

filled_polymer_stress = np.append(filled_base_df["c_pxx_polymer[1]"]/Ntotal, filled_df["c_pxx_polymer[1]"]/Ntotal)
filled_filler_stress  = np.append(filled_base_df["c_pxx_filler[1]"]/Ntotal, filled_df["c_pxx_filler[1]"]/Ntotal)

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(6, 10), constrained_layout=test)

ax = axs[0]
ax.plot(neat_shear  , neat_stress  , marker=neat_marker, alpha=0.3, label=r"$\sigma_{0}$")
ax.plot(filled_shear, filled_stress        , marker=filled_marker, color='b', label=r"$\sigma$")
ax.plot(filled_shear, filled_filler_pi , marker='s', color='c', label=r"$\pi_{f}$")
ax.plot(filled_shear, filled_polymer_pi, marker='s', color='r', label=r"$\pi_{p}$")
ax.plot(filled_shear, filled_filler_stress , marker=filled_marker, color='c', label=r"$\sigma_{f}$")
ax.plot(filled_shear, filled_polymer_stress, marker=filled_marker, color='r', label=r"$\sigma_{p}$")
ax.plot(filled_shear, (filled_filler_stress+filled_polymer_stress), ls='', marker='^', color='b', label=r"$\sigma_{p}+\sigma_{f}$")
ax.set_ylabel(r"Stress, $\sigma$", fontsize=fontsize)

ax = axs[1]
ax.plot(filled_shear, filled_stress-neat_stress, marker=filled_marker, label=r"$\Delta \sigma$")
#ax.plot(filled_shear, filled_filler_stress , marker=filled_marker, label=r"$\sigma_{f}$")
#ax.plot(filled_shear, filled_polymer_stress-neat_stress, marker=filled_marker, label=r"$\Delta \sigma_{p} = \sigma_{p} - \sigma_{0}$")
ax.plot(filled_shear, filled_polymer_pi-neat_stress, marker='^', color='r', label=r"$\Delta \bar{\sigma_{p}} = \bar{\sigma_{p}} - \sigma_{0}$")
ax.plot(filled_shear, filled_polymer_pi, marker=filled_marker, color='r', label=r"$\bar{\sigma_{p}}$")
ax.plot(filled_shear, filled_filler_pi, marker=filled_marker, color='b', label=r"$\bar{\sigma_{f}}$")
ax.plot(filled_shear, (filled_filler_pi - filled_polymer_pi) * Nfiller/Ntotal, marker='^', color='b', label=r"$(\bar{\sigma_{f}} - \bar{\sigma_{p}})\frac{N_{\mathrm{filler}}}{N}$")
ax.set_ylabel(r"$\sigma_{\mathrm{Filled}} - \sigma_{\mathrm{Neat}}$", fontsize=fontsize)

ax = axs[2]
ax.plot(filled_shear, filled_polymer_pi-neat_stress, marker=filled_marker, color='b', label=r"How filler modifies polymer response $\Delta \pi_{p} = \pi_{p} - \sigma_{0}$")
ax.plot(filled_shear, (filled_filler_pi - filled_polymer_pi) * Nfiller/Ntotal, marker=filled_marker, color='r', label=r"Direct excess contribution of filler $(\pi_{f} - \pi_{p})\frac{N_{\mathrm{filler}}}{N}$")
ax.set_ylabel(r"$\sigma_{\mathrm{Filled}} - \sigma_{\mathrm{Neat}}$", fontsize=fontsize)

for ax in axs.flat:
  ax.tick_params(axis='both', labelsize=labelsize, pad=0)
  ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)
  ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)

if test == 0:
  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)

import os
basename = os.path.splitext(__file__)[0]
plt.savefig(basename+".pdf")
