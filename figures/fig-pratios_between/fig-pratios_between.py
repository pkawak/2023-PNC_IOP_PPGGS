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
from scipy import stats
from scipy.integrate import quad

# figure settings
test = 1
fontsize  = 9
labelsize = 8
figsize   = (2.65, 2.6)
bin       = 0
bin2      = 1
num_bins  = 10
pt_to_in  = 72
filter_min_shear = 0.45
filter_shear = 500
figsize   = (150/pt_to_in, 110/pt_to_in)
#139.42093pt
neat_marker = "^"
filled_marker = "o"
bw_marker     = "s"

color_polymer = "r"
color_filler  = "#808080"
color_all     = "c"

# smarter binning algorithm
shear_delta = 0.1

# constants
Npolymer = 102500
Nfiller  = 51450
Ntotal   = Npolymer+Nfiller

analyses_glob = ["all", "polymer", "filler", "polymer_full", "filler_full",
"bw_fillers_4.0_0.0", "bw_fillers_4.5_0.0", "bw_fillers_5.0_0.0", "bw_fillers_5.5_0.0",
"bw_fillers_4.0_0.0_1", "bw_fillers_4.5_0.0_1", "bw_fillers_5.0_0.0_1", "bw_fillers_5.5_0.0_1",
"bw_fillers_wxlink_4.0_0.0_0", "bw_fillers_wxlink_4.5_0.0_0", "bw_fillers_wxlink_5.0_0.0_0", "bw_fillers_wxlink_5.5_0.0_0",
"bw_shells_3.5_0.0_1", "bw_shells_3.75_0.0_1", "bw_shells_4.0_0.0_1", "bw_shells_4.25_0.0_1"]
names = ["all", "polymer", "filler", "polymer_full", "filler_full",
"polymer_4.0", "polymer_4.5", "polymer_5.0", "polymer_5.5",
"polymer_4.0_diff", "polymer_4.5_diff", "polymer_5.0_diff", "polymer_5.5_diff",
"polymer_4.0_wxlink", "polymer_4.5_wxlink", "polymer_5.0_wxlink", "polymer_5.5_wxlink",
"shells_4.0_diff", "shells_4.5_diff", "shells_5.0_diff", "shells_5.5_diff"]

ke       = ["ke"]
pe       = ["pe"]
press    = ["pxx", "pyy", "pzz"]
stress   = ["pxy", "pxz", "pyz"]
values   = ["ke", "pe", "pxx", "pyy", "pzz", "pxy", "pxz", "pyz"]
#print(analyses_glob)

def get_file_name(anal, val):
  return("../particles_between_data_files/" + anal+"."+val+".stats")

def get_df(anal, val):
  file_name = get_file_name(anal, val)
  df = pd.read_csv(file_name, delimiter=' ', skiprows=1, header=0, index_col=False)
#  df["pnorm"] = df["pyy_Mean"] + df["pzz_Mean"]
  return(df)
#print(dfs[0])
#sys.exit()

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

def get_all_df(anal):
  df_ = []
  for val, i in zip(values, np.arange(len(values))):
    df = get_df(anal, val)
    df["Time"] = np.arange(len(df))
    df.columns = [ val+"_"+col if col!="Time" and col!="Count" else col for col in df.columns ]
    if i == 0:
      df_ = df
    else:
      df = df.drop(columns=["Count"])
      df_ = pd.merge(df_, df, on="Time")#left_index=True, right_index=True)
  df_["press_Mean"]    = (df_["pxx_Mean"]+df_["pyy_Mean"]+df_["pzz_Mean"])/3.0
  df_["press_StdDev"]  = (df_["pxx_StdDev"]+df_["pyy_StdDev"]+df_["pzz_StdDev"])/3.0
  df_["stress_Mean"]   = -(df_["pxy_Mean"]+df_["pxz_Mean"]+df_["pyz_Mean"])/3.0
  df_["stress_StdDev"] = -(df_["pxy_StdDev"]+df_["pxz_StdDev"]+df_["pyz_StdDev"])/3.0
  df_["Time"]          = filled_df["Time"]
  df_["v_shear_shift"] = filled_df["v_shear_shift"]
  df_["pnorm"]         = df_["pyy_Mean"] + df_["pzz_Mean"]

  df_          = df_[df_["v_shear_shift"] < filter_shear]
  df_          = df_[df_["v_shear_shift"] > filter_min_shear]

  if bin:
    bins = np.linspace(np.amin(df_["Time"]), np.amax(df_["Time"])+1, int(len(df_["Time"])/num_bins))
    df_["Time_Binned"] = pd.cut(df_["Time"], bins, right=False)
    df_ = df_.groupby("Time_Binned", as_index=False).mean()

  return(df_) 

dfs = []
for anal in analyses_glob:
  dfs.append(get_all_df(anal))

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
neat_stress_norm   = np.append(neat_base_df["pnorm"], neat_df["pnorm"])
filled_stress_norm = np.append(filled_base_df["pnorm"], filled_df["pnorm"])

filled_polymer_stress = np.append(filled_base_df["c_pxx_polymer[1]"]/Npolymer, filled_df["c_pxx_polymer[1]"]/Npolymer)
filled_filler_stress  = np.append(filled_base_df["c_pxx_filler[1]"]/Nfiller, filled_df["c_pxx_filler[1]"]/Nfiller)

filled_polymer_stress_norm = np.append(filled_base_df["c_pxx_polymer[2]"]/Npolymer + filled_base_df["c_pxx_polymer[3]"]/Npolymer, filled_df["c_pxx_polymer[2]"]/Npolymer + filled_df["c_pxx_polymer[3]"]/Npolymer)
filled_filler_stress_norm  = np.append(filled_base_df["c_pxx_filler[2]"]/Nfiller   + filled_base_df["c_pxx_filler[3]"]/Nfiller  , filled_df["c_pxx_filler[2]"]/Nfiller   + filled_df["c_pxx_filler[3]"]/Nfiller  )

fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)

polymer_ratio = dfs[analyses_glob.index("bw_fillers_wxlink_4.0_0.0_0")]
polymer_ratio["ratio_norm"] = polymer_ratio["pnorm"]/filled_polymer_stress_norm[1:]
polymer_ratio["ratio"] = polymer_ratio["pxx_Mean"]/filled_polymer_stress[1:]
filler_ratio  = dfs[analyses_glob.index("bw_shells_4.0_0.0_1")]
filler_ratio["ratio_norm"] = filler_ratio["pnorm"]/filled_filler_stress_norm[1:]
filler_ratio["ratio"] = filler_ratio["pxx_Mean"]/filled_filler_stress[1:]
if bin2:
  bins = np.linspace(np.amin(filled_df["Time"]), np.amax(filled_df["Time"])+1, int(len(filled_df["Time"])/num_bins))
  polymer_ratio["Time_Binned"] = pd.cut(polymer_ratio["Time"], bins, right=False)
  filler_ratio["Time_Binned"] = pd.cut(filler_ratio["Time"], bins, right=False)
  polymer_ratio = polymer_ratio.groupby("Time_Binned", as_index=False).mean()
  filler_ratio = filler_ratio.groupby("Time_Binned", as_index=False).mean()

# fit to line
x_fit = np.linspace(filter_min_shear, 2, 10) 
slope, intercept, r_value, p_value, std_err = stats.linregress(filler_ratio["v_shear_shift"], filler_ratio["ratio"])
ratio_fit_line = x_fit*slope + intercept
model          = np.poly1d(np.polyfit(filler_ratio["v_shear_shift"], filler_ratio["ratio"], 2))
ratio_fit_quad = model(x_fit)
slope, intercept, r_value, p_value, std_err = stats.linregress(filler_ratio["v_shear_shift"], filler_ratio["ratio_norm"])
ratio_norm_fit_line = x_fit*slope + intercept
model               = np.poly1d(np.polyfit(filler_ratio["v_shear_shift"], filler_ratio["ratio_norm"], 2))
ratio_norm_fit_quad = model(x_fit)

#ax.plot(filler_ratio["v_shear_shift"], filler_ratio["ratio"]     , marker=bw_marker, color=color_polymer, alpha=0.1, label=r"$\sigma_{\mathrm{extension}}$")
#ax.plot(filler_ratio["v_shear_shift"], filler_ratio["ratio_norm"], marker=bw_marker, color=color_filler , alpha=0.1, label=r"$\sigma_{\mathrm{normal}}$")

ax.plot(x_fit, ratio_fit_line     , color=color_polymer, label=r"$\sigma_{\mathrm{extension}}$")
ax.plot(x_fit, ratio_norm_fit_line, color=color_filler , label=r"$\sigma_{\mathrm{normal}}$")
#ax.plot(x_fit, ratio_fit_quad     , color=color_polymer, label=r"$\sigma_{\mathrm{extension}}$")
#ax.plot(x_fit, ratio_norm_fit_quad, color=color_filler , label=r"$\sigma_{\mathrm{normal}}$")

ax.tick_params(axis='both', labelsize=labelsize, pad=0)
ax.set_ylabel(r"$\left[ \sigma_{\mathrm{Friction}}/\sigma \right]_{\mathrm{Filler}}$", fontsize=fontsize)#, labelpad=10)
ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)#, labelpad=10)
ax.legend(loc='best', fontsize=labelsize, handlelength=1.0)

#ax.set_xlim(0, 450)
#ax.set_ylim(0, 23)
if test == 0:
  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)
import os
basename = os.path.splitext(__file__)[0]
plt.savefig(basename+".ratio.pdf")
#
#fig, ax = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
#
#polymer_ratio = dfs[analyses_glob.index("bw_fillers_wxlink_4.0_0.0_0")]
#polymer_ratio["ratio_norm"] = (polymer_ratio["pnorm"]-filled_polymer_stress_norm[1:])/filled_polymer_stress_norm[1:]
#polymer_ratio["ratio"] = (polymer_ratio["pxx_Mean"]-filled_polymer_stress[1:])/filled_polymer_stress[1:]
#filler_ratio  = dfs[analyses_glob.index("bw_shells_4.0_0.0_1")]
#filler_ratio["ratio_norm"] = (filler_ratio["pnorm"]-filled_filler_stress_norm[1:])/filled_filler_stress_norm[1:]
#filler_ratio["ratio"] = (filler_ratio["pxx_Mean"]-filled_filler_stress[1:])/filled_filler_stress[1:]
#if bin2:
#  bins = np.linspace(np.amin(filled_df["Time"]), np.amax(filled_df["Time"])+1, int(len(filled_df["Time"])/num_bins))
#  polymer_ratio["Time_Binned"] = pd.cut(polymer_ratio["Time"], bins, right=False)
#  filler_ratio["Time_Binned"] = pd.cut(filler_ratio["Time"], bins, right=False)
#  polymer_ratio = polymer_ratio.groupby("Time_Binned", as_index=False).mean()
#  filler_ratio = filler_ratio.groupby("Time_Binned", as_index=False).mean()
#
#ax.plot(filler_ratio["v_shear_shift"], filler_ratio["ratio"]     , marker=bw_marker, color=color_polymer, label=r"$\sigma_{\mathrm{extension}}$")
#ax.plot(filler_ratio["v_shear_shift"], filler_ratio["ratio_norm"], marker=bw_marker, color=color_filler , label=r"$\sigma_{\mathrm{normal}}$")
#
#ax.tick_params(axis='both', labelsize=labelsize, pad=0)
##ax.set_ylabel(r"$\Delta\sigma/\sigma_{\mathrm{All}}$", fontsize=fontsize)
#ax.set_ylabel(r"$\left[ \Delta \sigma_{\mathrm{Friction}}/\sigma \right]_{\mathrm{Filler}}$", fontsize=fontsize)
#ax.set_xlabel(r"Strain, $\gamma$", fontsize=fontsize)#, labelpad=10)
#ax.legend(loc='best', fontsize=labelsize, handlelength=1.5)
#
##ax.set_xlim(0, 450)
##ax.set_ylim(0, 23)
#if test == 0:
#  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)
#plt.savefig(basename+".diff_ratio.pdf")
