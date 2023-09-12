#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
created on Fri March 3 05:48:00 2023

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
figsize   = (2.65, 2.6)
pt_to_in         = 72
neat_marker      = "^"
filled_marker    = "o"
neat_color       = "r"
filler_color     = "b"
filled_color     = "#808080"
polymer_color    = "r"
other_color      = "c"

# binning settings
binIt            = 1
num_bins         = 20
bin_by           = "v_shear_shift"

# filter data settings
filter_min       = 0.0
filter_max       = 0.35
filterIt         = 0
filter_by        = "v_shear_shift"

data_dir         = "../../data/"

# read params
Ntotal        = int(os.popen("grep '[0-9] atoms' " + data_dir + "/deform_1.log | sed 's/  //g' | head -n1 | cut -d' ' -f1").read())
Npolymer      = int(os.popen("grep '[0-9] atoms' " + data_dir + "/deform_1.log | sed 's/  //g' | grep 'polymer' | cut -d' ' -f1").read())
Nfiller       = Ntotal - Npolymer

filled_file_names     = os.listdir(data_dir)
filled_file_names     = [ data_dir + f for f in filled_file_names if f.startswith("deform") and "dataframe" in f ]
filled_base_file_name = data_dir+"/quiescent-dataframe.out"
neat_file_names       = os.listdir(data_dir)
neat_file_names       = [ data_dir + f for f in neat_file_names if f.startswith("neat-deform") and "dataframe" in f ]
neat_base_file_name   = data_dir+"/neat-quiescent-dataframe.out"

neat_dfs              = [ pd.read_csv(neat_file_name) for neat_file_name in neat_file_names ]
neat_df_concat        = pd.concat(neat_dfs)
neat_by_row_index     = neat_df_concat.groupby(neat_df_concat.index)
neat_df               = neat_by_row_index.mean()
#neat_df.to_csv("neat_dataframe.out")

filled_dfs            = [ pd.read_csv(filled_file_name) for filled_file_name in filled_file_names ]
df_concat             = pd.concat(filled_dfs)
by_row_index          = df_concat.groupby(df_concat.index)
filled_df             = by_row_index.mean()
#filled_df.to_csv("dataframe.out")

neat_base_df          = pd.read_csv(neat_base_file_name)
filled_base_df        = pd.read_csv(filled_base_file_name)

# copied in from voronoi calc in ovito for out.npt.data in rep1
# TODO have to fix these to also be for each run
Vpolymer           = 112731.026
Vfiller            = 39519.715

def process_df(df, base_df, Lx0, Ntotal=Ntotal):
  Vtotal = Lx0**3
  df     ["pnorm"] = (df["v_pyy"] + df["v_pzz"])/2
  base_df["pnorm"] = (base_df["v_pyy"] + base_df["v_pzz"])/2

  df     ["tshear"] = np.log(1+df["v_shear_shift"])
  if filterIt:
    df        = df[df[filter_by] < filter_max]
    df        = df[df[filter_by] > filter_min]

  if binIt:
    bins = np.linspace(np.amin(df[bin_by]), np.amax(df[bin_by])+(df[bin_by].iat[2]-df[bin_by].iat[1]), num_bins+1)
    df["Time_Binned"] = pd.cut(df[bin_by], bins, right=False)
    df = df.groupby("Time_Binned", as_index=False).mean()

  df2 = pd.Series({"v_shear_shift": 0})
  df1 = df["v_shear_shift"]
  df_result = df1.append(df2, ignore_index=True)
  df_result = pd.concat([df2, df1]).reset_index(drop=True)
  df_result = df_result.to_frame()
  df_result.columns = ["shear"]
  df_result["stress"]      = np.append(base_df["v_pxx"], df["v_pxx"])*Ntotal/Vtotal
  df_result["stress_norm"] = -np.append(base_df["pnorm"], df["pnorm"])*Ntotal/Vtotal
  df_result["pi"]          = np.append(base_df["v_pxx"], df["v_pxx"])
  df_result["pi_norm"]     = -np.append(base_df["pnorm"], df["pnorm"])
  df_result["tshear"]      = np.log(1+df_result["shear"])
  df_result["temp"]        = np.append(base_df["c_ke_total"]  , df["c_ke_total"])/1.5/Ntotal

  df_result["polymer_stress_x_phi"]  = np.append(base_df["c_pxx_polymer[1]"], df["c_pxx_polymer[1]"])/Vtotal
  df_result["polymer_stress_x_phi_norm"] = -np.append(base_df["c_pxx_polymer[2]"] + base_df["c_pxx_polymer[3]"], df["c_pxx_polymer[2]"] + df["c_pxx_polymer[3]"])/2/Vtotal
  df_result["polymer_pi"]  = np.append(base_df["c_pxx_polymer[1]"], df["c_pxx_polymer[1]"])/Npolymer
  df_result["polymer_pi_norm"] = -np.append(base_df["c_pxx_polymer[2]"] + base_df["c_pxx_polymer[3]"], df["c_pxx_polymer[2]"] + df["c_pxx_polymer[3]"])/2/Npolymer
  df_result["polymer_stress"] = np.append(base_df["c_pxx_polymer[1]"], df["c_pxx_polymer[1]"])/Ntotal
  df_result["polymer_stress_norm"] = -np.append(base_df["c_pxx_polymer[2]"] + base_df["c_pxx_polymer[3]"], df["c_pxx_polymer[2]"] + df["c_pxx_polymer[3]"])/2/Ntotal
  df_result["polymer_stress_vol"] = np.append(base_df["c_pxx_polymer[1]"], df["c_pxx_polymer[1]"])/Vpolymer
  df_result["polymer_stress_vol_norm"] = -np.append(base_df["c_pxx_polymer[2]"] + base_df["c_pxx_polymer[3]"], df["c_pxx_polymer[2]"] + df["c_pxx_polymer[3]"])/2/Vpolymer
  df_result["polymer_temp"] = np.append(base_df["c_ke_polymer"], df["c_ke_polymer"])/1.5/Npolymer

  if "c_pxx_filler[1]" in df:
    df_result["filler_stress_x_phi"]   = np.append(base_df["c_pxx_filler[1]"] , df["c_pxx_filler[1]"] )/Vtotal
    df_result["filler_stress_x_phi_norm"]  = -np.append(base_df["c_pxx_filler[2]"]  + base_df["c_pxx_filler[3]"] , df["c_pxx_filler[2]"]  + df["c_pxx_filler[3]"] )/2/Vtotal
    df_result["filler_pi"]   = np.append(base_df["c_pxx_filler[1]"] , df["c_pxx_filler[1]"] )/Nfiller
    df_result["filler_pi_norm"]  = -np.append(base_df["c_pxx_filler[2]"]  + base_df["c_pxx_filler[3]"] , df["c_pxx_filler[2]"]  + df["c_pxx_filler[3]"] )/2/Nfiller
    df_result["filler_stress"]  = np.append(base_df["c_pxx_filler[1]"] , df["c_pxx_filler[1]"] )/Ntotal
    df_result["filler_stress_norm"]  = -np.append(base_df["c_pxx_filler[2]"]  + base_df["c_pxx_filler[3]"] , df["c_pxx_filler[2]"]  + df["c_pxx_filler[3]"] )/2/Ntotal
    df_result["filler_stress_vol"]  = np.append(base_df["c_pxx_filler[1]"] , df["c_pxx_filler[1]"] )/Vfiller
    df_result["filler_stress_vol_norm"]  = -np.append(base_df["c_pxx_filler[2]"]  + base_df["c_pxx_filler[3]"] , df["c_pxx_filler[2]"]  + df["c_pxx_filler[3]"] )/2/Vfiller
    df_result["filler_temp"]  = np.append(base_df["c_ke_filler"] , df["c_ke_filler"])/1.5/Nfiller

  return(df_result)

neat_base_df  ["pnorm"] = (neat_base_df["v_pyy"] + neat_base_df["v_pzz"])/2
filled_base_df["pnorm"] = (filled_base_df["v_pyy"] + filled_base_df["v_pzz"])/2
neat_base_df     = neat_base_df.mean()
filled_base_df   = filled_base_df.mean()

filled_dfs        = [ pd.read_csv(filled_file_name) for filled_file_name in filled_file_names ]
filled_dfs        = [filled_df] + filled_dfs
name_filled_dfs   = np.append(["Avg"], [ filled_file_name.split('_')[1].split('-data')[0] for filled_file_name in filled_file_names ])
filled_log_file_names = [ "deform_" + filled_file_name.split('_')[1].split('-data')[0] + ".log" for filled_file_name in filled_file_names ]
Lxtotals          = [ int(os.popen("grep \"Lx\" " + data_dir + "/" + f + " | sed 's/\( \)\+/ /g' | sed 's/^ //g' | sed 's/ $//g' | tr -s ' ' '\n' | nl -nln | grep \"Lx\" | cut -f1").read()) for f in filled_log_file_names ]
Lxtotals          = [ float(os.popen("grep -A1 \"Lx\" " + data_dir + "/" + f + " | tail -n1 | sed 's/\( \)\+/ /g' | sed 's/^ //g' | sed 's/ $//g' | cut -d' ' -f" + str(Lxtotal)).read()) for Lxtotal, f in zip(Lxtotals, filled_log_file_names) ]
Lxtotals          = np.append(np.mean(Lxtotals), Lxtotals)
Vtotals           = Lxtotals**3
filled_result_dfs = [ process_df(filled_df_, filled_base_df, Lx0) for filled_df_, Lx0 in zip(filled_dfs, Lxtotals) ]

# get mean of filled_result_dfs but only reps
filled_result_dfs_reps = filled_result_dfs[1:]
df_concat             = pd.concat(filled_result_dfs_reps)
by_row_index          = df_concat.groupby(df_concat.index)
filled_result_df_mean = by_row_index.mean()
filled_result_df_std  = by_row_index.std()
#print(filled_result_df_mean)

neat_dfs          = [ pd.read_csv(neat_file_name) for neat_file_name in neat_file_names ]
neat_dfs          = [neat_df] + neat_dfs
name_neat_dfs     = np.append(["Avg"], [ neat_file_name.split('_')[1].split('-data')[0] for neat_file_name in neat_file_names ])
neat_log_file_names = [ "neat-deform_" + neat_file_name.split('_')[1].split('-data')[0] + ".log" for neat_file_name in neat_file_names ]
neat_Lxtotals     = [ int(os.popen("grep \"Lx\" " + data_dir + "/" + f + " | sed 's/\( \)\+/ /g' | sed 's/^ //g' | sed 's/ $//g' | tr -s ' ' '\n' | nl -nln | grep \"Lx\" | cut -f1").read()) for f in neat_log_file_names ]
neat_Lxtotals     = [ float(os.popen("grep -A1 \"Lx\" " + data_dir + "/" + f + " | tail -n1 | sed 's/\( \)\+/ /g' | sed 's/^ //g' | sed 's/ $//g' | cut -d' ' -f" + str(Lxtotal)).read()) for Lxtotal, f in zip(neat_Lxtotals, neat_log_file_names) ]
neat_Lxtotals     = np.append(np.mean(neat_Lxtotals), neat_Lxtotals)
neat_Vtotals      = neat_Lxtotals**3
neat_result_dfs   = [ process_df(neat_df_, neat_base_df, Lx0, Ntotal=Npolymer) for neat_df_, Lx0 in zip(neat_dfs, neat_Lxtotals) ]

# get mean of neat_result_dfs but only reps
neat_result_dfs_reps = neat_result_dfs[1:]
df_concat             = pd.concat(neat_result_dfs_reps)
by_row_index          = df_concat.groupby(df_concat.index)
neat_result_df_mean = by_row_index.mean()
neat_result_df_std  = by_row_index.std()

def get_diff(y):
  return(y-y[0])

alpha = 1.0
filled_result_df = filled_result_dfs[0]
neat_result_df = neat_result_dfs[0]

plotted_df    = pd.DataFrame()
plotted_df["shear"] = filled_result_df["shear"]
plotted_df["Neat_shear"] = neat_result_df["shear"]

def plot_err(ax, y, marker=filled_marker, color=filled_color, alpha=alpha, label='none', x=filled_result_df["shear"], markerfacecolor='none', name=""):
  ymean = np.mean(y, axis=0)
  ystd  = np.std(y, axis=0)
  yse   = ystd/np.sqrt(len(y))
  ax.errorbar(x, ymean, yerr=yse, marker=marker, color=color, markerfacecolor=markerfacecolor, alpha=alpha, label=label)
  plotted_df[name] = ymean
  plotted_df[name+"_se"] = yse

fig, ax = plt.subplots(1, figsize=figsize, sharex=True, constrained_layout=True)

# stress vs shear
y     = [ get_diff(yi["stress"]).tolist() for yi in neat_result_dfs_reps ]
plot_err(ax, y, label=r"Neat", marker=neat_marker, color=neat_color, alpha=alpha, x=neat_result_df["shear"], markerfacecolor=neat_color, name="Neat")
y     = [ get_diff(yi["stress"]).tolist() for yi in filled_result_dfs_reps ]
plot_err(ax, y, marker=filled_marker, markerfacecolor=filled_color , color=filled_color , alpha=alpha, label=r"Filled", name="Filled")
ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)
common_format(ax)
ax.set_ylabel(r"Extensional Stress, $\sigma$", fontsize=fontsize)
ax.set_xlabel(r"Extensional Strain, $\varepsilon$", fontsize=fontsize)

out_file_name = os.path.basename(sys.argv[0]).split('.')[0]+".pdf"
fig.savefig(out_file_name)

fig, ax = plt.subplots(1, figsize=figsize, sharex=True, constrained_layout=True)
y     = [ get_diff(yf["stress"]-yn["stress"]).tolist() for yf, yn in zip(filled_result_dfs_reps, neat_result_dfs_reps) ]
#y     = [ (get_diff(yf["stress"])-get_diff(yn["stress"])).tolist() for yf, yn in zip(filled_result_dfs_reps, neat_result_dfs_reps) ]
#y     = [ ((yf["stress"])-(yn["stress"])).tolist() for yf, yn in zip(filled_result_dfs_reps, neat_result_dfs_reps) ]
plot_err(ax, y, label=r"Neat", marker=neat_marker, color=neat_color, alpha=alpha, x=neat_result_df["shear"], markerfacecolor=neat_color, name="Neat")

common_format(ax)
ax.set_ylabel(r"$\sigma_{\mathrm{Filled}}-\sigma_{\mathrm{Neat}}$", fontsize=fontsize)
ax.set_xlabel(r"Extensional Strain, $\varepsilon$", fontsize=fontsize)

out_file_name = os.path.basename(sys.argv[0]).split('.')[0]+".diff.pdf"
fig.savefig(out_file_name)

fig, ax = plt.subplots(1, figsize=figsize, sharex=True, constrained_layout=True)
y     = [ get_diff(yi["stress"]).tolist() for yi in filled_result_dfs_reps ]
plot_err(ax, y, marker=filled_marker, markerfacecolor=filled_color , color=filled_color , alpha=alpha, label=r"$i=$Total", name="Total")
y     = [ get_diff(yi["polymer_stress_x_phi"]).tolist() for yi in filled_result_dfs_reps ]
plot_err(ax, y, marker=filled_marker, markerfacecolor='none'       , color=polymer_color , alpha=alpha, label=r"$i=$Polymer", name="Polymer")
y     = [ get_diff(yi["filler_stress_x_phi"]).tolist() for yi in filled_result_dfs_reps ]
plot_err(ax, y, marker=filled_marker, markerfacecolor='none'       , color=filler_color , alpha=alpha, label=r"$i=$Filler", name="Filler")

ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)
common_format(ax)
ax.set_ylabel(r"Stress $\times$ Volume Fraction, $\sigma_{i}*\phi_{i}$", fontsize=fontsize)
ax.set_xlabel(r"Extensional Strain, $\varepsilon$", fontsize=fontsize)

out_file_name = os.path.basename(sys.argv[0]).split('.')[0]+".groups.pdf"
fig.savefig(out_file_name)

fig, ax = plt.subplots(1, figsize=figsize, sharex=True, constrained_layout=True)
y     = [ get_diff(yi["stress_norm"]).tolist() for yi in filled_result_dfs_reps ]
plot_err(ax, y, marker=filled_marker, markerfacecolor=filled_color , color=filled_color , alpha=alpha, label=r"Total", name="Total")
y     = [ get_diff(yi["polymer_stress_norm"]).tolist() for yi in filled_result_dfs_reps ]
plot_err(ax, y, marker=filled_marker, markerfacecolor='none'       , color=polymer_color , alpha=alpha, label=r"Polymer", name="Polymer")
y     = [ get_diff(yi["filler_stress_norm"]).tolist() for yi in filled_result_dfs_reps ]
plot_err(ax, y, marker=filled_marker, markerfacecolor='none'       , color=filler_color , alpha=alpha, label=r"Filler", name="Filler")

#ax.legend(loc='best', fontsize=fontsize, handlelength=1.0)
common_format(ax)
ax.set_ylabel(r"$P_{\mathrm{normal}}$", fontsize=fontsize)
ax.set_xlabel(r"Extensional Strain, $\varepsilon$", fontsize=fontsize)

out_file_name = os.path.basename(sys.argv[0]).split('.')[0]+".norm.groups.pdf"
fig.savefig(out_file_name)

#y     = [ (phi_F * get_diff(yi["filler_stress_vol"] - yi["polymer_stress_vol"])).tolist() for yi in filled_result_dfs_reps ]
#plot_err(ax, y, filled_marker, filler_color , 1.0, r"$\phi_{F} \left( \sigma_{F} - \sigma_{P} \right)$", name="filler diff")
#
#y     = [ get_diff(yi["polymer_stress_vol"] - neat_result_df["stress"]).tolist() for yi in filled_result_dfs_reps ]
#plot_err(ax, y, filled_marker, polymer_color , 1.0, r"$\sigma_{P} - \sigma_{0}$", name="polymer_diff")

plotted_df.to_csv("plot.csv")
