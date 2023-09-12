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

# figure settings
test = 1
fontsize         = 9
labelsize        = 8
figsize          = (2.65, 2.6)
pt_to_in         = 72
neat_marker      = "^"
filled_marker    = "o"
bw_marker        = "s"
filler_color     = "b"
filled_color     = "#808080"
polymer_color    = "r"

# binning settings
binIt            = 1
num_bins         = 20
bin_by           = "Time"
binIt2           = 0
num_bins2        = 20

# filter data settings
filter_min       = 0.0
filter_max       = 0.35
filterIt         = 0
filter_by        = "v_shear_shift"

data_dir         = "../../../Draft-PNC_Letter/data/"

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
filled_log_file_names = [ "deform_" + filled_file_name.split('_')[-1].split('-data')[0] + ".log" for filled_file_name in filled_file_names ]
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
name_neat_dfs     = np.append(["Avg"], [ neat_file_name.split('_')[2].split('-data')[0] for neat_file_name in neat_file_names ])
neat_log_file_names = [ "neat-deform_" + neat_file_name.split('_')[2].split('-data')[0] + ".log" for neat_file_name in neat_file_names ]
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

rho = Ntotal/Vtotals[0]
rho0 = Npolymer/neat_Vtotals[0]
delrho = rho - rho0
phi_F = Vfiller/Vtotals[0]

analyses_glob = ["all", "polymer", "filler",
"bw_fillers_wxlink_4.0_0.0_0", "bw_fillers_wxlink_4.5_0.0_0",
"bw_shells_3.5_0.0_1", "bw_shells_3.75_0.0_1"]
names = ["all", "polymer", "filler",
"polymer_4.0_wxlink", "polymer_4.5_wxlink",
"shells_3.5_diff", "shells_3.75_diff"]

ke       = ["ke"]
pe       = ["pe"]
press    = ["pxx", "pyy", "pzz"]
stress   = ["pxy", "pxz", "pyz"]
values   = ["ke", "pe", "pxx", "pyy", "pzz", "pxy", "pxz", "pyz"]

def get_file_name(anal, val, repnum=1):
  return(data_dir + "/particles_between_data_files/rep" + str(repnum) + "." + anal+"."+val+".stats")

def get_df(anal, val, repnum=1):
  file_name = get_file_name(anal, val, repnum)
  df = pd.read_csv(file_name, delimiter=' ', skiprows=1, header=0, index_col=False)
  return(df)

def get_all_df(anal, repnum=1):
  df_ = []
  for val, i in zip(values, np.arange(len(values))):
    df = get_df(anal, val, repnum)
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
  df_["pnorm_Mean"]    = - (df_["pyy_Mean"]   + df_["pzz_Mean"]  )/2
  df_["pnorm_StdDev"]  = - (df_["pyy_StdDev"] + df_["pzz_StdDev"])/2

  if binIt:
    bins = np.linspace(np.amin(df_[bin_by]), np.amax(df_[bin_by])+(df_[bin_by].iat[2]-df_[bin_by].iat[1]), num_bins+1)
    bins = np.linspace(np.amin(neat_df[bin_by]), np.amax(neat_df[bin_by])+(neat_df[bin_by].iat[2]-neat_df[bin_by].iat[1]), num_bins+1)
    df_["Time_Binned"] = pd.cut(df_[bin_by], bins, right=False)
    df_ = df_.groupby("Time_Binned", as_index=False).mean()

  return(df_) 

dfs_reps = [ [ get_all_df(anal, repnum) for anal in analyses_glob ] for repnum in range(1,6) ]
dfs = dfs_reps[0]
num_reps = len(dfs_reps)

def compute_stats_dfs(dfs_reps):
  dfs_mean = []
  dfs_std  = []
  for anal in analyses_glob:
    df_subset = [ dfs_reps[repnum][analyses_glob.index(anal)] for repnum in range(num_reps) ]
    df_concat = pd.concat(df_subset)
    by_row_index          = df_concat.groupby(df_concat.index)
    dfs_mean.append(by_row_index.mean())
    dfs_std.append(by_row_index.std())
  
  return dfs_mean, dfs_std

dfs_mean, dfs_std = compute_stats_dfs(dfs_reps)
#print("dfs:")
#print([ dfs_reps[repnum][analyses_glob.index("polymer")]["pyz_Mean"].to_string(index=False) for repnum in range(len(dfs_reps)) ])
#print("dfs_mean:")
#print(dfs_mean[analyses_glob.index("polymer")]["pyz_Mean"].to_string(index=False))
#print(dfs_std[analyses_glob.index("polymer")]["pyz_Mean"].to_string(index=False))

filled_result_df = filled_result_dfs[0]
neat_result_df = neat_result_dfs[0]
alpha = 1.0
plotted_df    = pd.DataFrame()
plotted_df["shear"] = filled_result_df["shear"]

def plot_err(ax, y, marker="None", color=filled_color, alpha=alpha, label='none', x=filled_result_df["shear"], markerfacecolor='none', name="", ls="-", is_add=1):
  ymean = np.mean(y, axis=0)
  ystd  = np.std(y, axis=0)
  yse   = ystd/np.sqrt(len(y))
  ax.errorbar(x, ymean, yerr=yse, marker=marker, color=color, markerfacecolor=markerfacecolor, alpha=alpha, ls=ls, label=label)
  if is_add:
    plotted_df[name] = ymean
    plotted_df[name+"_se"] = yse

# other types of plots follow
def plot_it(polymer_ratio, filler_ratio, ylabel, out_file, xname = "v_shear_shift", y1name = "ratio", y2name = "ratio_norm", ylabel_norm='', plot_base=0):

  fig1, ax1 = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
  fig2, ax2 = plt.subplots(1, 1, sharex=True, figsize=figsize, constrained_layout=test)
  ls = '-'
  if ylabel_norm=='':
    ylabel_norm = ylabel
  
  ax = ax1
  if plot_base:
    y = [ (yi["filler_pi"]).tolist() for yi in filled_result_dfs_reps ]
    plot_err(ax, y, color=filler_color, label="Filler")
    y = [ (yi["polymer_pi"]).tolist() for yi in filled_result_dfs_reps ]
    plot_err(ax, y, color=polymer_color, label="Polymer")
    ls = ''
  y = [ (yi[y1name]).tolist() for yi in filler_ratio ]
  plot_err(ax, y, x=filler_ratio[0][xname], color=filler_color , marker=bw_marker, ls=ls, is_add=0, label="Filler Contact")
  y = [ (yi[y1name]).tolist() for yi in polymer_ratio ]
  plot_err(ax, y, x=filler_ratio[0][xname], color=polymer_color, marker=bw_marker, ls=ls, is_add=0, label="Interfiller Polymer")
#  ax.plot(filler_ratio[xname] , filler_ratio[y1name] , ls=ls, marker=bw_marker, color=filler_color, label="Filler Contact")
#  ax.plot(polymer_ratio[xname], polymer_ratio[y1name], ls=ls, marker=bw_marker, color=polymer_color, label="Interfiller Polymer")

  ax.tick_params(axis='both', labelsize=labelsize, pad=0)
  ax.set_ylabel(ylabel, fontsize=fontsize)#, labelpad=10)
  ax.legend(loc='best', fontsize=labelsize, handlelength=1.5)
  
  ax = ax2
  if plot_base:
    y = [ (yi["filler_pi_norm"]).tolist() for yi in filled_result_dfs_reps ]
    plot_err(ax, y, color=filler_color, label="Filler")
    y = [ (yi["polymer_pi_norm"]).tolist() for yi in filled_result_dfs_reps ]
    plot_err(ax, y, color=polymer_color, label="Polymer")
  y = [ (yi[y2name]).tolist() for yi in filler_ratio ]
  plot_err(ax, y, x=filler_ratio[0][xname], color=filler_color , marker=bw_marker, ls=ls, is_add=0, label="Filler Contact")
  y = [ (yi[y2name]).tolist() for yi in polymer_ratio ]
  plot_err(ax, y, x=filler_ratio[0][xname], color=polymer_color, marker=bw_marker, ls=ls, is_add=0, label="Interfiller Polymer")
  
  ax.tick_params(axis='both', labelsize=labelsize, pad=0)
  ax.set_ylabel(ylabel_norm, fontsize=fontsize)#, labelpad=10)
  ax.set_xlabel(r"Extensional Strain, $\varepsilon$", fontsize=fontsize)#, labelpad=10)
#  ax.legend(loc='best', fontsize=labelsize, handlelength=1.5)

  if test == 0:
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.999, top=0.999)
  return fig1, ax1, fig2, ax2

# process particles_between data
polymer_ratio = [ dfs[analyses_glob.index("bw_fillers_wxlink_4.0_0.0_0")] for dfs in dfs_reps ]
filler_ratio  = [ dfs[analyses_glob.index("bw_shells_3.75_0.0_1")] for dfs in dfs_reps ]
for repnum in range(num_reps):
  polymer_ratio[repnum]["ratio"]      = polymer_ratio[repnum]["pxx_Mean"]
  filler_ratio [repnum]["ratio"]      = filler_ratio [repnum]["pxx_Mean"]
  polymer_ratio[repnum]["ratio_norm"] = polymer_ratio[repnum]["pnorm_Mean"]
  filler_ratio [repnum]["ratio_norm"] = filler_ratio [repnum]["pnorm_Mean"]
ylabel      = r"$\langle S_{xx}^{i} \rangle_{i}$"
ylabel_norm = r"$\left( \langle S_{yy}^{i} \rangle_{i} + \langle S_{zz}^{i} \rangle_{i} \right)/2$"
ylabel      = r"Per Bead Ext. Stress"
ylabel_norm = r"Per Bead Normal Pressure"

basename = os.path.splitext(__file__)[0]
out_file1 = basename+".x.pdf"
out_file2 = basename+".y.pdf"
fig1, ax1, fig2, ax2 = plot_it(polymer_ratio, filler_ratio, ylabel, out_file1, ylabel_norm=ylabel_norm, plot_base=1)
ax2.set_yticks([-0.05, 0.0, 0.05, 0.1, 0.15, 0.2])
ax2.set_yticklabels(["", "0.0", "", "0.1", "", "0.2"])
print("Saving ", out_file1)
fig1.savefig(out_file1)
print("Saving ", out_file2)
fig2.savefig(out_file2)
plt.close()
