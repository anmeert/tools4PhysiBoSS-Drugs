#!/usr/bin/env python3
# coding: utf-8

# based on plot_time_course.py from tools4physicell github repo

import os, re, sys
import glob, json

import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy.io import loadmat

import matplotlib.pyplot as plt
import seaborn as sns
import statistics


modules_path = os.path.dirname(os.path.realpath(__file__))
modules_path = os.path.join(modules_path, 'modules')
sys.path.append(modules_path)

import multicellds 

def get_time_for_cell_limit(timecourse_df, cell_limit):
    # get all times 
    time_list = timecourse_df.time.unique()
    for index, time in enumerate(time_list):
        time_df = timecourse_df[timecourse_df["time"] == time]
        if (sum(time_df["live"]) > cell_limit):
            return time_list[index-1]
    return 0

def simtype(row):
    if "output" in row["simulation_name"]:
        return "no-drug"
    else:
        return "drug"

def main():
    np.set_printoptions(threshold=sys.maxsize)

    data_path = "output/time_course_cell_layers_50.csv"
    wildtype_name = "output_LNCaP"
    drug_simulation_name = "LNCaP_mut_RNA_00_3_5_Ipatasertib_Pictilisib_0_0_0_0"
    cell_limit = 1700
    
    cell_layer_df = pd.read_csv(data_path)
    cell_layer_df["sim_type"] = cell_layer_df.apply (lambda row: simtype(row), axis=1)
    grouped_cell_layer_df = cell_layer_df.groupby(["cell_layer","sim_type", "time"])["live"].median().reset_index()
    ending_time = get_time_for_cell_limit(grouped_cell_layer_df[grouped_cell_layer_df["sim_type"] == "no-drug"], cell_limit)
    print(ending_time)
    num_cell_layers = 5
    # calculate the median wildtype AUC for each cell layer
    wildtype_median_AUCs = {"1": 0, "2": 0, "3": 0, "4":0, "5": 0}
    wildtype_df = cell_layer_df.loc[cell_layer_df["simulation_name"].str.contains(wildtype_name)]
    for cell_layer_num in range(1,6):
        replicate_aucs = []
        for replicate_num in range(1,11):
            wildtype_layer_df = wildtype_df[(wildtype_df["cell_layer"] == cell_layer_num) & (wildtype_df["simulation_name"].str.endswith(str(replicate_num)))]
            wildtype_layer_df_time_cut = wildtype_layer_df[wildtype_layer_df["time"] <= ending_time]
            auc = np.trapz(wildtype_layer_df_time_cut["live"].values, wildtype_layer_df_time_cut["time"].values)
            replicate_aucs.append(auc) 
        wildtype_median_AUCs[str(cell_layer_num)] = statistics.median(replicate_aucs)
    print(wildtype_median_AUCs)

    # calculate the median drug simulation AUC for each cell layer
    drug_median_AUCs = {"1": 0, "2": 0, "3": 0, "4":0, "5": 0}
    drug_df = cell_layer_df.loc[cell_layer_df["simulation_name"].str.contains(drug_simulation_name)]
    for cell_layer_num in range(1,6):
        replicate_aucs = []
        for replicate_num in range(1,11):
            drug_layer_df = drug_df[(drug_df["cell_layer"] == cell_layer_num) & (drug_df["simulation_name"].str.endswith(str(replicate_num)))]
            drug_layer_df_time_cut = drug_layer_df[drug_layer_df["time"] <= ending_time]
            auc = np.trapz(drug_layer_df_time_cut["live"].values, drug_layer_df_time_cut["time"].values)
            replicate_aucs.append(auc)  
        drug_median_AUCs[str(cell_layer_num)] = statistics.median(replicate_aucs)
    print(drug_median_AUCs)

    # calculate the growth index for each layer
    growth_index_layers = {"1": 0, "2": 0, "3": 0, "4":0, "5": 0}
    for layer_num in range(1,6):
        growth_index = np.log2(drug_median_AUCs[str(layer_num)] / wildtype_median_AUCs[str(layer_num)])
        growth_index_layers[str(layer_num)] = growth_index
    print(growth_index_layers)


    # plot the growth curves for each layer
   
    # Set time column as the dataframe index
    sns.set_context('paper')
    patch_color = "lightgrey"
    
    print("Creating figure")
    # Create a figure
    fig, axs = plt.subplots(1, 5, figsize=(8,3), dpi=300)
    
    # plot Alive vs Time
    curve_params = {}
    # cellline_color_dict = {'BPH1': '#f50707', 'DU145':'#f58a07', '22Rv1': '#0b03ff', 'PC3': '#12de19', 'LNCaP': '#0ffaf2', 'VCaP': '#a210eb'}
    cellline_color_dict = {'output': '#f50707', 'Pictilisib':'#0b03ff'}
    curve_params['live'] = {'color': '#75db75', 'label': 'Alive'}
    curve_params['apoptotic'] = {'color': '#ef4242', 'label': 'Apoptotic'}
    curve_params['necrotic'] = {'color':'#97723d', 'label': 'Necrotic'}
    line_width = 1.

    for cell_layer in range(1,6):
        drug_cell_layer_df = cell_layer_df[(cell_layer_df["cell_layer"] == cell_layer) & (cell_layer_df["simulation_name"].str.contains("Ipatasertib_Pictilisib"))]
        drug_cell_layer_df = drug_cell_layer_df.groupby(['time'])['live'].median().reset_index()
        print(drug_cell_layer_df)
        drug_cell_layer_df.to_csv("test.csv")
        axs[cell_layer-1].plot(drug_cell_layer_df.time, drug_cell_layer_df["live"], "-", c='#0b03ff', label="Ipatasertib: IC50, Pictilisib: IC90", linewidth=line_width)
        no_drug_cell_layer_df = cell_layer_df[(cell_layer_df["cell_layer"] == cell_layer) & (cell_layer_df["simulation_name"].str.contains("output"))]
        no_drug_cell_layer_df = no_drug_cell_layer_df.groupby(['time'])['live'].median().reset_index()
        axs[cell_layer-1].plot(no_drug_cell_layer_df.time, no_drug_cell_layer_df["live"], "-", c= '#f50707', label="no drug", linewidth=line_width)
        axs[cell_layer-1].set(ylim=(0,7000))
        axs[cell_layer-1].title.set_text("Cell layer " + str(cell_layer))
    # setting axes labels
    handles, labels = axs[4].get_legend_handles_labels()
    fig.legend(handles, labels, loc = 'lower right')
    axs[0].set_ylabel('Nº of cells')
    axs[2].set_xlabel('Time (min)')
    
    
    # axs[cell_layer-1].tick_params(axis='x', labelsize=12)
    # axs[cell_layer-1].tick_params(axis='y', labelsize=12)

    
    
    # Saving fig
    fig.tight_layout(pad=2.4, w_pad=0.5, h_pad=1.4)
    # different figure name for each data folder 
    base_name = "SuppMat_Fig12"
    fig_fname = "output/cell_layer_timecourse.png"
    if not os.path.exists('output'):
        os.makedirs('output')
    fig.savefig(fig_fname)
    print("Saving fig as %s" % fig_fname)

main() 

# %%
