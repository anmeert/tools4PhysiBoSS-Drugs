#!/usr/bin/env python3
# coding: utf-8

# based on plot_time_course.py from tools4physicell github repo

import os, re, sys
import glob, json
import argparse

import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy.io import loadmat

import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations, product


modules_path = os.path.dirname(os.path.realpath(__file__))
modules_path = os.path.join(modules_path, 'modules')
sys.path.append(modules_path)

import multicellds 


sns.set(style="ticks", palette="Paired")
sns.set_context('paper')

def create_parser():
    parser = argparse.ArgumentParser(description="Plot total cell grouped as Alive/Necrotic/Apoptotic vs Time")
    
    parser.add_argument("cellline_data_folder", action="store", help="folder were the folders with drug data for all cell lines are stored (e.g. output) ")

    return parser
    

def main():
   
    parser = create_parser()
    args = parser.parse_args()
    
    phases_dict = multicellds.default_phases_dict
    phase_grouping = multicellds.default_phase_grouping
    
    celllines_sims_paths = []
    for directory in os.listdir(args.cellline_data_folder):
        if "output_" in directory:
            print(directory)
            celllines_sims_paths.append(args.cellline_data_folder + directory)

    for el in celllines_sims_paths:
        print(el)
    # Set time column as the dataframe index
    sns.set_context('paper')
    patch_color = "lightgrey"
    
    print("Creating figure")
    # Create a figure
    fig, ax = plt.subplots(1, 1, figsize=(8,3), dpi=300)
    
    # plot Alive vs Time
    curve_params = {}
    cellline_color_dict = {'BPH1': '#f50707', 'DU145':'#f58a07', '22Rv1': '#0b03ff', 'PC3': '#12de19', 'LNCaP': '#0ffaf2', 'VCaP': '#a210eb'}
    curve_params['live'] = {'color': '#75db75', 'label': 'Alive'}
    curve_params['apoptotic'] = {'color': '#ef4242', 'label': 'Apoptotic'}
    curve_params['necrotic'] = {'color':'#97723d', 'label': 'Necrotic'}
    line_width = 1.

    for data_folder in celllines_sims_paths:
        # Globing output files according to the output format specified
        phase_col = "current_phase"

        mcds = multicellds.MultiCellDS(output_folder=data_folder)
        df_iterator = mcds.cells_as_frames_iterator()
        num_of_files = mcds.cells_file_count()
        
        # Initializing a Pandas Databrafe to store the data
        columns = ["time", "live", "apoptotic", "necrotic"]
        data = np.zeros((num_of_files, 4), dtype=int)
        df_time_course = pd.DataFrame(columns=columns, data=data)

        print("Reading cell_output files from %i input files from %s" % (num_of_files, data_folder))
        # Iterating over all cell_output files
        for i, (t, df) in enumerate(df_iterator):
            print("\tProcessing time step: %.0f" % t)

            # Rename the phases integer codes using the phases_dict as the mapping
            s = df[phase_col]
            s.replace(to_replace=phases_dict, value=None, inplace=True)
            
            # Count the number of cells in each phase
            counts = s.value_counts()
        
            df_time_course.loc[i, 'time'] = t
            print("time:", t)
            # group the previous phases count into the three general classes:
            # Alive, Apoptotic, Necrotic
            for k, v in counts.to_dict().items():
                if k not in phase_grouping:
                    continue
                df_time_course.loc[i, phase_grouping[k]] += v
                
        for k,v in cellline_color_dict.items():
            if k in data_folder:
                c = v
                l = k
        ax.plot(df_time_course.time, df_time_course["live"], "-", c=c, label=l, linewidth=line_width)
        
    # setting axes labels
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Nº of cells")

    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

    ax.set(ylim=(0,50000))
    
    # Showing legend
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.grid(linestyle='dotted')
    ax.legend()
    
    # Saving fig
    fig.tight_layout()
    # different figure name for each data folder 
    base_name = "all_celllines"
    fig_fname = "./timecourse" + "_" + base_name + ".png"
    fig.savefig(fig_fname)
    print("Saving fig as %s" % fig_fname)
    
main() 

# %%
