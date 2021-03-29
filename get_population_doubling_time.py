#!/usr/bin/env python3
# coding: utf-8

import os, re, sys
import glob, json
import argparse

import numpy as np
import pandas as pd
from scipy.io import loadmat

import matplotlib.pyplot as plt
import seaborn as sns


modules_path = os.path.dirname(os.path.realpath(__file__))
modules_path = os.path.join(modules_path, 'modules')
sys.path.append(modules_path)

import multicellds 


sns.set(style="ticks", palette="Paired")
sns.set_context('paper')

def read_csv():
    # reading a cell_output file (plain text ; separated columns)
    # any function can be used here, using pandas is just a shorcut
    if args.format == 'mat':
        df = process_mat(f)
        phase_col = "current_phase"
        # This should be changed regardin time stamp in MultiCellDS
        time *= 60 
        
def create_parser():
    parser = argparse.ArgumentParser(description="Plot total cell grouped as Alive/Necrotic/Apoptotic vs Time")
    
    parser.add_argument("data_folder", action="store", help="folder were the data is stored")

    parser.add_argument("--figout", action="store", dest="fig_fname", default="./cell_vs_time.png",
                        help="File name to save the plot")

    parser.add_argument("--motout", action="store", dest="mot_fname", default="./mot_vs_time.png",
                    help="File name to save the plot")
    return parser

def pb_output_iterator(data_folder, sep=";"):
    globing = os.path.join(data_folder, "cells_[0-9]*.txt")
    for fname in sorted(glob.glob(globing)):
        df = pd.read_csv(fname, sep=sep)
        t = df.Time[0]
        yield (t, df)

def count_pb_files(data_folder):
    globing = os.path.join(data_folder, "cells_[0-9]*.txt")
    return len(glob.glob(globing))

def main():
   
    parser = create_parser()
    args = parser.parse_args()
    
    phases_dict = multicellds.default_phases_dict
    phase_grouping = multicellds.default_phase_grouping
    print(phase_grouping)

    # Globing output files according to the output format specified
    
    phase_col = "current_phase"
    mcds = multicellds.MultiCellDS(output_folder=args.data_folder)
    df_iterator = mcds.cells_as_frames_iterator()
    num_of_files = mcds.cells_file_count()

    # Initializing a Pandas Databrafe to store the data
    columns = ["time", "live", "apoptotic", "necrotic"]
    data = np.zeros((num_of_files, 4), dtype=int)
    df_time_course = pd.DataFrame(columns=columns, data=data)

    print("Reading cell_output files from %i input files from %s" % (num_of_files, args.data_folder))
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
    
    # get the two data-points in which the population doubling happens 
    # intial number of cells = 1138, so doubling = 2276

    # get the time where live cells are 
    print(df_time_course)
    doubling_index = 0
    for index, row in df_time_course.iterrows():
        if df_time_course.loc[index, "live"] > 2276:
            break
        doubling_index += 1
    
    # retrieve sub dataframe with just 
    doubling_df = df_time_course.iloc[[doubling_index-1, doubling_index], :]
    cellnum_diff = doubling_df.at[doubling_index, "live"] - doubling_df.at[doubling_index-1, "live"]
    time_diff = doubling_df.at[doubling_index, "time"] - doubling_df.at[doubling_index-1, "time"]
    slope = cellnum_diff / time_diff
    # b = y - ax
    b = doubling_df.at[doubling_index, "live"] - slope * doubling_df.at[doubling_index, "time"]
    # calculate doubling time
    # x = (y - b)/a
    double_numcell = 2276
    doubling_min = (double_numcell - b) / slope
    doubling_hours = doubling_min / 60
    print(doubling_min)
    print(doubling_hours)
main()
    
   