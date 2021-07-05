#!/usr/bin/env python3
# coding: utf-8


import os, re, sys
import glob, json

import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy.io import loadmat

import matplotlib.pyplot as plt
import seaborn as sns


modules_path = os.path.dirname(os.path.realpath(__file__))
modules_path = os.path.join(modules_path, 'modules')
sys.path.append(modules_path)

import multicellds 


def main():
    wildtype_folder = ".../PhysiBoSSv2/output_LNCaP/"   
    drug_simulation_folder = ".../PhysiBoSSv2/output/double/"
    drug_simulation_name = "LNCaP_mut_RNA_00_3_5_Ipatasertib_Pictilisib_0_0_0_0"
    
    wildtype_paths = []
    for folder in os.listdir(wildtype_folder):
        wildtype_paths.append(wildtype_folder + folder)
        
    drug_simulation_paths = []
    for folder in os.listdir(drug_simulation_folder):
        if drug_simulation_name in folder:
            drug_simulation_paths.append(drug_simulation_folder + folder)
        
    

    layer_dataframe = pd.DataFrame()

    for data_folder in wildtype_paths + drug_simulation_paths:

        mcds = multicellds.MultiCellDS(output_folder=data_folder)
        df_iterator = mcds.cells_as_frames_iterator()
        num_of_files = mcds.cells_file_count()

        phases_dict = multicellds.default_phases_dict
        phase_grouping = multicellds.default_phase_grouping

        # Initializing a Pandas Dataframe to store the data
        columns = ["time", "live", "apoptotic", "necrotic", "cell_layer"]
        data = np.zeros((num_of_files*5, 5), dtype=int)
        df_time_course = pd.DataFrame(columns=columns,data= data)
   
        print("Reading cell_output files from %i input files from %s" % (num_of_files, data_folder))
        # Iterating over all cell_output files
        for i, (t, df) in enumerate(df_iterator):
            print("\tProcessing time step: %.0f" % t)
            
            # separate the dataframe in different cell layers according to the coordinates
            
            center = np.array((0,0,0))
            cell_layers = {'1': pd.Interval(left=0, right=50), '2': pd.Interval(left=50, right=100), '3': pd.Interval(left=100, right= 150), '4': pd.Interval(left=150, right=200), '5': pd.Interval(left=200, right=250)}

            df["eucl_distance"] = df.apply(lambda row: np.linalg.norm(np.array((row["x_position"], row["y_position"], row["z_position"]) - center)), axis=1)
        

            # assign a layer to every cell
            
            def get_cell_layer (row):
                for k,v in cell_layers.items():
                    if row["eucl_distance"] in v:
                        return k
    
            df["cell_layer"] = df.apply(lambda row: get_cell_layer(row), axis=1)
            # count number of cells in every layer

            # Rename the phases integer codes using the phases_dict as the mapping
           
            # Count the number of cells in each phase
            for layer,interval in cell_layers.items():
                current_index = len(cell_layers)*i + int(layer) -1
                sub_df = df[df["cell_layer"] == layer]
                phase_column = sub_df["current_phase"]
                phase_column.replace(to_replace=phases_dict, value=None, inplace=True)
                if len(phase_column) != 0:
                    counts = phase_column.value_counts()
                    for phase, phase_count in counts.to_dict().items():
                        if phase not in phase_grouping:
                            continue
                        df_time_course.loc[current_index, phase_grouping[phase]] += phase_count
                df_time_course.loc[current_index, "cell_layer"] = layer
                df_time_course.loc[current_index, "time"] = t
            

            print("time:", t)
        df_time_course["simulation_name"] = os.path.basename(os.path.normpath(data_folder))
        frames = [layer_dataframe, df_time_course]
        layer_dataframe = pd.concat(frames)
    layer_dataframe.to_csv("output/time_course_cell_layers_50.csv")
    
main() 

# %%
