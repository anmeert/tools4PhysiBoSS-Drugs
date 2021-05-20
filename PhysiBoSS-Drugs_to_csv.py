#!/usr/bin/env python3
# coding: utf-8

# based on plot_time_course.py from tools4physicell github repo
# author Annika Meert 

import os, re, sys
import argparse

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from itertools import compress


modules_path = os.path.dirname(os.path.realpath(__file__))
modules_path = os.path.join(modules_path, 'modules')
sys.path.append(modules_path)

import multicellds 

def create_parser():
    parser = argparse.ArgumentParser(description="Plot total cell grouped as Alive/Necrotic/Apoptotic vs Time")
    
    parser.add_argument("-double","--double_folder", nargs='?', action="store", help="folder were all the folders with double drug data for one cell line is stored (e.g. double) ")
    parser.add_argument("-single","--single_folder", nargs='?',  action="store", help="folder were all the folders with single drug data for one cell line is stored (e.g. single) ")
    parser.add_argument("-untreated","--untreated_folder", nargs='?', action="store", help="folder that contains all replicates of the untreated simulation (e.g. output_LNCaP")
    
    return parser
    


def drug_level_to_IC(drug_level):
    if (drug_level == "1"):
        return "IC10"
    elif (drug_level == "2"):
        return "IC30"
    elif (drug_level == "3"):
        return "IC50"
    elif (drug_level == "4"):
        return "IC70"
    elif (drug_level == "5"):
        return "IC90"
    else:
        return ""


def main():
   
    parser = create_parser()
    args = parser.parse_args()
    
    phases_dict = multicellds.default_phases_dict
    phase_grouping = multicellds.default_phase_grouping
    
    single_sims_paths = []
    double_sims_paths = []
    untreated_sims_paths = []
    dir_list = [args.double_folder,args.single_folder] 
    
    if args.double_folder:
        subfolders = os.listdir(args.double_folder)
        for folder in subfolders:
            double_sims_paths.append(args.double_folder + folder)
    
    if args.single_folder:
        subfolders = os.listdir(args.single_folder)
        for folder in subfolders:
            single_sims_paths.append(args.single_folder + folder)
    
    if args.untreated_folder:
        subfolders = os.listdir(args.untreated_folder)
        for folder in subfolders:
            if os.path.isdir(args.untreated_folder + folder):
                untreated_sims_paths.append(args.untreated_folder + folder)
    
    all_sims_paths = single_sims_paths + double_sims_paths + untreated_sims_paths
   
    # for el in all_sims_paths:
    #     print(el)

    simulation_names = []
    prolif_auc = []
    apop_auc = []
    drug_1 = []
    drug_2 = []
    conc_1 = []
    conc_2 = []

    for data_folder in all_sims_paths:
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
                
        # get the drug names and the drug concentrations for the current simulation 
        drug_names = ["Ipatasertib", "Afatinib", "Ulixertinib", "Luminespib", "Selumetinib", "Pictilisib"]
        simulation_name = os.path.basename(os.path.normpath(data_folder))
        # has to be modified so it's in the right order
        used_drugs = [drug for drug in drug_names if drug in simulation_name]
        first_drug = second_drug = ""
        first_conc = second_conc = ""
        if len(used_drugs) == 0:
            #wildtype 
            first_drug = second_drug = "wildtype"
            first_conc = second_conc = "None"
        elif len(used_drugs) == 1:
            first_drug = used_drugs.pop(0)
            # replace drug levels with IC values 
            first_conc = drug_level_to_IC(simulation_name[17])
            second_drug = first_drug
            second_conc = first_conc
        else:
            first_drug = used_drugs.pop(0)
            first_conc = drug_level_to_IC(simulation_name[17])
            second_drug = used_drugs.pop(0)
            second_conc  = drug_level_to_IC(simulation_name[19])
               
        drug_1.append(first_drug)
        drug_2.append(second_drug)
        conc_1.append(first_conc)
        conc_2.append(second_conc)
        simulation_names.append(simulation_name)

        # calculate the auc's for the current simulation 
        auc_live = np.trapz(df_time_course['live'], df_time_course.time)
        prolif_auc.append(auc_live)
        auc_apoptotic =  np.trapz(df_time_course['apoptotic'], df_time_course.time)
        apop_auc.append(auc_apoptotic)

    # store the information in a dataframe
    drug_dataframe = pd.DataFrame({'simulation_name': simulation_name, 'drug_1': drug_1, 'drug_2': drug_2, 'conc_1': conc_1, 'conc_2': conc_2, 'auc_live': auc_live, 'auc_apop': auc_apoptotic})

    output_name = "simulation_data.csv"
    if not os.path.exists('output'):
        os.makedirs('output')
    drug_dataframe.to_csv('output/' + output_name) 
main()