#!/usr/bin/env python3
# coding: utf-8
# 

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
    
    parser.add_argument("drug_data_folder", action="store", help="folder were all the folders with drug data for one cell line is stored (e.g. output) ")

    parser.add_argument("--drug_name", action="store", help="name of the drug that should be tested")

    parser.add_argument("--drug_level", nargs='?', action="store", help="IC level that the boxplots should be created from (e.g. IC90)")
    
    parser.add_argument("--wildtype_folder", nargs='?')
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

    simulation_paths = []

    for simulation in os.listdir(args.drug_data_folder):
        if args.drug_name in simulation:
            if args.drug_level:
                if simulation[17] == args.drug_level:
                    simulation_paths.append(args.drug_data_folder + simulation)
            else:
                simulation_paths.append(args.drug_data_folder + simulation)
    for directory in os.listdir(args.wildtype_folder):
        print(directory)
        if os.path.isdir(os.path.join(args.wildtype_folder, directory)):
            simulation_paths.append(args.wildtype_folder + directory)


    for el in simulation_paths:
        print(el)

    sim_auc_live_pairs = {}
    sim_auc_apoptotic_pairs = {}
    drug_1 = []
    drug_2 = []
    conc_1 = []
    conc_2 = []
    
    for data_folder in simulation_paths:
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
        
        # calculate the auc for the current simulation 
        auc_live = np.trapz(df_time_course['live'], df_time_course.time)
        auc_live_sims = integrate.simps(df_time_course['live'], df_time_course.time)
        print("AUC trapz: " + str(auc_live))
        print("AUC simps: " + str(auc_live_sims))
        sim_auc_live_pairs[simulation_name] = auc_live

        auc_apoptotic =  np.trapz(df_time_course['apoptotic'], df_time_course.time)
        sim_auc_apoptotic_pairs[simulation_name] = auc_apoptotic


    drug_dataframe = pd.DataFrame(list(sim_auc_live_pairs.items()), columns=["simulation", "auc_live"]) 
    drug_dataframe["auc_apoptotic"] = list(sim_auc_apoptotic_pairs.values())
    drug_dataframe["drug_1"] = drug_1
    drug_dataframe["drug_2"] = drug_2
    drug_dataframe["conc_1"] = conc_1
    drug_dataframe["conc_2"] = conc_2

    print(drug_dataframe)

    #####################
    # create the boxplots
    #####################

    # the x-axis describes the drug concentration while the y-axis describes AUC 
    sns.set(style="whitegrid")
    ax = sns.boxplot(x='conc_1', y='auc_live', data=drug_dataframe, order=["None", "IC10", "IC30", "IC50", "IC70", "IC90"])
    ax.set_xlabel("Drug concentration")
    ax.set_ylabel("AUC")
    plt.savefig('boxplot_' + args.drug_name + ".png")

main()     