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
    
    parser.add_argument("double_folder", action="store", help="folder were all the folders with double drug data for one cell line is stored (e.g. double) ")
    parser.add_argument("single_folder", action="store", help="folder were all the folders with single drug data for one cell line is stored (e.g. single) ")
    parser.add_argument("wildtype_folder", action="store", help="folder that contains all replicates of the wildtype simulation (e.g. output_LNCaP")
    
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
    dir_list = [args.double_folder,args.single_folder] 
    for directory in dir_list:
        all_subfolders = os.listdir(directory)
        for folder in all_subfolders:
            if "double" in directory:
                double_sims_paths.append(directory + folder)
            else:
                single_sims_paths.append(directory + folder)
    
    all_sims_paths = single_sims_paths + double_sims_paths
    for directory in os.listdir(args.wildtype_folder):
        if os.path.isdir(args.wildtype_folder + directory):
            all_sims_paths.append(args.wildtype_folder + directory)

    for el in all_sims_paths:
        print(el)


    sim_auc_live_pairs = {}
    sim_auc_apoptotic_pairs = {}
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
        
        # calculate the auc for the current simulation 
        auc_live = np.trapz(df_time_course['live'], df_time_course.time)
        auc_live_sims = integrate.simps(df_time_course['live'], df_time_course.time)
        print("AUC trapz: " + str(auc_live))
        print("AUC simps: " + str(auc_live_sims))
        sim_auc_live_pairs[simulation_name] = auc_live

        auc_apoptotic =  np.trapz(df_time_course['apoptotic'], df_time_course.time)
        sim_auc_apoptotic_pairs[simulation_name] = auc_apoptotic

        # Set time column as the dataframe index
        sns.set_context('paper')
        patch_color = "lightgrey"
        
        print("Creating figure")
        # Create a figure
        fig, ax = plt.subplots(1, 1, figsize=(8,3), dpi=300)
        
        # plot Alive vs Time
        curve_params = {}
        curve_params['live'] = {'color': '#75db75', 'label': 'Alive'}
        curve_params['apoptotic'] = {'color': '#ef4242', 'label': 'Apoptotic'}
        curve_params['necrotic'] = {'color':'#97723d', 'label': 'Necrotic'}
        line_width = 2.
        for k,pdict in curve_params.items():
            c = pdict['color']
            l = pdict['label']
            ax.plot(df_time_course.time, df_time_course[k], "-", c=c, label=l, linewidth=line_width)
        
        # setting axes labels
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Nº of cells")

        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)

        ax.set(ylim=(0,20000))
        
        # Showing legend
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.grid(linestyle='dotted')
        ax.legend()
        
        # Saving fig
        fig.tight_layout()
        # different figure name for each data folder 
        base_name = os.path.basename(os.path.normpath(data_folder))
        # fig_fname = "./timecourse" + "_" + base_name + ".png"
        # fig.savefig(fig_fname)
        # print("Saving fig as %s" % fig_fname)
    

    # with sim_auc pair list calculate now the list with differences in AUC to the wildtype 
    wildtype_name = "output_LNCaP"
    wildtype_auc_live_dict = {k:v for k,v in sim_auc_live_pairs.items() if wildtype_name in k}
    wildtype_auc_live = sum(wildtype_auc_live_dict.values()) / len(wildtype_auc_live_dict.values())
    auc_live_differences = {k: (v - wildtype_auc_live) for k,v in sim_auc_live_pairs.items()}
    log2_auc_live_ratio = {k: (np.log2(v/wildtype_auc_live)) for k,v in sim_auc_live_pairs.items()}

    # with sim_auc pair list calculate now the list with differences in AUC to the wildtype 
    wildtype_auc_apoptotic_dict = {k:v for k,v in sim_auc_apoptotic_pairs.items() if wildtype_name in k}
    wildtype_auc_apoptotic = sum(wildtype_auc_apoptotic_dict.values()) / len(wildtype_auc_apoptotic_dict.values())
    auc_apoptotic_differences = {k: (v - wildtype_auc_apoptotic) for k,v in sim_auc_apoptotic_pairs.items()}
    log2_auc_apoptotic_ratio = {k: (np.log2(v/wildtype_auc_apoptotic)) for k,v in sim_auc_apoptotic_pairs.items()} 

    # store the information in a dataframe
    drug_dataframe = pd.DataFrame(list(sim_auc_live_pairs.items()), columns=["simulation", "auc_live"]) 
    drug_dataframe["auc_live_difference"] = list(auc_live_differences.values())
    drug_dataframe["log2_auc_live_ratio"] = list(log2_auc_live_ratio.values())
    drug_dataframe["auc_apoptotic"] = list(sim_auc_apoptotic_pairs.values())
    drug_dataframe["log2_auc_apoptotic_ratio"] = list(log2_auc_apoptotic_ratio.values())
    drug_dataframe["drug_1"] = drug_1
    drug_dataframe["drug_2"] = drug_2
    print(list(sim_auc_apoptotic_pairs.values()))
    print(list(log2_auc_apoptotic_ratio.values()))
    drug_dataframe["conc_1"] = conc_1
    drug_dataframe["conc_2"] = conc_2

    # list of dictionaries where each dictionary contains the same drug-pairs

    print(drug_dataframe)
    
    drug_names = ["Ipatasertib", "Afatinib", "Ulixertinib", "Luminespib", "Selumetinib", "Pictilisib"]
    # print(test_drug_dataframe.head())
    # print(test_drug_dataframe.count())
    # create a figure that contains all the subplots 
    fig, axes = plt.subplots(nrows = 6, ncols = 6, sharex='col', sharey='row', figsize=(12,9))
    cbar_ax = fig.add_axes([0.89, 0.45, 0.05, 0.5])
  
    for row in range(6):
        for col in range(6):
            if row < col:
                axes[row, col].axis('off')
            # get the current drugs
            drug_a = drug_names[row]
            drug_b = drug_names[col]
            # get the data for the drugs used 
            # print(drug_a)
            # print(drug_b)
            drug_pair_data = drug_dataframe.loc[(drug_dataframe["drug_1"] == drug_a) & (drug_dataframe["drug_2"] == drug_b)]
            # print(drug_pair_data.head())
            if not drug_pair_data.empty:
                df_wide = drug_pair_data.pivot_table( index= 'conc_2', columns='conc_1', values='log2_auc_live_ratio', aggfunc='first')
                print(df_wide)
                # cmap = sns.diverging_palette(220,20, as_cmap=True)
                ax = sns.heatmap(ax=axes[col, row], data=df_wide, cbar_ax=cbar_ax, cbar_kws={'label': 'log2 (drug_AUC / wildtype_AUC)'}, cmap="RdBu", vmin=-0.5, vmax=0.5)
                ax.invert_yaxis()
                ax.tick_params(axis='x', labelrotation=45)
                axes[col, row].set(xlabel= "", ylabel= "")

    for ax, col in zip(axes[5,:], drug_names):
        ax.set_xlabel(col)

    for ax, row in zip(axes[:,0], drug_names):
        ax.set_ylabel(row, rotation=90, size='large')

    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    fig.suptitle('Growth behaviour of LNCaP upon drug administration with respect to wildtype LNCaP', y=0.98)
    # plt.show()
    plt.savefig('heatmap_live_multiple_wildtype' + '.png')

    # print apoptosis heatmap
    fig, axes = plt.subplots(nrows = 6, ncols = 6, sharex='col', sharey='row', figsize=(12,9))
    cbar_ax = fig.add_axes([0.89, 0.45, 0.05, 0.5])
  
    for row in range(6):
        for col in range(6):
            if row < col:
                axes[row, col].axis('off')
            # get the current drugs
            drug_a = drug_names[row]
            drug_b = drug_names[col]
            # get the data for the drugs used 
            # print(drug_a)
            # print(drug_b)
            drug_pair_data = drug_dataframe.loc[(drug_dataframe["drug_1"] == drug_a) & (drug_dataframe["drug_2"] == drug_b)]
            # print(drug_pair_data.head())
            if not drug_pair_data.empty:
                df_wide = drug_pair_data.pivot_table( index= 'conc_2', columns='conc_1', values='log2_auc_apoptotic_ratio', aggfunc='first')
                print(df_wide)
                # cmap = sns.diverging_palette(220,20, as_cmap=True)
                ax = sns.heatmap(ax=axes[col, row], data=df_wide, cbar_ax=cbar_ax, cbar_kws={'label': 'log2 (drug_AUC / wildtype_AUC)'}, cmap="RdBu", vmin=-0.5, vmax=0.5)
                ax.invert_yaxis()
                ax.tick_params(axis='x', labelrotation=45)
                axes[col, row].set(xlabel= "", ylabel= "")

    for ax, col in zip(axes[5,:], drug_names):
        ax.set_xlabel(col)

    for ax, row in zip(axes[:,0], drug_names):
        ax.set_ylabel(row, rotation=90, size='large')

    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    fig.suptitle('Apoptosis variation of LNCaP upon drug administration with respect to wildtype LNCaP', y=0.98)
    # plt.show()
    plt.savefig('heatmap_apoptosis_multiple_wildtype' + '.png')


    # ############################################################################################
    # # Synergy calculations 
    # ############################################################################################

    # # calculate the bliss independence reference model 
    # double_drugs = drug_dataframe[drug_dataframe["drug_1"] != drug_dataframe["drug_2"]]
    # single_drugs = drug_dataframe[drug_dataframe["drug_1"] == drug_dataframe["drug_2"]]
    # double_drugs["CI"] = 1
    # double_drugs["bliss_independence"] = 0
    # for index, row in double_drugs.iterrows():
    #     print(row['drug_1'])
    #     drug_a = single_drugs[(single_drugs["drug_1"] == row["drug_1"]) & (single_drugs["conc_1"] == row["conc_1"])]
    #     drug_b = single_drugs[(single_drugs["drug_1"] == row["drug_2"]) & (single_drugs["conc_1"] == row["conc_2"])]
    #     print(drug_b)
    #     E_a = drug_a["log2_auc_live_ratio"].iloc[0]
    #     E_b = drug_b["log2_auc_live_ratio"].iloc[0]
    #     # Bliss independence can only be calculated with positive values, but we are interested in the inhibition (negative values)
    #     # so set the positive values to zero and take the absolute values of the negative inhibition values
    #     if E_a > 0:
    #         E_a = 0
    #     else:
    #         abs(E_a)
    #     if E_b > 0:
    #         E_b = 0
    #     else:
    #         abs(E_b)
    #     print(E_a)
    #     print(E_b)
    #     bliss_independence = E_a + E_b - E_a * E_b
    #     print(bliss_independence)
    #     double_drugs.at[index, "bliss_independence"] = bliss_independence
    #      # calculate the combination index 
    #     CI = row["bliss_independence"] / row["log2_auc_live_ratio"]
    #     print(CI)
    #     double_drugs.at[index, "CI"] = CI
    # print(double_drugs)

    # # plot the same heatmap as above just with the combination index values 

    # fig, axes = plt.subplots(nrows = 6, ncols = 6, sharex='col', sharey='row', figsize=(12,9))
    # cbar_ax = fig.add_axes([0.89, 0.45, 0.05, 0.5])
  
    # for row in range(6):
    #     for col in range(6):
    #         if row < col:
    #             axes[row, col].axis('off')
    #         # get the current drugs
    #         drug_a = drug_names[row]
    #         drug_b = drug_names[col]
    #         # get the data for the drugs used 
    #         # print(drug_a)
    #         # print(drug_b)
    #         drug_pair_data = double_drugs.loc[(double_drugs["drug_1"] == drug_a) & (double_drugs["drug_2"] == drug_b)]
    #         # print(drug_pair_data.head())
    #         if not drug_pair_data.empty:
    #             df_wide = drug_pair_data.pivot_table( index= 'conc_2', columns='conc_1', values='CI', aggfunc='first')
    #             print(df_wide)
    #             # cmap = sns.diverging_palette(220,20, as_cmap=True)
    #             ax = sns.heatmap(ax=axes[col, row], data=df_wide, cbar_ax=cbar_ax, cbar_kws={'label': 'Combination index (CI)'}, cmap="RdBu", vmin=-0.5, vmax=0.5)
    #             ax.invert_yaxis()
    #             ax.tick_params(axis='x', labelrotation=45)
    #             axes[col, row].set(xlabel= "", ylabel= "")

    # for ax, col in zip(axes[5,:], drug_names):
    #     ax.set_xlabel(col)

    # for ax, row in zip(axes[:,0], drug_names):
    #     ax.set_ylabel(row, rotation=90, size='large')

    # fig.tight_layout()
    # fig.subplots_adjust(top=0.95)
    # fig.suptitle('Bliss independence - combination index of each drug-pair in LNCaP', y=0.98)
    # # plt.show()
    # plt.savefig('Bliss_independence' + '.png')


main() 

# %%
