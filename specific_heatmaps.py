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
# from matplotlib.colors import DivergingNorm


modules_path = os.path.dirname(os.path.realpath(__file__))
modules_path = os.path.join(modules_path, 'modules')
sys.path.append(modules_path)

import multicellds 


sns.set(style="ticks", palette="Paired")
sns.set_context('paper')


def create_parser():
    parser = argparse.ArgumentParser(description="Plot total cell grouped as Alive/Necrotic/Apoptotic vs Time")
    
    parser.add_argument("drug_data_folder", action="store", help="folder were all the folders with drug data for one cell line is stored (e.g. output) ")

    parser.add_argument("wildtype_folder", action="store", help="folder that contains all replicates of the wildtype simulation (e.g. output_LNCaP")
    
    parser.add_argument("--drug1")
    parser.add_argument("--drug2")
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

    sims_paths = []
    for directory in os.listdir(args.drug_data_folder):
        if args.drug2 is None and  "single" in directory:
            # we just search for the heatmap of a single drug
            for sub_dir in args.drug_data_folder + directory:
                if args.drug1 in sub_dir:
                    sims_paths.append(args.drug_data_folder + directory + "/" + sub_dir)
        else:
            # we search for the heatmap of a pair of drugs 
            if "double" == directory:
                for sub_dir in os.listdir(args.drug_data_folder + directory):
                    if args.drug1 in sub_dir and args.drug2 in sub_dir:
                        sims_paths.append(args.drug_data_folder + directory + "/" + sub_dir)
            if "single" == directory:
                for sub_dir in os.listdir(args.drug_data_folder + directory):
                    if args.drug1 in sub_dir or args.drug2 in sub_dir:
                        sims_paths.append(args.drug_data_folder + directory + "/" + sub_dir)
  
    for directory in os.listdir(args.wildtype_folder):
        if os.path.isdir(args.wildtype_folder + directory):
            sims_paths.append(args.wildtype_folder + directory)

    for el in sims_paths:
        print(el)


    sim_auc_live_pairs = {}
    sim_auc_apoptotic_pairs = {}
    drug_1 = []
    drug_2 = []
    conc_1 = []
    conc_2 = []

    for data_folder in sims_paths:
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
            first_drug = args.drug1
            second_drug = args.drug2
            first_conc = second_conc = "None"
        elif len(used_drugs) == 1:
            if used_drugs.pop(0) == args.drug1:
                first_drug = args.drug1
                # replace drug levels with IC values 
                first_conc = drug_level_to_IC(simulation_name[17])
                second_drug = args.drug2
                second_conc = "None"
            else:
                first_drug = args.drug1
                first_conc = "None"
                second_drug = args.drug2
                second_conc = drug_level_to_IC(simulation_name[17])
                
                
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
    
    drug_names = ["Ipatasertib", "Afatinib", "Ulixertinib", "Luminespib", "Selumetinib", "Pictilisib"]
    # print(test_drug_dataframe.head())
    # print(test_drug_dataframe.count())
    # create a figure that contains all the subplots 
    fig, axes = plt.subplots(sharex='col', sharey='row', figsize=(12,9))
    # cbar_ax = fig.add_axes([0.89, 0.45, 0.05, 0.5])
  
    # print(drug_pair_data.head())
    # drug_df_without_wildtype = drug_dataframe[(drug_dataframe.drug_1 != "wildtype")]
    drug_df_filtered = drug_dataframe[['drug_1', 'drug_2', 'conc_1', 'conc_2', 'log2_auc_live_ratio']]
    # drug_df_doubles = drug_df_filtered[drug_df_filtered["drug_1"] != drug_df_filtered["drug_2"]]
    drug_df_averages = drug_df_filtered.groupby(['drug_1', 'drug_2', 'conc_1', 'conc_2']).mean()
    drug_df_averages = drug_df_averages.reset_index()
    print(drug_df_averages)
    drug_df_averages['conc_1'] = pd.Categorical(drug_df_averages['conc_1'], ordered=True, categories=['None', 'IC10', 'IC30', 'IC50', 'IC70', 'IC90'])
    drug_df_averages['conc_2'] = pd.Categorical(drug_df_averages['conc_2'], ordered=True, categories=['None', 'IC10', 'IC30', 'IC50', 'IC70', 'IC90'])
    output_name = ""
    if args.drug2 is None:
        output_name = "dataframe_" + args.drug1 + ".csv"
    else:
        output_name = "dataframe_" + args.drug1 + "_" + args.drug2 + ".csv"
    drug_df_averages.to_csv(output_name)
    drug_dataframe.to_csv("full_data" + output_name)
    df_wide = drug_df_averages.pivot_table( index= 'conc_2', columns='conc_1', values='log2_auc_live_ratio')
    print(df_wide)
    # cmap = sns.diverging_palette(220,20, as_cmap=True)
    ax = sns.heatmap(data=df_wide, cbar_kws={'label': 'log2 (drug_AUC / wildtype_AUC)'}, cmap="RdBu", vmin=-0.5, vmax=0.5)
    plt.xlabel(args.drug1)
    plt.ylabel(args.drug2)
    ax.invert_yaxis()
    
    # ax.tick_params(axis='x', labelrotation=45)
 
    # for ax, col in zip(axes[5,:], drug_names):
    #     ax.set_xlabel(col)

    # for ax, row in zip(axes[:,0], drug_names):
    #     ax.set_ylabel(row, rotation=90, size='large')

    fig.tight_layout()
    # fig.subplots_adjust(top=0.95)
    # fig.suptitle('Growth behaviour of LNCaP upon drug administration with respect to wildtype LNCaP', y=0.98)
    # plt.show()
    if args.drug2 is None:
        plt.savefig('heatmap_' + args.drug1 + '.png')
    else:
        plt.savefig('heatmap_' + args.drug1 + "_" + args.drug2 + '.png')

    # # print apoptosis heatmap
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
    #         drug_pair_data = drug_dataframe.loc[(drug_dataframe["drug_1"] == drug_a) & (drug_dataframe["drug_2"] == drug_b)]
    #         # print(drug_pair_data.head())
    #         if not drug_pair_data.empty:
    #             df_wide = drug_pair_data.pivot_table( index= 'conc_2', columns='conc_1', values='log2_auc_apoptotic_ratio', aggfunc='first')
    #             print(df_wide)
    #             # cmap = sns.diverging_palette(220,20, as_cmap=True)
    #             ax = sns.heatmap(ax=axes[col, row], data=df_wide, cbar_ax=cbar_ax, cbar_kws={'label': 'log2 (drug_AUC / wildtype_AUC)'}, cmap="RdBu", vmin=-0.5, vmax=0.5)
    #             ax.invert_yaxis()
    #             ax.tick_params(axis='x', labelrotation=45)
    #             axes[col, row].set(xlabel= "", ylabel= "")

    # for ax, col in zip(axes[5,:], drug_names):
    #     ax.set_xlabel(col)

    # for ax, row in zip(axes[:,0], drug_names):
    #     ax.set_ylabel(row, rotation=90, size='large')

    # fig.tight_layout()
    # fig.subplots_adjust(top=0.95)
    # fig.suptitle('Apoptosis variation of LNCaP upon drug administration with respect to wildtype LNCaP', y=0.98)
    # # plt.show()
    # plt.savefig('heatmap_apoptosis_multiple_wildtype' + '.png')


    ############################################################################################
    # Synergy calculations 
    ############################################################################################

    # only calculate bliss in the case that we have two drugs 
    if(args.drug2 != None):
        # calculate the bliss independence reference model 
        double_drugs = drug_df_averages.loc[(drug_df_averages["conc_1"] != "None") & (drug_df_averages["conc_2"] != "None")]
        single_drugs = drug_df_averages.loc[(drug_df_averages["conc_1"] == "None") | (drug_df_averages["conc_2"] == "None")]
        double_drugs["CI"] = 1.0
        double_drugs["bliss_independence"] = 0.0
        print("single drugs:")
        print(single_drugs)
        for index, row in double_drugs.iterrows():
            print("row: ")
            print(row)
            print(row["drug_1"])
            drug_a = single_drugs[(single_drugs["drug_1"] == row["drug_1"]) & (single_drugs["conc_1"] == row["conc_1"])]
            print("Drug a:")
            print(drug_a)
            drug_b = single_drugs[(single_drugs["drug_2"] == row["drug_2"]) & (single_drugs["conc_2"] == row["conc_2"])]
            print("drug B:")
            print(drug_b)
            E_a = drug_a["log2_auc_live_ratio"].iloc[0]
            E_b = drug_b["log2_auc_live_ratio"].iloc[0]
            print(E_a)
            print(E_b)
            
            # Bliss independence can only be calculated with positive values, but we are interested in the inhibition (negative values)
            # so set the positive values to zero and take the absolute values of the negative inhibition values
            if E_a > 0:
                E_a = 0
            else:
                E_a = abs(E_a)
            if E_b > 0:
                E_b = 0
            else:
                E_b = abs(E_b)

            bliss_independence = E_a + E_b - E_a * E_b
            print("Bliss independence:")
            print(bliss_independence)
            double_drugs.at[index, "bliss_independence"] = bliss_independence

            # calculate the combination index 
            CI = double_drugs.at[index, "bliss_independence"] / row["log2_auc_live_ratio"]
            double_drugs.at[index, "CI"] = CI
        print(double_drugs)

        double_df_filtered = double_drugs[['drug_1', 'drug_2', 'conc_1', 'conc_2', 'CI']]
        # double_df_averages = drug_df_filtered.groupby(['drug_1', 'drug_2', 'conc_1', 'conc_2']).mean()
        df_wide = double_df_filtered.pivot_table( index= 'conc_2', columns='conc_1', values='CI', aggfunc='first')
        print(df_wide)
        if args.drug2 is None:
            df_wide.to_csv('bliss_dataframe_' + args.drug1 + '.csv')
        else:
            df_wide.to_csv('bliss_dataframe_' + args.drug1 + "_" + args.drug2 + '.csv')
        
        fig, axes = plt.subplots(sharex='col', sharey='row', figsize=(12,9))
        # cmap = sns.diverging_palette(220,20, as_cmap=True)
        # norm = DivergingNorm(vmin=0, vcenter=1, vmax=10)
        ax = sns.heatmap(data=df_wide, cbar_kws={'label': 'Combination Index (CI)', 'ticks': [0,1,2,3,4,5,6,7,8,9,10]}, cmap="RdBu", vmin=0, center=1, vmax=10)
        ax.invert_yaxis()
        plt.xlabel(args.drug1)
        plt.ylabel(args.drug2)
        
        # ax.tick_params(axis='x', labelrotation=45)
    
        # for ax, col in zip(axes[5,:], drug_names):
        #     ax.set_xlabel(col)

        # for ax, row in zip(axes[:,0], drug_names):
        #     ax.set_ylabel(row, rotation=90, size='large')

        fig.tight_layout()
        # fig.subplots_adjust(top=0.95)
        # fig.suptitle('Growth behaviour of LNCaP upon drug administration with respect to wildtype LNCaP', y=0.98)
        # plt.show()
        if args.drug2 is None:
            plt.savefig('bliss_' + args.drug1 + '.png')
        else:
            plt.savefig('bliss_' + args.drug1 + "_" + args.drug2 + '.png')


main() 

# %%
