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
from matplotlib.colors import DivergingNorm

drug_df_averages = pd.read_csv("/home/bscuser/scripts/PhysiBoSS_data_analysis/dataframe_Luminespib_Pictilisib.csv")
print(drug_df_averages)
drug1 = "Luminespib"
drug2 = "Pictilisib"

############################################################################################
# Synergy calculations 
############################################################################################



# calculate the bliss independence reference model 
double_drugs = drug_df_averages.loc[(drug_df_averages["conc_1"] != "None") & (drug_df_averages["conc_2"] != "None")]
single_drugs = drug_df_averages.loc[(drug_df_averages["conc_1"] == "None") ^ (drug_df_averages["conc_2"] == "None")]
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
    CI = double_drugs.at[index, "bliss_independence"] / abs(row["log2_auc_live_ratio"])
    double_drugs.at[index, "CI"] = CI
print(double_drugs)

double_df_filtered = double_drugs[['drug_1', 'drug_2', 'conc_1', 'conc_2', 'CI']]
# double_df_averages = drug_df_filtered.groupby(['drug_1', 'drug_2', 'conc_1', 'conc_2']).mean()
df_wide = double_df_filtered.pivot_table( index= 'conc_2', columns='conc_1', values='CI', aggfunc='first')
print(df_wide)
if drug2 is None:
    df_wide.to_csv('bliss_dataframe_' + drug1 + '.csv')
else:
    df_wide.to_csv('bliss_dataframe_' + drug1 + "_" + drug2 + '.csv')

fig, axes = plt.subplots(sharex='col', sharey='row', figsize=(12,9))
# cmap = sns.diverging_palette(220,20, as_cmap=True)
norm = DivergingNorm(vmin=0, vcenter=1, vmax=5)
ax = sns.heatmap(data=df_wide, cbar_kws={'label': 'Combination Index (CI)', 'ticks': [0,1,2,3,4,5,6,7,8,9,10]}, cmap="BrBG_r", norm=norm)
ax.invert_yaxis()
plt.xlabel(drug1)
plt.ylabel(drug2)
    
# ax.tick_params(axis='x', labelrotation=45)

# for ax, col in zip(axes[5,:], drug_names):
#     ax.set_xlabel(col)

# for ax, row in zip(axes[:,0], drug_names):
#     ax.set_ylabel(row, rotation=90, size='large')

fig.tight_layout()
# fig.subplots_adjust(top=0.95)
# fig.suptitle('Growth behaviour of LNCaP upon drug administration with respect to wildtype LNCaP', y=0.98)
# plt.show()
if drug2 is None:
    plt.savefig('bliss_' + drug1 + '.png')
else:
    plt.savefig('bliss_' + drug1 + "_" + drug2 + '.png')



