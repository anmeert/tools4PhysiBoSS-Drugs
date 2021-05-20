#!/usr/bin/env python3
# coding: utf-8
# author: Annika Meert

import os, re, sys
import glob, json
import argparse

import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy.io import loadmat

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from itertools import combinations, product
from matplotlib.colors import DivergingNorm

class MidpointNormalize(mcolors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        v_ext = np.max( [ np.abs(self.vmin), np.abs(self.vmax) ] )
        # for bliss independence
        x, y = [-v_ext, self.midpoint, v_ext], [-5, 0.5, 0.7]
        return np.ma.masked_array(np.interp(value, x, y))

drug_dataframe = pd.read_csv("./data/LNCaP_simulation_data.csv")
drug1 = "Ipatasertib"
drug2 = "Pictilisib"

# calculate the medians for all replicates
drug_df_medians = drug_dataframe.groupby(['drug_1', 'drug_2', 'conc_1', 'conc_2'])['auc_live'].median().reset_index()

# calculate the log2_ratio
median_no_drug_auc = drug_df_medians.loc[(drug_df_medians['conc_1'] == "None") & (drug_df_medians['conc_2'] == "None")].iloc[0]['auc_live']
drug_df_medians['log2_ratio'] = np.log2(drug_df_medians['auc_live'] / median_no_drug_auc)
drug_df_medians['scaled_auc'] = 1 - ((drug_df_medians['auc_live'] - drug_df_medians['auc_live'].min()) / (drug_df_medians['auc_live'].max() - drug_df_medians['auc_live'].min()))

# filter out double, single and wildtype for the wanted drugs
drug_df_medians = drug_df_medians.loc[((drug_df_medians['drug_1'] == drug1) & (drug_df_medians['drug_2'] == drug2)) ^ (drug_df_medians['drug_1'] == 'wildtype') \
    ^ (drug_df_medians['drug_1'] == drug1) & (drug_df_medians['drug_2'] == drug1) ^(drug_df_medians['drug_1'] == drug2) & (drug_df_medians['drug_2'] == drug2)]

# rename single concentrations to None 
drug_df_medians.conc_2[((drug_df_medians['drug_1'] == drug1) & (drug_df_medians['drug_2'] == drug1) ^ (drug_df_medians['drug_1'] == 'wildtype'))] = "None"
drug_df_medians.conc_1[((drug_df_medians['drug_1'] == drug2) & (drug_df_medians['drug_2'] == drug2) ^ (drug_df_medians['drug_1'] == 'wildtype'))] = "None"

#reorder concentration values 
drug_df_medians['conc_1'] = pd.Categorical(drug_df_medians['conc_1'], ordered=True, categories=['None', 'IC10', 'IC30', 'IC50', 'IC70', 'IC90'])
drug_df_medians['conc_2'] = pd.Categorical(drug_df_medians['conc_2'], ordered=True, categories=['None', 'IC10', 'IC30', 'IC50', 'IC70', 'IC90'])

############################################################################################
# Synergy calculations 
############################################################################################

# calculate the bliss independence reference model 
double_drugs = drug_df_medians.loc[(drug_df_medians["conc_1"] != "None") & (drug_df_medians["conc_2"] != "None")]
single_drugs = drug_df_medians.loc[(drug_df_medians["conc_1"] == "None") ^ (drug_df_medians["conc_2"] == "None")]
double_drugs["CI"] = 1.0
double_drugs["bliss_independence"] = 0.0
for index, row in double_drugs.iterrows():
    drug_a = single_drugs[(single_drugs["drug_1"] == row["drug_1"]) & (single_drugs["conc_1"] == row["conc_1"])]
    drug_b = single_drugs[(single_drugs["drug_2"] == row["drug_2"]) & (single_drugs["conc_2"] == row["conc_2"])]
    E_a = drug_a["scaled_auc"].iloc[0]
    E_b = drug_b["scaled_auc"].iloc[0]
    
    bliss_independence = E_a + E_b - E_a * E_b
    double_drugs.at[index, "bliss_independence"] = bliss_independence

    # calculate the combination index 
    CI = double_drugs.at[index, "bliss_independence"] / abs(row["scaled_auc"])
    double_drugs.at[index, "CI"] = CI

double_df_filtered = double_drugs[['drug_1', 'drug_2', 'conc_1', 'conc_2', 'CI']]
# double_df_averages = drug_df_filtered.groupby(['drug_1', 'drug_2', 'conc_1', 'conc_2']).mean()
df_wide = double_df_filtered.pivot_table( index= 'conc_2', columns='conc_1', values='CI', aggfunc='first')

fig, axes = plt.subplots(sharex='col', sharey='row', figsize=(12,9))
# cmap = sns.diverging_palette(220,20, as_cmap=True)
cmap = plt.get_cmap("BrBG_r")
x = np.arange( 0, 1, 1e-1 )
xlen = x.shape[ 0 ]
z = np.random.random( xlen**2 )*12 - 2
norm = MidpointNormalize( midpoint = 1 , vmin= 0.5 , vmax=5.0)
ax = sns.heatmap(data=df_wide, cbar_kws={'label': 'Combination Index (CI)'}, cmap=cmap, norm=norm)
ax.invert_yaxis()
plt.xlabel(drug1, labelpad=15)
plt.ylabel(drug2, labelpad=15)

fig.tight_layout()
fig.subplots_adjust(top=0.95)
fig.suptitle('Growth behaviour of LNCaP upon drug administration with respect to wildtype LNCaP', y=0.98)
if not os.path.exists('output'):
    os.makedirs('output')
plt.savefig('output/bliss_pair_' + drug1 + "_" + drug2 + '.png')



