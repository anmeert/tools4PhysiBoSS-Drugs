#!/usr/bin/env python3
# coding: utf-8
# author: Annika Meert

import os, re, sys
import argparse

import numpy as np
import pandas as pd
import scipy.integrate as integrate

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

class MidpointNormalize(mcolors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        v_ext = np.max( [ np.abs(self.vmin), np.abs(self.vmax) ] )
        # for heatmaps
        x, y = [-v_ext, self.midpoint, v_ext], [0, 0.5, 0.9]
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

############################################
# create heatmaps 
#############################################

font = {'size': 25}
drug_df_medians['conc_1'] = pd.Categorical(drug_df_medians['conc_1'], ordered=True, categories=['None', 'IC10', 'IC30', 'IC50', 'IC70', 'IC90'])
drug_df_medians['conc_2'] = pd.Categorical(drug_df_medians['conc_2'], ordered=True, categories=['None', 'IC10', 'IC30', 'IC50', 'IC70', 'IC90'])
df_wide = drug_df_medians.pivot_table( index= 'conc_2', columns='conc_1', values='log2_ratio')
sns.set(font_scale=2.0)
fig, axes = plt.subplots(sharex='col', sharey='row', figsize=(12,9))
cmap = "RdBu"
norm = MidpointNormalize( midpoint = 0, vmin=-0.32, vmax=0.2)
ax = sns.heatmap(data=df_wide, cbar_kws={'label': 'growth index'}, cmap=cmap, norm=norm)
ax.invert_yaxis()
plt.xlabel(drug1, labelpad=15)
plt.ylabel(drug2, labelpad=15)


fig.tight_layout()
plt.savefig('output/heatmap_pair_' + drug1 + "_" + drug2 + '.png', dpi=300)