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
        # for heatmaps
        x, y = [-v_ext, self.midpoint, v_ext], [0, 0.5, 0.9]
        return np.ma.masked_array(np.interp(value, x, y))
 
drug_dataframe = pd.read_csv("./data/LNCaP_simulation_data.csv")

drug_names = ["Ipatasertib", "Afatinib", "Ulixertinib", "Luminespib", "Selumetinib", "Pictilisib"]

# create a figure that contains all the subplots 
fig, axes = plt.subplots(nrows = 6, ncols = 6, sharex='col', sharey='row', figsize=(12,9))
cmap = "RdBu"
cbar_ax = fig.add_axes([0.89, 0.45, 0.05, 0.5])

dataframe_medians = drug_dataframe.groupby(['drug_1', 'drug_2', 'conc_1', 'conc_2'])['auc_live'].median().reset_index()
median_no_drug_auc = dataframe_medians.loc[(dataframe_medians['conc_1'] == "None") & (dataframe_medians['conc_2'] == "None")].iloc[0]['auc_live']
dataframe_medians['log2_ratio'] = np.log2(dataframe_medians['auc_live'] / median_no_drug_auc)
dataframe_medians['scaled_auc'] = 1 - ((dataframe_medians['auc_live'] - dataframe_medians['auc_live'].min()) / (dataframe_medians['auc_live'].max() - dataframe_medians['auc_live'].min()))
norm =  MidpointNormalize( midpoint = 0, vmin=dataframe_medians["log2_ratio"].min(), vmax=dataframe_medians["log2_ratio"].max())

for row in range(6):
    for col in range(6):
        if row < col:
            axes[row, col].axis('off')
        # get the current drugs
        drug_a = drug_names[row]
        drug_b = drug_names[col]
        # get the data for the drugs used 
        drug_pair_data = dataframe_medians.loc[(dataframe_medians["drug_1"] == drug_a) & (dataframe_medians["drug_2"] == drug_b)]
        if not drug_pair_data.empty:
            df_wide = drug_pair_data.pivot_table( index= 'conc_2', columns='conc_1', values='log2_ratio', aggfunc='first')
            ax = sns.heatmap(ax=axes[col, row], data=df_wide, cbar_ax=cbar_ax, cbar_kws={'label': 'growth index'}, cmap=cmap, norm=norm)
            ax.invert_yaxis()
            ax.tick_params(axis='x', labelrotation=45)
            axes[col, row].set(xlabel= "", ylabel= "")

for ax, col in zip(axes[5,:], drug_names):
    ax.set_xlabel(col)

for ax, row in zip(axes[:,0], drug_names):
    ax.set_ylabel(row, rotation=90, size='large')

fig.tight_layout()
fig.subplots_adjust(top=0.95)
fig.suptitle('Growth behaviour of LNCaP upon drug treatment with respect to no-drug LNCaP', y=0.98)
plt.savefig('output/full_growth_heatmap' + '.png')


