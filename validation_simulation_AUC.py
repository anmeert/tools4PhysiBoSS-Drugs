#!/usr/bin/env python3
# coding: utf-8
# author: Annika Meert

import os, re, sys
import numpy as np
import pandas as pd

drug_dataframe = pd.read_csv("./data/validation_Pictilisib_10uM.csv")
drug_df_filtered = drug_dataframe[['drug_1',  'conc_1', 'conc_2', 'auc_live']]
drug_df_medians = drug_df_filtered.groupby(['drug_1', 'conc_1', 'conc_2']).median()
drug_df_medians = drug_df_medians.reset_index()
print(drug_df_medians)
drug_auc = drug_df_medians.loc[drug_df_medians["conc_1"] != "None"].iloc[0]["auc_live"]
no_drug_auc = drug_df_medians.loc[drug_df_medians["conc_1"] == "None"].iloc[0]["auc_live"]
log2_auc_ratio = np.log2(drug_auc/no_drug_auc)
print(log2_auc_ratio)
