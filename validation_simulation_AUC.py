#!/usr/bin/env python3
# coding: utf-8
# author: Annika Meert

import os, re, sys
import numpy as np
import pandas as pd

drug_dataframe = pd.read_csv("./data/validation_Selumetinib_30uM.csv")
drug_df_filtered = drug_dataframe[['drug_1',  'conc_1', 'conc_2', 'log2_auc_live_ratio']]
drug_df_averages = drug_df_filtered.groupby(['drug_1', 'conc_1', 'conc_2']).mean()
drug_df_averages = drug_df_averages.reset_index()
print(drug_df_averages)
