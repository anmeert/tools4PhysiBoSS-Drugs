#!/usr/bin/env python3
# coding: utf-8
# 
# Calculate Kruskal Wallis for all drug simulations

import pandas as pd
from scipy.stats import stats
from itertools import product, combinations, permutations
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations, product

# load the csv 
drug_dataframe = pd.read_csv("./data/LNCaP_simulation_data.csv")

# go through all drug combinations
drug_names = ["Ipatasertib", "Afatinib", "Ulixertinib", "Luminespib", "Selumetinib", "Pictilisib"]
drug_combs = combinations(drug_names, 2)
drug_concs = ["IC10", "IC30", "IC50", "IC70", "IC90"]
conc_combs = product(drug_concs, repeat=2)
conc_combs_list = list(conc_combs)
significance_list = []
drug_1 = []
drug_2 = []

for drug1, drug2 in drug_combs:
        drug_pair_data = drug_dataframe.loc[(drug_dataframe["drug_1"] == drug1) & (drug_dataframe["drug_2"] == drug2)]
        # add the no-drug replicates to the replicate dictionary 
        no_drug_name = "output_LNCaP"
        no_drug_data = drug_dataframe.loc[(drug_dataframe['conc_1'] == "None") & (drug_dataframe['conc_2'] == "None")]
        AUC_no_drug = no_drug_data['auc_live']
        conc_replicates = [AUC_no_drug]
        # for each drug combination group the simulations with the same drug concentration and append to the replicate dictionary
        for conc1, conc2 in conc_combs_list:
            replicates = drug_pair_data.loc[(drug_pair_data['conc_1'] == conc1) & (drug_pair_data['conc_2'] == conc2)]
            AUC_replicate_array = replicates['auc_live']
            conc_replicates = np.append(conc_replicates, [AUC_replicate_array], 0)
        H, pval = stats.kruskal(*conc_replicates)
        significance_list.append(pval)
        drug_1.append(drug1)
        drug_2.append(drug2)

significance = pd.DataFrame(zip(significance_list, drug_1, drug_2), columns=['significance','drug_1','drug_2'])
print(significance)

