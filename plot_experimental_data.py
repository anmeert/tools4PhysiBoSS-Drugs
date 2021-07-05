#!/usr/bin/env python3
# coding: utf-8
# author: Annika Meert 
# plot experimental LNCaP data 
import os
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt

def normalize (row, norm_3_3, norm_10, norm_30):
    if row["uM"] == 3.3:
        return row["mean"] + norm_3_3
    if row["uM"] == 10:
        return row["mean"] + norm_10
    else: 
        return row["mean"] + norm_30

# load the csv 
drug_dataframe = pd.read_csv("./data/LNCaP_experimental_data.csv")
print(drug_dataframe.head())

selum_data = drug_dataframe.loc[(drug_dataframe['drug'] == "Selumetinib") & (drug_dataframe['time2'] > 0) & (drug_dataframe['time2'] < 14)]
selum_3_3_norm = 1 - selum_data[selum_data["uM"] == 3.3]["mean"].iloc[0]
selum_10_norm = 1 - selum_data[selum_data["uM"] == 10]["mean"].iloc[0]
selum_30_norm = 1 - selum_data[selum_data["uM"] == 30]["mean"].iloc[0]
selum_data["mean_norm"] = selum_data.apply (lambda row: normalize(row, selum_3_3_norm, selum_10_norm, selum_30_norm), axis=1)
print(selum_data)
picti_data = drug_dataframe.loc[(drug_dataframe['drug'] == "Pictilisib") & (drug_dataframe['time2'] > 0) & (drug_dataframe['time2'] < 14)]
picti_3_3_norm = 1 - picti_data[picti_data["uM"] == 3.3]["mean"].iloc[0]
picti_10_norm = 1 - picti_data[picti_data["uM"] == 10]["mean"].iloc[0]
picti_30_norm = 1 - picti_data[picti_data["uM"] == 30]["mean"].iloc[0]
picti_data["mean_norm"] = picti_data.apply (lambda row: normalize(row, picti_3_3_norm, picti_10_norm, picti_30_norm), axis=1)
lumen_data = drug_dataframe.loc[(drug_dataframe['drug'] == "NMS-E973") & (drug_dataframe['time2'] > 0) & (drug_dataframe['time2'] < 14)]
lumen_3_3_norm = 1 - lumen_data[lumen_data["uM"] == 3.3]["mean"].iloc[0]
lumen_10_norm = 1 - lumen_data[lumen_data["uM"] == 10]["mean"].iloc[0]
lumen_30_norm = 1 - lumen_data[lumen_data["uM"] == 30]["mean"].iloc[0]
lumen_data["mean_norm"] = lumen_data.apply (lambda row: normalize(row, lumen_3_3_norm, lumen_10_norm, lumen_30_norm), axis=1)
control_data = drug_dataframe.loc[(drug_dataframe['drug'] == "Control") & (drug_dataframe['time2'] > 0) & (drug_dataframe['time2'] < 14)]
control_data["mean_norm"] = control_data["mean"]

sns.color_palette("hls", 8)
fig1, ax1 = plt.subplots()
ax1 = sns.lineplot(data=selum_data.append(control_data), x='time2', y='mean_norm', hue="uM", palette="tab10")
plt.xlabel("Time (hours)")
plt.ylabel("Centered cell Index")
leg=ax1.get_legend()
title = 'Drug dose (uM)'
leg.set_title(title)
if not os.path.exists('output'):
    os.makedirs('output')
fig1.savefig("output/Selum_experimental_data.png", dpi= 300)

fig2, ax2 = plt.subplots()
ax2 = sns.lineplot(data=picti_data.append(control_data), x='time2', y='mean_norm', hue="uM", palette="tab10")
plt.xlabel("Time (hours)")
plt.ylabel("Centered cell Index")
leg=ax2.get_legend()
title = 'Drug dose (uM)'
leg.set_title(title)
fig2.savefig("output/Picti_experimental_data.png", dpi= 300)

fig3, ax3 = plt.subplots()
ax3 = sns.lineplot(data=lumen_data.append(control_data), x='time2', y='mean_norm', hue="uM", palette="tab10")
plt.xlabel("Time (hours)")
plt.ylabel("Centered cell Index")
leg=ax3.get_legend()
title = 'Drug dose (uM)'
leg.set_title(title)
fig3.savefig("output/Lumen_experimental_data.png", dpi= 300)