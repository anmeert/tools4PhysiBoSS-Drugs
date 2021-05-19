#!/usr/bin/env python3
# coding: utf-8
# author: Annika Meert 
# plot experimental LNCaP data 

import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt

# load the csv 
drug_dataframe = pd.read_csv("./LNCaP_experimental_data.csv")
print(drug_dataframe.head())

selum_data = drug_dataframe.loc[(drug_dataframe['drug'] == "Selumetinib") & (drug_dataframe['time2'] > 0) & (drug_dataframe['time2'] < 15)]
print(selum_data)
picti_data = drug_dataframe.loc[(drug_dataframe['drug'] == "Pictilisib") & (drug_dataframe['time2'] > 0) & (drug_dataframe['time2'] < 15)]
lumen_data = drug_dataframe.loc[(drug_dataframe['drug'] == "NMS-E973") & (drug_dataframe['time2'] > 0) & (drug_dataframe['time2'] < 15)]
control_data = drug_dataframe.loc[(drug_dataframe['drug'] == "Control") & (drug_dataframe['time2'] > 0) & (drug_dataframe['time2'] < 15)]

sns.color_palette("hls", 8)
fig1, ax1 = plt.subplots()
ax1 = sns.lineplot(data=selum_data.append(control_data), x='time2', y='mean', hue="uM", palette="tab10")
plt.xlabel("Time (hours)")
plt.ylabel("Cell Index (a.u.)")
leg=ax1.get_legend()
title = 'Drug dose (uM)'
leg.set_title(title)
fig1.savefig("output/Selum_experimental_data.png", dpi= 300)

fig2, ax2 = plt.subplots()
ax2 = sns.lineplot(data=picti_data.append(control_data), x='time2', y='mean', hue="uM", palette="tab10")
plt.xlabel("Time (hours)")
plt.ylabel("Cell Index (a.u.)")
leg=ax2.get_legend()
title = 'Drug dose (uM)'
leg.set_title(title)
fig2.savefig("output/Picti_experimental_data.png", dpi= 300)

fig3, ax3 = plt.subplots()
ax3 = sns.lineplot(data=lumen_data.append(control_data), x='time2', y='mean', hue="uM", palette="tab10")
plt.xlabel("Time (hours)")
plt.ylabel("Cell Index (a.u.)")
leg=ax3.get_legend()
title = 'Drug dose (uM)'
leg.set_title(title)
fig3.savefig("output/Lumen_experimental_data.png", dpi= 300)