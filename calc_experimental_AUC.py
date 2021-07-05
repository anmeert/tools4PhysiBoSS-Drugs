# author Annika Meert
#%%
import pandas as pd 
import numpy as np

#%%
data = pd.read_csv("./data/LNCaP_experimental_data.csv")
data.head()
# %%
pd.set_option("Display.max_rows", None)
# print(list(data.columns))
# retrieve the wildtype dataframe
wildtype_df = data.loc[(data["drug"] == "Control") & (data["time2"] < 15.0) & (data["time2"] > 0)]
#print(wildtype_df)

# retrieve the dataframe for Pictilisib
picti_df = data.loc[(data["drug"] == "Pictilisib") & (data["uM"] == 10.0) & (data["time2"] < 15.0) & (data["time2"] > 0)]
picti_norm = 1 - picti_df["mean"].iloc[0]
picti_df["mean_norm"] = picti_df["mean"] + picti_norm
# retrieve dataframe for HSPs target 
HSPs_df = data.loc[(data["drug"] == "NMS-E973") & (data["uM"] == 3.3) & (data["time2"] < 15.0) & (data["time2"] > 0)]
HSPs_norm = 1 - HSPs_df["mean"].iloc[0]
HSPs_df["mean_norm"] = HSPs_df["mean"] + HSPs_norm
# retrieve data for Selumetinib
selum_df = data.loc[(data["drug"] == "Selumetinib") & (data["uM"] == 10.0) & (data["time2"] < 15.0) & (data["time2"] > 0)]
selum_norm = 1 - selum_df["mean"].iloc[0]
selum_df["mean_norm"] = selum_df["mean"] + selum_norm

# %%
# calculate the log2(auc_ratio) for each drug 
auc_wildtype = np.trapz(wildtype_df["mean"], wildtype_df["time2"])
print(auc_wildtype)
auc_picti = np.trapz(picti_df["mean_norm"], picti_df["time2"])
print("Pictilisib AUC:" + str(auc_picti))
auc_HSPs = np.trapz(HSPs_df["mean_norm"], HSPs_df["time2"])
print("Luminespib AUC:" + str(auc_HSPs))
auc_selum = np.trapz(selum_df["mean_norm"], selum_df["time2"])
print("Selumetinib AUC:" + str(auc_selum))

log2_auc_ratio_picti = np.log2(auc_picti/auc_wildtype)
print("Pictilisib log2 auc ratio:" + str(log2_auc_ratio_picti))
log2_auc_ratio_HSPs = np.log2(auc_HSPs/auc_wildtype)
print("Luminespib log2 auc ratio:" + str(log2_auc_ratio_HSPs))
log2_auc_ratio_selum = np.log2(auc_selum/auc_wildtype)
print("Selumetinib log2 auc ratio:" + str(log2_auc_ratio_selum))

# %%
