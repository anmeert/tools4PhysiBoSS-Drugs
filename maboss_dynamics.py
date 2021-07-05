 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

lncap_data = pd.read_csv("data/MaBoSS_data_LNCaP.tsv", sep='\t')
nodes = ["AKT", "EGFR", "ERK", "HSPs", "MEK1_2", "PI3K"]

single = lncap_data[lncap_data["TYPE"] == "SINGLE"]
single_drug = single[(single["gene1"].isin(nodes)) | (single["gene2"].isin(nodes))]
single_drug_ko = single_drug[(single_drug["mutant1"] == "KO") | (single_drug["mutant2"] == "KO")]
single_melted = single_drug_ko.melt("gene1", value_vars=["Apoptosis", "Proliferation"], var_name="Apop&Prolif", value_name="increase")
print(single_melted)
single_melted["increase"] = single_melted["increase"].astype("float")
ax = sns.barplot(x="gene1", y = "increase",hue="Apop&Prolif", data=single_melted)
ax.set(xlabel="Node", ylabel="Phenotype variation")
ax.legend(title = "Phenotype")
plt.savefig("output/Phenotype_vatiation_MaBoSS.png")