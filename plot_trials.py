#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations, product
# trial for plots 

#%%
values = [[1 ,2 ,4], [5, 6, 8], [3,1,0]]
columns = ["Bla", "Bla1", "Bla2"]
rows = ["row1", "row2", "row3"]
# plot the differences for each drug pair 
df = pd.DataFrame(data=values, columns = columns, index=rows)

# fig, ax = plt.subplots(figsize=(11,9))
# sns.heatmap(df, cmap="Blues", linewidth=0.3, cbar_kws={"shrink": 0.8})
# plt.set_xticklabels(rotation=30)
# fig_heat = sns_heatmap.get_figure()
# fig_heat.savefig("test.png")

heatmap = sns.heatmap(values)
heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=30)
plt.xlabel("drug 1 with different concentrations")
plt.ylabel("drug 2 with different concentrations")
plt.show()

# %%

iris = sns.load_dataset("iris")
g = sns.PairGrid(iris)
g.map(sns.scatterplot)

# %%
iris.head()
print(type(iris))


#%%
g = sns.PairGrid(iris, diag_sharey=False)
g.map_upper(sns.heatmap)
g.map_lower(sns.kdeplot)
g.map_diag(sns.kdeplot)



# %%
#### old code 
# # for each simulation build up the dataframe with correct column and row names for plotting it with seaborn 
    # for drug_1, drug_2 in drug_name_combinations:
    #     # all simulations for one drug pair 
    #     sims = { key:value for key,value in auc_differences.items() if all([drug_1 in key, drug_2 in key])}
    #     columns = []
    #     rows = []
    #     values = np.zeros((5,5))
    #     if sims:
    #         for conc_1, conc_2 in drug_conc_combinations:
    #             if drug_1 + "_" + str(conc_1) not in columns:
    #                 columns.append(drug_1 + "_" + str(conc_1))
    #             if  drug_2 + "_" +  str(conc_2) not in rows:
    #                 rows.append(drug_2 + "_" + str(conc_2))
    #             # set the auc difference data for the current drug pair and concentrations if available
    #             current_sim = []
    #             for k,v in sims.items():
    #                 base_name = os.path.basename(os.path.normpath(k))
    #                 first_conc_location = 17
    #                 second_conc_location = 19
    #                 if all([base_name[first_conc_location] == str(conc_1), base_name[second_conc_location] == str(conc_2)]):
    #                     values[conc_2 - 1 ,conc_1 - 1] = auc_differences.get(k)
        
    #         # plot the differences for each drug pair 
    #         df = pd.DataFrame(data=values, columns = columns, index=rows)
            
    #         fig, ax = plt.subplots(figsize=(11,9))
    #         heatmap = sns.heatmap(df, cmap="Blues", linewidth=0.3, cbar_kws={"shrink": 0.8})
    #         heatmap.set_xticklabels(heatmap.get_xticklabels(),rotation=30)
    #         fig_heat = heatmap.get_figure()
    #         print("double_" + str(conc_1) + "_" + str(drug_1) + "_" + str(conc_2) + "_" + str(drug_2) + ".png")
    #         fig_heat.savefig("double_" + str(conc_1) + "_" + str(drug_1) + "_" + str(conc_2) + "_" + str(drug_2) + ".png")
            

        

        # plot all together in one plot 
     
#%%
        
drug_names = ["Ipatasertib", "Afatinib", "Ulixertinib", "Luminespib", "Selumetinib", "Pictilisib"]
combs = combinations(drug_names, 2)
conc = [1,2,3,4,5]
conc_combs = product(conc, repeat=2)
list_concs = list(conc_combs)
#len(list(combs))*len(list(conc_combs))
test_drug_dataframe = pd.DataFrame(columns=['drug_1', 'drug_2', 'conc_1', 'conc_2', 'auc_difference'])

for index, pair in enumerate(combs):
    # print(pair) 
    for i, conc in enumerate(list_concs):
        # print(conc)
        new_list = {'drug_1':pair[0], 'drug_2': pair[1], 'conc_1':conc[0], 'conc_2':conc[1], 'auc_difference':1}
        # print(i)
        test_drug_dataframe.loc[len(test_drug_dataframe.index)] = new_list
        # test_drug_dataframe.append(new_list,ignore_index=True,sort=False)
# print(test_drug_dataframe.head())
# print(test_drug_dataframe.count())
# create a figure that contains all the subplots 
fig, axes = plt.subplots(nrows = 6, ncols = 6, sharex='col', sharey='row', figsize=(12,8))
cbar_ax = fig.add_axes([0.9, 0.45, 0.05, 0.5])

for row in range(6):
    for col in range(6):
        if row < col:
            axes[row, col].axis('off')
        # get the current drugs
        drug_a = drug_names[row]
        drug_b = drug_names[col]
        # get the data for the drugs used 
        # print(drug_a)
        # print(drug_b)
        drug_pair_data = test_drug_dataframe.loc[(test_drug_dataframe["drug_1"] == drug_a) & (test_drug_dataframe["drug_2"] == drug_b)]
        # print(drug_pair_data.head())
        if not drug_pair_data.empty:
            df_wide = drug_pair_data.pivot_table( index= 'conc_1', columns='conc_2', values='auc_difference', aggfunc='first')
            # print(df_wide)
            sns.heatmap(ax=axes[col,row], data=df_wide, cbar_ax=cbar_ax, cmap="vlag")
            axes[col,row].set(xlabel= "", ylabel= "")

for ax, col in zip(axes[5,:], drug_names):
    ax.set_xlabel(col)

for ax, row in zip(axes[:,0], drug_names):
    ax.set_ylabel(row, rotation=90, size='large')

fig.tight_layout()
plt.show()
# %%
