
# Introduction

This repository contains the scripts that are necessary to perform drug simulations with PhysiBoSS-Drugs and the data analysis of my master thesis "PhysiBoSS-Drugs: A personalized multiscale simulation framework that enables predicting drug synergies in cell populations”.

# Available scripts

## Before running simulations with PhysiBoSS-Drugs 

### Fitting GDSC data

The raw GDSC data ([http://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_public_raw_data_25Feb20.csv](http://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_public_raw_data_25Feb20.csv)) needs to be fitted to a multi-level mixed effect model (*Vis et al, 2016*). This can be done with the help of the R package *gdscIC50* ([https://github.com/CancerRxGene/gdscIC50](https://github.com/CancerRxGene/gdscIC50)). The script **fit_GDSC_data.R** contains the necessary code to obtain dose response curves and the CSV file that is necessary to integrate into PhysiBoSS-Drugs. The raw GDSC data has to be downloaded beforehand and stored to the data folder. <p>&nbsp;</p>

## Data-analysis after performing drug simulations

### Convert PhysiBoSS-Drugs output files into a CSV file

The script **PhysiBoSS-Drugs_to_csv.py** which is based on the script **plot_time_course.py** from the Github repository: [https://github.com/migp11/tools4physicell](https://github.com/migp11/tools4physicell) enables the conversion of single-, double- and untreated drug simulation data that was obtained with PhysiBoSS-Drugs into a CSV file.
The output CSV file will be stored in the output folder of tools4PhysiBoSS-Drugs. The CSV file is needed to perform the following data analysis and to obtain corresponding plots.

The CSV file used for the analysis in my master thesis can be found in the data folder with the name **LNCaP_simulation_data.csv**.

#### Usage of PhysiBoSS-Drugs_to_csv.py

`PhysiBoSS-Drugs_to_csv.py --single data/single_data_folder --double data/double_data_folder --untreated data/untreated_data_folder` <p>&nbsp;</p>

### Plotting heatmaps (Fig 2, Fig 3 and SuppMat, Fig 13, Fig 14, Fig 15)

For plotting the pair growth and synergy heatmaps (Fig 2 and 3 of the Main text) the two scripts **growth_pair_heatmap.py** and **synergy_pair_plot.py** can be used. The scripts require the specification of the path that leads to the CSV file containing the simulation data. The data used in my master thesis can be found in **"/data/LNCaP_simulation_data.csv**.

For plotting the complete growth and synergy heatmaps (SuppMat, Fig 13, Fig 14, Fig 15) the two scripts **full_growth_heatmap.py** and **full_synergy_plot.py** can be used. The same path specifications as above have to be done. Furthermore, the two drugs that should be analyzed in the heatmaps have to be specified within the script. <p>&nbsp;</p>

### Kruskal-Wallis

To test which of the drug combinations is significant a Kruskal-Wallis rank sum test can be performed with the script **kruskal-wallis.py**. The path to the simulation data and the used drugs have to be specified inside the script. <p>&nbsp;</p>

### Experimental validation

The experimental data used in my thesis is stored in **data/LNCaP_experimental_data.csv** Plots of the experimental data can be obtained with the script **plot_experimental_data.py**. The output will be saved in the output folder of the repository. \
To calculate the experimental AUCs the script **calc_experimental_AUC.py** can be used.
To validate the results additional drug simulations were performed for *Luminespib*, *Pictilisib* and *Selumetinib* with matching drug concentrations. The corresponding CSV files can be found in **data/validation_Luminespib_3_3uM.csv**,**data/validation_Pictilisib_10uM.csv** and **data/validation_Selumetinib_30uM.csv**. The growth indices for those simulations can be calculated with the script **validation_simulation_AUC.py**. <p>&nbsp;</p>

### Plotting the time course

The script to plot the timecourse of drug simulations as in SuppMat, Fig 12 is **timecourse_SuppMat_section7_fig12.py** which is based on the **plot_time_course.py** script from [https://github.com/migp11/tools4physicell](https://github.com/migp11/tools4physicell).