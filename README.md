
# Introduction

This repository contains the scripts that are necessary for the data analysis performed in my master thesis "PhysiBoSS-Drugs: A personalized multiscale simulation framework that enables predicting drug synergies in cell populations”.

# Available scripts

## Convert PhysiBoSS-Drugs output files into a CSV file

The script **PhysiBoSS-Drugs_to_csv.py** which is based on the script **plot_time_course.py** from the Github repository: [https://github.com/migp11/tools4physicell](https://github.com/migp11/tools4physicell) enables the conversion of single-, double- and untreated drug simulation data that was obtained with PhysiBoSS-Drugs into a CSV file.
The output CSV file is stored in the output folder of tools4PhysiBoSS-Drugs. The CSV file is needed to perform the following data analysis and corresponding plots.

The CSV file used for the analysis in my master thesis can be found in the data folder with the name **LNCaP_simulation_data.csv**.

### Usage of PhysiBoSS-Drugs_to_csv.py

`PhysiBoSS-Drugs_to_csv.py --single data/single_data_folder --double data/double_data_folder --untreated data/untreated_data_folder`

## Ploting time course: plot_time_course.py

The script plot_time_course.py plot number of cells vs time grouping cell by phase (alive, necrotic, apoptotic) <br>
The color mapping can be easily customized to represent other cell-agent variables (eg. color mutants or other cell states)
	
~~~~
usage: plot_time_course.py [-h] [--format {physicell,physiboss}]
                           [--figout FIG_FNAME] [--csvout CSV_FNAME]
                           data_folder

Plot total cell grouped as Alive/Necrotic/Apoptotic vs Time

positional arguments:
  data_folder           folder were the data is stored

optional arguments:
  -h, --help            show this help message and exit
  --format {physicell,physiboss}
                        Format of the input data
  --figout FIG_FNAME    File name to save the plot
  --csvout CSV_FNAME    File name to store the summary table used for the plot	
~~~~

#### Examples

`plot_time_course.py output_test --format physiboss --figout physibos_time_plot.png`

`plot_time_course.py output_test --format physicell --figout physicell_time_plot.png`
