
# load data
library(dplyr)
library(ggplot2)
library(gdscIC50)
gdsc_raw_data_org <- read.csv("./data/GDSC2_public_raw_data_25Feb20.csv")


# remove failed drug treatments and missing drugs
gdsc_raw_data <- removeFailedDrugs(gdsc_raw_data_org)
gdsc_raw_data <- removeMissingDrugs(gdsc_raw_data)

normalized_gdsc_raw_data <- normalizeData(gdsc_raw_data,
                                         trim = T,
                                         neg_control = "NC-1",
                                         pos_control = "B")


# this dataframe now has a single row for each treated well on the plate in the data
# each well has a normalized_intensity value - a cell viability normalized to the control wells on the plate
head(normalized_gdsc_raw_data, 2)

# the nlme model used to fit the GDSC data needs every drug treatment to be compared on a single scale no matter
# what concentrations were used for a particular drug (this scale is set with the maximum concentration at 9 
# because original GDSC data used a 9 point 2 fold dilution for each drug treatment)

scaled_gdsc_raw_data <- setConcsForNlme(normalized_gdsc_raw_data, group_conc_ranges = F)

unique(subset(normalized_gdsc_raw_data, DRUG_ID_lib == 1510 & dose == "D1",
              select =  c("DRUGSET_ID", "lib_drug", "dose", "CONC")))

set_concs_test <- setConcsForNlme(normalized_gdsc_raw_data, group_conc_ranges = T) 

unique(subset(set_concs_test, DRUG_ID_lib == 1510,
              select =  c("lib_drug", "dose", "CONC", "maxc", "x")))

set_concs_test <- setConcsForNlme(normalized_gdsc_raw_data, group_conc_ranges = F) 
unique(subset(set_concs_test, DRUG_ID_lib == 1510, select =  c("lib_drug", "dose", "CONC", "maxc", "x")))


nlme_data <- prepNlmeData(scaled_gdsc_raw_data,
                          cl_id = "COSMIC_ID",
                          drug_specifiers = c("DRUG_ID_lib", "maxc"))

nlme_model <- fitModelNlmeData(nlme_data, isLargeData = F)

nlme_stats <- calcNlmeStats(nlme_model, nlme_data)
head(data.frame(nlme_stats), 2)

# select the relevant cell lines of the nlme_stats and filter out the important drugs
prostate_drugs = nlme_stats %>% 
                  filter(CL == 924100 | CL == 924105 | CL == 905935 | CL == 907788 | CL == 905934 | CL ==1299075) %>% 
                  select(xmid, scal, CL, maxc, DRUG_ID_lib, IC50) %>% 
                  filter(DRUG_ID_lib == 1924 | DRUG_ID_lib == 1912 | DRUG_ID_lib == 1032 | DRUG_ID_lib == 1168 | DRUG_ID_lib == 2047 | DRUG_ID_lib == 1559 | DRUG_ID_lib == 1372 | DRUG_ID_lib == 1736 | DRUG_ID_lib == 1058 | DRUG_ID_lib == 1560 | DRUG_ID_lib == 2043)
write.csv(prostate_drugs, "output/fitted_dose_response_curves.csv")
