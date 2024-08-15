# Lab Flammability Testing

## A GitHub Repository for: 

## A Community Led Approach to Plant Flammability Testing for Landscaping Defensible Space

## Kristina Fauss, Joe Celebrezze, Indra Boving, Robert Fitch, Rachel Dye, and Max Moritz

--------------------------------

## Introductory Statement
This repository is meant for the storage and sharing of data, scripts, figures, and mixed effects model result tables related to the paper titled *A Community Led Approach to Plant Flammability Testing for Landscaping Defensible Space* by Kristina Fauss, Joe Celebrezze, Indra Boving, Robert Fitch, Rachel Dye, and Max Moritz which showcases a participatory research approach to laboratory flammability testing and recommends fire-safe plants for landscaping in the wildland-urban-interface (WUI) in Santa Barbara county and surrounding regions. Otherwise, the analyses and data wrangling presented here could be used as a framework for future studies investigating similar questions such as the ongoing work in the Moritz Fire Lab at UCSB.

--------------------------------

## Table of Contents

[Breakdown of Folders](https://github.com/celebrezze/flam-methods-comparison#breakdown-of-folders)

[- Raw Data](https://github.com/celebrezze/flam-methods-comparison#raw-data)

[- GIS](https://github.com/celebrezze/flam-methods-comparison#raw-data)

[- Processed Data](https://github.com/celebrezze/flam-methods-comparison#processed-data)

[- Scripts](https://github.com/celebrezze/flam-methods-comparison#scripts)

[- Figures](https://github.com/celebrezze/flam-methods-comparison#figures)

[- Mixed Effects Model Selection Tables](https://github.com/celebrezze/flam-methods-comparison#mixed-effects-model-selection-tables)

[Metadata](https://github.com/celebrezze/flam-methods-comparison#metadata)

[Contact Information](https://github.com/celebrezze/flam-methods-comparison#contact-information)

--------------------------------

## Breakdown of Folders

### Data:
The **data** folder consists of three subfolders - **raw-data**, **processed-data**, and **GIS**. For more information regarding specific column names, see metadata

#### Raw Data:
The **raw-data** folder consists of four subfolders:

  **flamm**: includes the raw data for laboratory flammability testing, *burn_samples_flamm.csv*
  
  ----
  
  **plant-traits**: includes multiple datasets which describe various plant traits as shown below
  
  *ARTCAL_LMA.csv*: data necessary to calculate the average leaf mass per area for *Artemesia californica*
  
  *cup_weights.csv*: the weights of the tins used for the LFM measurements and corresponding numbers
  
  *leaf_area_flamm.csv*: raw output from scanning leaves
  
  *leaf_data_flamm.csv*: leaf thickness and necessary masses to calculate LFM
  
  *stem_data_flamm.csv*: necessary raw data to calculate stem-specific plant traits
  
  *stemleaf_massratio_flamm.csv*: necessary data to calculate stem-to-mass and leaf-to-mass ratios
  
  ----
  
  **survey**
  
  *survey_addresses.csv*: addresses of survey respondents
  
  *survey_data.csv*: compilation of all necessary data to analyze survey responses
  
  ----
  
  **thermocouplers**: raw time series data for temperature and heat flux from the datalogger; labelled with date (YYYYMMDD)
  
#### Processed Data:
The **processed-data** folder consists of datasets manipulated/wrangled at some stage from the **raw-data**, primarily in data wrangling scripts (*flamm_data_wrangling.Rmd* and *plant_traits_data_wrangling.Rmd*). See the data wrangling scripts for more information on how we cleaned raw data and produced the necessary processed dataframes.

### GIS:
The **GIS** folder includes various geospatial data including rasters and geodatabases which identify the extent of wildland-urban-interface in the study area (Santa Barbara county) to report how important it is to conduct community-focused flammability testing in Santa Barbara county .

### Scripts:
The **scripts** folder includes scripts for all of the code we used to wrangle data, complete analyses, and design tables and figures for the main body of the paper, the supplementary index, and for exploratory analyses. The scripts are numbered in a logical order which follows the order presented in the paper. Further details regarding each of the 9 main scripts follow:

  *1.1_flamm_data_wrangling.Rmd*: this takes the datasets from the **flamm** subfolder of the **raw-data** folder and cleans them up so that they could be combined into one main dataset for further analyses. 
  
  *1.2_plant_traits_data_wrangling.Rmd*: this takes the datasets from the **plant-traits** subfolder of the **raw-data** folder, cleans them up, combines them for futher analyses **and** sets up clean dataframes which combine flammability measurements with plant traits which were used in the bulk of other analyses.
  
  *2.1_exploratory_analyses.Rmd*: this involves misc. exploratory analyses which we used to inform our expectations for future analyses, understand patterns in our data, and identify cases where further data cleaning was necessary.
  
  *2.2_sample_weight_exploration.Rmd*: this involves a more pointed investigation into the importance of sample weight (or sample mass, dry weight, and wet weight) in driving relationships between plant traits and flammability and compares different ways to account for differences in sample weight.
  
  *3_survey.Rmd*: this involves all analyses and figures made associated with the community survey data.
  
  *4_interspecific_differences.Rmd*: this focuses on analyses and figures concerned with interspecific differences in flammability metrics and plant traits; more specifically, it contains statistical tests, visualizations and summary tables that relate to interspecific differences.
  
  *5.1_MEM_selection.Rmd*: this involves all code necessary for the linear mixed effects model selection relied upon in the manuscript.
  
  *5.2_model_table_examine.Rmd*: this examines model selection tables to select top-performing models.
  
  *6_conceptual_figures.Rmd*: this involves relating flammability results to the results of the community survey. More specifically, it includes code necessary to produce flammability triangle plots which investigate interspecific differences in combustibility, consumability, and sustainability scores and it includes scatterplots which compare species 'desirability scores' to flammability scores.
  
  **extra-analyses**: other analyses which were either not relied upon in the final product or only minorly relied upon, typically to inform other analyses, were placed in this folder
  
  **python_scripts**: linear mixed model selections were scrutinized throughout the development of this project, and model selections were investigated using both R/RStudio and Python; the scripts in this folder are code involved in running MEM selections in python

### Python Notebook Outputs
The **python-nb-outputs** folder involves outputs from the MEM selections run on Python.

### Model Summary Tables
The **model_summary_tables** folder involves summary tables from MEM selections run on either Python or R.

### Figures:
The **figures** folder includes all figures included in the main text of the paper and the supplementary index, as well as figures we did not end up presenting (mostly placed in the *extra-figures* folder).

## Metadata:
This is located in *METADATA.Rmd* and *METADATA.html* (made from knitting *METADATA.Rmd*).

## Contact Information

Kristina Fauss*: kfauss@ucsb.edu

Joe Celebrezze: celebrezze@ucsb.edu

**correspondence*
