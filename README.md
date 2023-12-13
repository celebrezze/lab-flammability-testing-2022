# Lab Flammability Testing

## A GitHub Repository for: 

## Landscaping in fire-prone neighborhoods: laboratory flammability testing in response to community survey responses

## Kristina Fauss, Joe Celebrezze, Indra Boving, Robert Fitch, Rachel Dye, and Max Moritz

--------------------------------

## Introductory Statement
This repository is meant for the storage and sharing of data, scripts, figures, and mixed effects model result tables related to the paper titled *Landscaping in fire-prone neighborhoods: laboratory flammability testing in response to community survey responses* by Kristina Fauss, Joe Celebrezze, Indra Boving, Robert Fitch, Rachel Dye, and Max Moritz which showcases a partipatory research approach to laboratory flammability testing and recommends fire-safe plants for landscaping in the wildland-urban-interface (WUI) in Santa Barbara county and surrounding regions. Otherwise, the analyses and data wrangling presented here could be used as a framework for future studies investigating similar questions such as the ongoing work in the Moritz Fire Lab at UCSB.

--------------------------------

## Table of Contents

[Breakdown of Folders](https://github.com/celebrezze/flam-methods-comparison#breakdown-of-folders)

[- Raw Data](https://github.com/celebrezze/flam-methods-comparison#raw-data)

[- Processed Data](https://github.com/celebrezze/flam-methods-comparison#processed-data)

[- Scripts](https://github.com/celebrezze/flam-methods-comparison#scripts)

[- Figures](https://github.com/celebrezze/flam-methods-comparison#figures)

[- Mixed Effects Model Selection Tables](https://github.com/celebrezze/flam-methods-comparison#mixed-effects-model-selection-tables)

[Metadata](https://github.com/celebrezze/flam-methods-comparison#metadata)

[Contact Information](https://github.com/celebrezze/flam-methods-comparison#contact-information)

--------------------------------

## Breakdown of Folders

### Data:
The **data** folder consists of two subfolders, **raw-data** and **processed-data**. For more information regarding specific column names, see metadata

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
The **processed-data** folder consists of datasets manipulated/wrangled at some stage from the **raw-data**, primarily in data wrangling scripts. 

### Processed Data:
The **processed-data** folder consists of datasets manipulated/wrangled at some stage from the *raw-data*. This folder primarily includes datasets of different iterations that include data in which both methods were used simultaneously. See the metadata for more information. Otherwise, see *1.data_wrangling_methods.Rmd* (for the combined datasets as described in the metadata or *data_wrangling_SEKI.Rmd* (in *extra-analyses*; for *seki_flam_data_all.csv*) for details on data wrangling.

### Scripts:
The **scripts** folder includes scripts for all of the code we used to wrangle data, complete analyses, and design tables and figures for the main body of the paper, the supplementary index, and for exploratory analyses (which are primarily located in the *extra-analyses* folder inside of the **scripts** folder). The scripts are numbered in a logical order which follows the order presented in the paper. Further details regarding each of the 6 main scripts follow:

  *1.data_wrangling_methods.Rmd*: this takes the datasets from the *raw-data* folder and cleans them up so that they could be combined into one main dataset for further analyses. It removes species with less than 6 ignitions for either of the methods (leaving *Adenostoma fasciculatum*, *Ceanothus megacarpus*, *Arctostaphylos patula*, and *Ceanothus cordulatus*), removes rows with NA values in certain variables, and moves any instances of manual ignitions (i.e., after 7 minutes elapsed, we manually ignited the samples by lifting them into the propane-fueled pilot flame) into a specific dataset, otherwise removing them from the bulk of datasets.
  
  *2.literature_review.Rmd*: this involves all analyses and figures relating to the literature review. This includes the map labelled Figure 2a in the paper. IMPORTANT NOTE: For this map, we used global ignitions data readily available from the Global Fire Atlas through ORNL DAAC, Distributed Active Archive Center for Biogeochemical Dynamics (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1642). This data was **not** stored in our GitHub due to the large file size. Due to this, data must be locally requested from the website and then proper urls must be manually inputted into the code by whoever may be running the script for it to run properly.
  
  *3.1.water_content_vs_flam.Rmd*: this involves code necessary for figures that compare water content (or dry weight) to flammability metrics and accompanying mixed effects model selections and statistical tests. It includes a wide variety of iterations to look at the data including iterations not included in the paper or supplementary materials.
  
  *3.2.mixed_effects_model_selections.Rmd*: this involves the primary mixed effects model selection discussed in the paper and accompanying mixed effects model summary tables.
  
  *4.interspecific_differences.Rmd*: this focuses on the section of the supplementary index concerned with interspecific differences and involves many iterations of similar analyses (a lot of which were not utilized); otherwise, it contains statistical tests, visualizations and summary tables that relate to interspecific differences.
  
  *5.PCA.Rmd*: this involves all code necessary for the principal component analyses included in the paper.
  
  *extra-analyses*: as previously alluded to, any exploratory analyses or scripts which were improved upon or elaborated on by the main 6 scripts described above were placed in the *extra-analyses* folder. This folder includes analyses not mentioned above such as variance decomposition, classification and regression trees, segmented regressions, using the flammability index developed in Essaghi et. al. 2017, and an investigation into manual ignitions (mostly for *Ceanothus cordulataus*). Importantly, it also contains the *data_wrangling_SEKI.Rmd* file dedicated to wrangling the *SEKI.flammability.csv* into a more usable format in *seki_flam_data_all.csv*.
  
### Figures:
The **figures** folder includes all figures included in the main text of the paper and the supplementary index, as well as figures we did not end up presenting (mostly placed in the *extra-figures* folder). All main and supplementary figures were labelled appropriately with FigX. or FigSX. preceding the description of the figure. Note that *Fig1.methods.images.png* is not made in any script, but instead consists of two pictures taken by Indra Boving and Joe Celebrezze.

### Mixed Effects Model Selection Tables:
These are placed in the **mem-model-selection** folder and informed our conclusions regarding this analysis.

## Metadata:
This is located in *METADATA.Rmd* and *METADATA.html* (made from knitting *METADATA.Rmd*).

## Contact Information

Joe Celebrezze*: celebrezze@ucsb.edu

Indra Boving: bovingi@ucsb.edu

**correspondence*
