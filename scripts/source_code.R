# Source Code

# packages
library(tidyverse)
library(ggsci) #for scale_color_lancet
library(Ternary) #for ternary plots
library(cowplot) #to arrange plots w/ plot_grid
library(ggpubr) #to arrange plots
library(lubridate) #to handle dates
library(hms) #to handle times
library(ggfortify) #to add things to ggplots
library(ggbiplot) #for nicer PCA plots
library(sjPlot) #for nice tables of eigenvalues
library(factoextra) #for k-means clustering
library(GGally) #for ggcorr()
library(psych) #for pairs.panels() to look at correlations
library(gam) #for generalized additive models
library(scales) #adjusting transparency (alpha) amongst other uses
library(kableExtra) #to make nice tables
library(remef) #for 'remef'
library(lme4) #for linear mixed effects models
library(caret) #for optimization routine
library(magrittr)
library(janitor) #for data wrangling
library(bayestestR) #for area_under_curve
#library(multcompView) # to display Tukey-Kramer results
library(factoextra) #for k-means clustering
library(emmeans) # for tukey-kramer tests (w/ random effects)
library(performance) # for multicollinearity check
library(car) # for qqPlot
library(gtools) # for combination()
library(glmmLasso) # specific to running lasso on generalized linear mixed models
library(glmnet) # 'best' package for running lasso on generalized linear models; also, has good function for cross-validation process
#library(HH)
library(devtools)
library(likert) # for likert scale
#library(hrbrthemes)


# defining functions
group_by = dplyr::group_by
summarise = dplyr::summarise
plot_grid = cowplot::plot_grid
ggbiplot = ggbiplot::ggbiplot
rename = dplyr::rename
select = dplyr::select
alpha = scales::alpha

# main dataset
flamm.df <- read.csv(here('data', 'processed-data', 'main_dataset.csv'))