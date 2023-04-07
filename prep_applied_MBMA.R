
# PRELIMINARIES ---------------------------------------------------------------

# This script uses renv to preserve the R environment specs (e.g., package versions.)
library(renv)
# run this if you want to reproduce results using the R environment we had:
# renv::restore()

# data-wrangling packages
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyverse)
library(fastDummies)
# meta-analysis packages
library(metafor)
library(robumeta)
library(PublicationBias)
library(MetaUtility)
# other
library(here)
library(xtable)
library(testthat)

# run this only if you want to update the R environment specs
# renv::snapshot()


# set working directories
code.dir = here()

data.dir = str_replace_all( string = here(),
                                    pattern = "Code",
                                    replacement = "Data" ) 

results.dir = str_replace_all( string = here(),
                            pattern = "Results",
                            replacement = "Data" ) 


# get helper fns
setwd(code.dir)
source("helper_applied_MBMA.R")

# no sci notation
options(scipen=999)



# EXAMPLE: MEAT CONSUMPTION -----------------------------------------

# using this meta-analysis because has a clear randomized vs. nonrandomized indicator
setwd("~/Dropbox/Personal computer/Reference sheets/Library of prepped example meta-analyses/MRM")

d = fread("mathur_awr_prepped.csv")

d = d %>% filter( exclude.main == FALSE )
expect_equal( nrow(d), 100 )
expect_equal( length(unique(d$authoryear)), 34 )

table(d$randomized)

d$yi = d$logRR
d$vi = d$varlogRR
d$sei = sqrt(d$vi)
d$pval = 2 * ( 1 - pnorm( abs( d$yi/sqrt(d$vi) ) ) )
d$affirm = d$yi > 0 & d$pval < 0.05


# no apparent clustering from the forest plot
d$cluster = d$authoryear


d = d %>% add_column(.before = 1,
                     meta.name = "mathur")

# save it
setwd(data.dir)
fwrite(d, "mathur_prepped.csv")




# EXAMPLE: ATRIAL FIBRILLATION & COGNITIVE IMPAIRMENT -----------------------------------------


# Figure 1
# k = 14

# authors reported no pub bias based on Egger's
RR = c(1.70,
       1.05,
       1.56,
       8.10,
       1.14,

       2.88,
       4.01,
       1.10,
       0.90,
       1.08,
       1.56,
       1.38,
       1.09,
       1.14)

hi = c(2.63,
       2.20,
       1.92,
       34.53,
       1.78,

       6.58,
       8.74,
       3.02,
       1.62,
       1.46,
       1.74,
       1.73,
       2.20,
       1.26)

study.type = c( rep("Cross-sectional", 5),
                rep("Prospective cohort", 9) )


d = scrape_meta(est = RR,
                hi = hi)

d = d %>% rename(vi = vyi)
d$sei = sqrt(d$vi)
d$pval = 2 * ( 1 - pnorm( abs( d$yi/sqrt(d$vi) ) ) )
d$affirm = d$yi > 0 & d$pval < 0.05

d$effect.measure = "log-rr"

# no apparent clustering from the forest plot
d$cluster = 1:nrow(d)

d$study.type = study.type

d = d %>% add_column(.before = 1,
                     meta.name = "kalantarian_all")

# sanity check: reproduce reported analysis
# reported (abstract): RR: 1.40 [CI, 1.19 to 1.64]).
mod = rma.uni(yi = d$yi,
              vi = d$vi,
              method  = "DL")
expect_equal( as.numeric( round( exp(mod$b), 2 ) ), 1.40 )
expect_equal( as.numeric( round( exp(mod$ci.lb), 2 ) ), 1.19 )
# works :)

# save it
setwd(data.dir)
fwrite(d, "kalantarian_all_prepped.csv")





