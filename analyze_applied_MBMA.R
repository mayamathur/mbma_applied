

# Goals:

# For fixed eta:
# - Get E-value for est and CI. Can't you just use existing PublicationBias::svalue after pre-shifting estimates? :)
# - Sanity check: Compare that to the E-value obtained using Eq. 4.1 in Tier 2 (SAPB est / lambda thing)

# Next: Try the dementia meta-analysis you saved :)

# Idea for presentation in paper:
#  - For fixed eta taken from SAPB-E
#  - For worst-case eta

# Notation in results files:
#  - "unadj" suffix in analysis name means unadjusted for confounding, but could be adjusted for pub bias

# To do:
# - If using AWR in paper, note that for that analysis, I used modelweights = HIER, which slightly affects naive meta-analysis results


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
library(weightr)
# other
library(here)
library(xtable)
library(testthat)
library(EValue)

# run this only if you want to update the R environment specs
# renv::snapshot()


# parameters to toggle
# should we redo the somewhat slow line plots?
redo.plots = TRUE

# set working directories
code.dir = here()

data.dir = str_replace_all( string = here(),
                            pattern = "Code",
                            replacement = "Data" ) 

results.dir = str_replace_all( string = here(),
                               pattern = "Code",
                               replacement = "Results" ) 

overleaf.dir = "/Users/mmathur/Dropbox/Apps/Overleaf/Multiple-bias meta-analysis Overleaf (MBMA)/R_objects"


# get helper fns
setwd(code.dir)
source("helper_applied_MBMA.R")

# no sci notation
options(scipen=999)



# ANALYZE BOTH METAS ---------------------------------------------------------------

meta.names = c("kalantarian_all",
               "mathur")

# to run only one
#meta.names = "mathur"

for ( i in 1:length(meta.names) ) {
  
  
  meta.name = meta.names[i]
  
  # get prepped data
  setwd(data.dir)
  d = fread( paste(meta.name, "prepped.csv", sep = "_"))
  
  # fixed eta
  # c.f. SAPB-E: for top med metas, mean etahat=1.02 and Q95 = 1.62
  eta = 4
  
  # ~ Set up parameters specific to this meta-analysis -----------
  
  if ( meta.name == "kalantarian_all" ) {
    
    # effect size transformation
    z.to.r = FALSE
    take.exp.inv = FALSE
    take.exp = TRUE
    transf.name = "exp"  # just used for reporting
    
    # threshold for effect sizes
    q = log(1)
    
    # decide which studies are confounded
    Ci = rep(1, nrow(d))
    
    # # vector of selection probabilities to use for all specifications
    # # more dense at the very small ones
    # eta.vec = rev( seq(1, 200, 0.25) )
    # #eta.vec = c( 20, rev(seq(1,15,1)) )
    # # as a list
    # el = as.list( eta.vec )
  }
  
  
  if ( meta.name == "mathur" ) {
    
    # effect size transformation
    z.to.r = FALSE
    take.exp.inv = FALSE
    take.exp = TRUE
    transf.name = "exp"  # just used for reporting
    
    # threshold for effect sizes
    q = log(1)
    
    # decide which studies are confounded
    Ci = ( d$randomized == FALSE )
    
    # # vector of selection probabilities to use for all specifications
    # # more dense at the very small ones
    # eta.vec = rev( seq(1, 200, 0.25) )
    # #eta.vec = c( 20, rev(seq(1,15,1)) )
    # # as a list
    # el = as.list( eta.vec )
  }
  
  
  # ~ Transformation for meta-analysis -----------
  # by default, identity transformation
  transf = function(x) x
  if ( z.to.r == TRUE ) transf = function(x) z_to_r(x)
  if ( take.exp.inv == TRUE ) transf = function(x) exp_inv(x)
  if ( take.exp == TRUE ) transf = function(x) exp(x)
  
  
  # ~ Naive model ----------- 
  
  ( meta.naive = robu( yi ~ 1, 
                       data = d, 
                       studynum = cluster,
                       var.eff.size = vi,
                       small = TRUE ) )
  
  # from AWR (note use of hierarchical model weights):
  # ( meta.rob = robu( logRR ~ 1,
  #                    data = d,
  #                    studynum = as.factor(authoryear),
  #                    var.eff.size = varlogRR,
  #                    modelweights = "HIER",
  #                    small = TRUE) )
  
  
  res = report_meta(meta.naive,
                    .mod.type = "robu",
                    .analysis.label = "naive",
                    .transformation = transf,
                    .transformed.scale.name = transf.name)
  
  
  cat("\nFlag1")
  
  # ~ Worst-case SAPB without confounding adjustment -----------
  
  # fit weighted robust model
  meta.worst = robu( yi ~ 1,
                     studynum = cluster,
                     data = d[ d$affirm == FALSE, ],
                     var.eff.size = vi,
                     small = TRUE )
  
  
  new.row = report_meta(meta.worst,
                        .mod.type = "robu",
                        .analysis.label = "worst-unadj",
                        .transformation = transf,
                        .transformed.scale.name = transf.name)
  res = bind_rows(res, new.row)
  
  
  # ~ SAPB (without confounding adjustment) -----------
  
  meta.SAPB.row = corrected_meta_mbma(dat = d,
                                      cluster = d$cluster,
                                      Ci = Ci,
                                      
                                      # sens params
                                      EB.affirm.obs = 0,
                                      EB.nonaffirm.obs = 0,
                                      
                                      eta = eta,
                                      
                                      # effect-size transformations for reporting
                                      # see report_meta for info
                                      transformation = transf,
                                      transformed.scale.name = transf.name,
                                      
                                      suffix = "unadj")
  
  res = bind_rows(res, meta.SAPB.row)
  
  
  # **sanity check vs. existing R package:
  # MBMA should agree with SAPB when there's no confounding
  meta.SAPB.check = PublicationBias::corrected_meta(yi = d$yi,
                                                    vi = d$vi,
                                                    cluster = d$cluster,
                                                    eta = eta,
                                                    model = "robust",
                                                    favor.positive = TRUE)
  
  expect_equal( meta.SAPB.row$Mhat, transf(meta.SAPB.check$est) )
  expect_equal( meta.SAPB.row$MLo, transf(meta.SAPB.check$lo) )
  expect_equal( meta.SAPB.row$MHi, transf(meta.SAPB.check$hi) )
  expect_equal( meta.SAPB.row$MPval, meta.SAPB.check$pval )
  
  
  
  # ~ MBMA: fixed eta, EB.affirm.obs, EB.nonaffirm.obs -----------
  
  ( meta.mbma.row = corrected_meta_mbma(dat = d,
                                        cluster = d$cluster,
                                        Ci = Ci,
                                        
                                        # sens params
                                        EB.affirm.obs = log(1.5),
                                        EB.nonaffirm.obs = log(1.1),
                                        
                                        # EB.affirm.obs = 0,
                                        # EB.nonaffirm.obs = 0,
                                        eta = eta,
                                        
                                        # effect-size transformations for reporting
                                        # see report_meta for info
                                        transformation = transf,
                                        transformed.scale.name = transf.name) )
  
  
  res = bind_rows(res, meta.mbma.row)
  
  cat("\nFlag2")
  
  
  # ~ MBMA: E-values for fixed etas -----------
  
  # start new results df for E-values
  res2 = evalue_mbma(dat = d,
                     Ci = Ci,
                     eta = eta,
                     EB_grid_hi = log(10),
                     q = 0)
  
  
  
  # and with eta = 1
  new.row = evalue_mbma(dat = d,
                        Ci = Ci,
                        eta = 1,
                        EB_grid_hi = log(10),
                        q = 0)
  res2 = res2 %>% add_row(new.row)
  
  # sanity check: should be close (but maybe not equal) to the E-value above for eta = 1
  evalue( est = RR( res$Mhat[ res$Analysis == "naive" ] ) )
  
  # **sanity check: 
  # if all studies are confounded, E-value for estimate should just be equal to Mhat from 
  #  SAPB since lambda = 1
  if ( all(Ci == 1) ) {
    expect_equal( exp( res2$EB_est[ res2$eta_assumed == eta ] ),
                  res$Mhat[ res$Analysis == "mbma-unadj" ],
                  tol = 0.001 )
    
    expect_equal( exp( res2$EB_ci[ res2$eta_assumed == eta ] ),
                  res$MLo[ res$Analysis == "mbma-unadj" ],
                  tol = 0.001 )
    
    
  }
  # for many more sanity checks, see 2022-7-3 check theory vs. R
  
  
  cat("\nFlag3")
  
  
  # ~ MBMA: E-values for worst-case eta -----------
  
  # for this one, can just directly apply E-value to worst-case estimate on RR scale
  temp = evalue( RR(res$Mhat[res$Analysis == "worst-unadj"] ),
                 RR(res$MLo[res$Analysis == "worst-unadj"] ) )
  
  # not putting Eb_est, etc. for this one because we're not reporting on that scale anyway
  res2 = res2 %>% add_row( EB_est = NA, EB_ci = NA, eval_est = temp["E-values", "point"], eval_ci = temp["E-values", "lower"], eta_assumed = Inf )
  
  
  
  
  # ~ Subadditivity/superadditivity plots as in R01 ---------------------------------------------------------------
  
  if ( meta.name == "mathur" & redo.plots == TRUE ) {
    
    
    # plotting dataframe
    .eta_running = c( seq(1, 20, 1) )  # "running" to indicate this is the running var for x-axis, not the eta used in analysis
    .group = c("None (uncorrected)", "Confounding only", "Publication bias only", "Both")
    .EB.obs.scen = c("EB.obs.same", "EB.obs.different")
    
    dp = expand.grid( .eta_running,
                      .group,
                      .EB.obs.scen )
    names(dp) = c("eta_running", "group", "EB_obs_scen")
    
    # eta to use in analysis depends on group
    dp$eta_assumed = dp$eta_running
    dp$eta_assumed[ dp$group %in% c("None (uncorrected)", "Confounding only") ] = 1
    
    # EB.obs settings to use in analysis depends on group
    conf_scens = c("Confounding only", "Both")
    dp$EB.nonaffirm.obs = 0
    dp$EB.affirm.obs = 0
    
    dp$EB.nonaffirm.obs[ dp$group %in% conf_scens & dp$EB_obs_scen == "EB.obs.same" ] = log(1.25)
    dp$EB.affirm.obs[ dp$group %in% conf_scens & dp$EB_obs_scen == "EB.obs.same" ] = log(1.25)
    
    dp$EB.nonaffirm.obs[ dp$group %in% conf_scens & dp$EB_obs_scen == "EB.obs.different" ] = log(1)
    dp$EB.affirm.obs[ dp$group %in% conf_scens & dp$EB_obs_scen == "EB.obs.different" ] = log(1.25*2)
    
    # dp$EB.nonaffirm.obs[ dp$group %in% conf_scens & dp$EB_obs_scen == "EB.obs.same" ] = log(2)
    # dp$EB.affirm.obs[ dp$group %in% conf_scens & dp$EB_obs_scen == "EB.obs.same" ] = log(2)
    # 
    # dp$EB.nonaffirm.obs[ dp$group %in% conf_scens & dp$EB_obs_scen == "EB.obs.different" ] = log(1)
    # dp$EB.affirm.obs[ dp$group %in% conf_scens & dp$EB_obs_scen == "EB.obs.different" ] = log(4)
    
    # sanity check
    dp %>% group_by(group, EB_obs_scen) %>%
      summarise( mean(exp(EB.affirm.obs)),
                 mean(exp(EB.nonaffirm.obs)) )
    
    
    dp = suppressMessages( dp %>% rowwise() %>%
                             mutate( corrected_meta_mbma(dat = d,
                                                         cluster = d$cluster,
                                                         Ci = Ci,  # this is a global var
                                                         
                                                         # sens params
                                                         # these vars are from dp, not d
                                                         EB.affirm.obs = EB.affirm.obs,
                                                         EB.nonaffirm.obs = EB.nonaffirm.obs,
                                                         
                                                         eta = eta_assumed,
                                                         
                                                         # effect-size transformations for reporting
                                                         # see report_meta for info
                                                         transformation = transf,
                                                         transformed.scale.name = transf.name) ) )
    # for plotting joy
    dp$EB_obs_scen_pretty = NA
    dp$EB_obs_scen_pretty[ dp$EB_obs_scen == "EB.obs.same" ] = "(A) Assume affirmatives and nonaffirmatives moderately confounded"
    dp$EB_obs_scen_pretty[ dp$EB_obs_scen == "EB.obs.different" ] = "(B) Assume affirmatives highly confounded; nonaffirmatives unconfounded"
    
    
    # ~~ Version #1: Mhat from each analysis type -----------
    
    # set line colors
    levels(dp$group)
    colors = c("black", "#ebae34", "#3499eb", "red")
    
    # for choosing plot limits
    min(dp$Mhat)
    ymin = 1
    
    max(dp$Mhat)
    ymax = 1.25
    
    p1 = ggplot( data = dp,
                 # needs to be eta_running rather than eta_assumed because latter is always 1 when not correcting pub bias
                 aes( x = eta_running,  
                      y = Mhat,
                      color = group ) ) +
      
      geom_hline( color = "black", lty = 2, yintercept = 1) +
      geom_line(size=1) +
      
      
      scale_y_continuous( limits = c(ymin, ymax),
                          breaks = seq(ymin, ymax, 0.05),
                          trans = "log10") +
      
      scale_x_continuous( limits = c(1, max(dp$eta_running)),
                          breaks = seq(1, max(dp$eta_running), 2)) +
      
      ylab("Corrected pooled estimate (RR)") +
      xlab( bquote( bold("Assumed") ~ bold(eta) ) ) +
      
      scale_color_manual( values = colors ) +
      guides( color = guide_legend("Biases corrected: ") ) +
      
      theme_classic() + 
      
      facet_grid( ~ EB_obs_scen_pretty ) +
      
      # base_size controls all text sizes; default is 11
      # https://ggplot2.tidyverse.org/reference/ggtheme.html
      theme_bw(base_size = 18) +
      theme( text = element_text(face = "bold"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             legend.position="bottom" ) 
    
    
    
    
    
    
    # ~~ Version #2: Ratio of Mhats -----------
    
    
    dp2 = dp %>% filter( group %in% c("Publication bias only", "Both") )
    
    
    dp2 = dp2 %>% group_by( eta_assumed, EB_obs_scen ) %>%
      mutate( Mhat_ratio = Mhat[ group == "Publication bias only" ] /
                Mhat[ group == "Both" ] ) %>%
      # only need one row per eta
      filter( !duplicated(eta_assumed) ) %>%
      select(eta_assumed,
             EB_obs_scen,
             Mhat_ratio)
    
    min(dp2$Mhat_ratio)
    ymin = 1
    
    max(dp2$Mhat_ratio)
    ymax = 1.09
    
    # ratio of Mhat when correcting for both biases vs. publication bias only
    
    p2 = ggplot( data = dp2,
                 aes( x = eta_assumed,
                      y = Mhat_ratio,
                      color = EB_obs_scen) ) +
      
      geom_line( size = 1.1) +
      
      
      geom_hline( color = "gray",
                  yintercept = 1,
                  lty = 0) +
      
      scale_y_continuous( limits = c(ymin, ymax),
                          breaks = seq(ymin, ymax, 0.01),
                          trans = "log10") +
      
      scale_x_continuous( limits = c(1, max(dp$eta_running)),
                          breaks = seq(1, max(dp$eta_running), 2)) +
      
      ylab("Ratio of Mhat corrected for pub bias / Mhat corrected for both (RR scale)") +
      xlab(bquote(eta)) +
      
      guides( color = guide_legend("Estimate type") ) +
      
      theme_classic()
    
    # ~ Write Plots -----------
    
    my_ggsave(name = paste(meta.name, "multiple_Mhats_plot.pdf", sep = "_"),
              .plot = p1,
              .width = 16,
              .height = 8)
    
    my_ggsave(name = paste(meta.name, "Mhat_ratio_plot.pdf", sep = "_"),
              .plot = p2,
              .width = 10,
              .height = 8)
    
  }  # end "if ( meta.name == "mathur" & redo.plots == TRUE )"
  
  # ~ Write Results Tables -----------
  
  fwrite(res, 
         paste(meta.name, "_corrected_metas.csv", sep = "") ) 
  
  fwrite(res2, 
         paste(meta.name, "_evalues.csv", sep = "") ) 
  
  # ~~ Write One-Off Stats ---------------------------
  setwd(results.dir)
  update_result_csv( name = paste( meta.name, " k" ),
                     value = nrow(d) )
  
  update_result_csv( name = paste( meta.name, " k affirm" ),
                     value = sum(d$affirm == TRUE) )
  
  update_result_csv( name = paste( meta.name, " k nonaffirm" ),
                     value = sum(d$affirm == FALSE) )
  
  update_result_csv( name = paste( meta.name, " k perc nonaffirm" ),
                     value = mean(d$affirm == FALSE)*100 )
  
  if ( meta.name == "mathur" ) {
    update_result_csv( name = paste( meta.name, " k randomized" ),
                       value = sum(d$randomized == TRUE) )
    
    update_result_csv( name = paste( meta.name, " k nonrandomized" ),
                       value = sum(d$randomized == FALSE) )
  }
  
  
  # lambdas
  # @CHECK THESE! Especially naive tau issue
  d$wi = 1/(d$vi + res$Shat[ res$Analysis == "naive" ]^2)
  
  if ( meta.name == "mathur" ) {
    lambdaAc = sum( d$wi[d$affirm == FALSE & d$randomized == FALSE] ) / sum( d$wi[ d$affirm == FALSE ] ) 
    lambdaA = sum( d$wi[d$affirm == TRUE & d$randomized == FALSE] ) / sum( d$wi[ d$affirm == TRUE ] ) 
  }
  
  if ( meta.name == "kalantarian_all" ) {
    # 1 becuase all studies confounded
    lambdaAc = 1
    lambdaA = 1 
  }
  
  update_result_csv( name = paste( meta.name, "lambdaAc" ),
                     value = lambdaAc )
  
  
  update_result_csv( name = paste( meta.name, "lambdaA" ),
                     value = lambdaA )
  
  
  # r ratio
  r = sum( d$wi[ d$affirm == TRUE ] ) / sum( d$wi[ d$affirm == FALSE ] )
  
  update_result_csv( name = paste( meta.name, "r ratio" ),
                     value = r )
  
  
  # attenuation % for uncorrected vs. fixed eta
  atten_perc = 100*( (res$Mhat[res$Analysis == "naive"] - res$Mhat[res$Analysis == "mbma-unadj"]) / res$Mhat[res$Analysis == "naive"] )
  update_result_csv( name = paste( meta.name, "attenutation perc" ),
                     value = atten_perc )
  
  # save all corrected metas
  update_result_csv( name = paste( meta.name, "Analysis=", res$Analysis, "Mhat" ),
                     value = res$Mhat )
  
  update_result_csv( name = paste( meta.name, "Analysis=", res$Analysis, "MLo" ),
                     value = res$MLo )
  
  update_result_csv( name = paste( meta.name, "Analysis=", res$Analysis, "MHi" ),
                     value = res$MHi )
  
  update_result_csv( name = paste( meta.name, "Analysis=", res$Analysis, "MPval" ),
                     value = res$MPval )
  
  
  # save all E-values
  update_result_csv( name = paste( meta.name, "eta=", res2$eta, "E-value est" ),
                     value = res2$eval_est )
  
  update_result_csv( name = paste( meta.name, "eta=", res2$eta, "E-value CI" ),
                     value = res2$eval_ci )
  
} # end loop over both metas






