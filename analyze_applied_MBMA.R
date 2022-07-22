

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



# ANALYZE THEM :) ---------------------------------------------------------------

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
  
  
  
  
  # ~ Save Results Locally -----------
  
  # ~~ One-off stats ---------------------------
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
  
  

  
  fwrite(res, 
         paste(meta.name, "_corrected_metas.csv", sep = "") ) 
  
  fwrite(res2, 
         paste(meta.name, "_evalues.csv", sep = "") ) 
  
  
  
  
  
  
  
  
  
  
  
  
  
  # . ----
  # . ----
  # UNMODIFIED CODE FROM SAPB-T ANALYSIS_GENERAL:  --------------------------------------------------------------
  # 
  # # Standard publication bias methods
  # # Egger test
  # Egger = regtest( x = d$yi,
  #                  vi = d$vi,
  #                  predictor = "sei")
  # 
  # # selection model
  # tryCatch({
  #   ( m1 = weightfunct( effect = d$yi,
  #                       v = d$vi,
  #                       steps = c(0.025, 1),
  #                       table = TRUE ) )
  # })
  # 
  # # S-Values
  # # to shift to 0
  # S0 = svalue( yi = d$yi,
  #              vi = d$vi,
  #              q = 0,
  #              clustervar = d$cluster,
  #              # always TRUE because we've flipped if necessary
  #              favor.positive = TRUE,
  #              model = "robust",
  #              return.worst.meta = TRUE )
  # 
  # # to shift to q
  # Sq = svalue( yi = d$yi,
  #              vi = d$vi,
  #              q = q,
  #              clustervar = d$cluster,
  #              # always TRUE because we've flipped if necessary
  #              favor.positive = TRUE,
  #              model = "robust")
  # 
  # 
  # 
  # #################### ~~ SAVE RESULTS #################### 
  # 
  # ##### ~~ Table for Paper #####
  # new.row = data.frame( Meta = meta.name,
  #                       k = nrow(d),
  #                       EstNaive = paste( my_round( estNaive, 2 ), format_CI( loNaive, hiNaive, 2 ), sep = " " ),
  #                       EstWorst = paste( my_round( estWorst, 2 ), format_CI( loWorst, hiWorst, 2 ), sep = " " ),
  #                       S0 = format_sval( S0$stats$sval.est, digits = 0 ),
  #                       S0.CI = format_sval( S0$stats$sval.ci, digits = 0 ),
  #                       Sq = format_sval( Sq$sval.est, digits = 0 ),
  #                       Sq.CI = format_sval( Sq$sval.ci, digits = 0 ),
  #                       Egger.pval = format_pval(Egger$pval) )
  # 
  # 
  # if ( i > 1 & exists("res") ) res = rbind( res, new.row ) else res = new.row
  # 
  # 
  # ##### One-Off Stats for Paper #####
  # 
  # update_result_csv( name = paste( meta.name, " k" ),
  #                    value = nrow(d) )
  # 
  # update_result_csv( name = paste( meta.name, " k affirm" ),
  #                    value = sum(d$affirm == TRUE) )
  # 
  # update_result_csv( name = paste( meta.name, " k nonaffirm" ),
  #                    value = sum(d$affirm == FALSE) )
  # 
  # update_result_csv( name = paste( meta.name, " k" ),
  #                    value = nrow(d) )
  # 
  # update_result_csv( name = paste( meta.name, "estNaive" ),
  #                    value = my_round( estNaive, 2 ) )
  # 
  # update_result_csv( name = paste( meta.name, "estNaive lo" ),
  #                    value = my_round( loNaive, 2 ) )
  # 
  # update_result_csv( name = paste( meta.name, "estNaive hi" ),
  #                    value = my_round( hiNaive, 2 ) )
  # 
  # update_result_csv( name = paste( meta.name, "estNaive pval" ),
  #                    value = format.pval( pvalNaive, eps = 0.0001 ) )
  # 
  # 
  # 
  # 
  # #################### ~~ LINE PLOT #################### 
  # 
  # # needed later even if we're not redoing line plot
  # axis.font.size = 16
  # axis.label.size = 20
  # 
  # if ( redo.line.plot == TRUE ) {
  #   
  #   # get estimates at each value
  #   res.list = lapply( el, 
  #                      function(x) {
  #                        
  #                        cat("\n Working on eta = ", x)
  #                        
  #                        return( corrected_meta( yi = d$yi, 
  #                                                vi = d$vi,
  #                                                eta = x,
  #                                                clustervar = d$cluster,
  #                                                model = "robust",
  #                                                favor.positive = TRUE) )
  #                      } )
  #   
  #   
  #   re.rob.ests = as.data.frame( do.call( "rbind", res.list ) )
  #   
  #   # save results because lapply above is slow
  #   # because each column is secretly a list, impeding write.csv
  #   re.rob.ests = as.data.frame( apply(re.rob.ests, 2, unlist) )
  #   
  #   ##### ~~ Make Plot #####
  #   
  #   # simplify breaks a little compared to eta
  #   breaks = c(200, 150, 100, 50, 40, 30, 20, 10, 5)
  #   axis.font.size = 16
  #   axis.label.size = 20
  #   
  #   if (meta.name == "Li") ylabel = "Corrected estimate (HR)"
  #   if (meta.name == "Ali") ylabel = "Corrected estimate (BMD % change)"
  #   
  #   p.rob = ggplot( ) +
  #     # "q" line
  #     geom_hline( yintercept = transf(q), lty = 2) +
  #     
  #     # point ests will be the same for worst-case, so doesn't matter
  #     geom_hline( yintercept = estWorst, color = "red", lty = 2 ) +
  #     
  #     geom_ribbon( data = re.rob.ests, aes( x = eta, ymin = transf(lo), ymax = transf(hi) ), fill = "black", alpha = .1) +
  #     
  #     geom_line( data = re.rob.ests, aes( x = eta, y = transf(est) ), color = "black", lwd = 1.2) +
  #     
  #     xlab( bquote( "Severity of hypothetical publication bias" ~ (eta) ) ) +
  #     ylab( ylabel ) + 
  #     
  #     #xlab( bquote( eta ) ) +
  #     #ylab( bquote(  hat(mu)[eta]^"'"  ) ) +
  #     
  #     #scale_x_continuous( limits = c(1, max(breaks)), breaks = breaks ) +
  #     
  #     theme_classic() + 
  #     
  #     theme(axis.title = element_text(size=axis.label.size),
  #           axis.text = element_text(size=axis.font.size) ) 
  #   
  #   
  #   # if (meta.name == "Li 2018") p.rob = p.rob + scale_y_continuous( limits = c( 0.8, 1.4 ),
  #   #                                                            breaks = seq( 0.8, 1.4, 0.1 ) )
  #   
  #   p.rob
  #   # LOOKS GREAT! 
  # }
  # 
  # 
  # 
  # #################### ~~ SIGNIFICANCE FUNNEL PLOT #################### 
  # 
  # 
  # if ( meta.name == "Li 2018" ) {
  #   xmin = -0.6  # these aren't being used
  #   xlabel = "Point estimate (log-HR; signs recoded)"
  # }
  # 
  # if (meta.name == "Alibhai 2017" ){
  #   xmin = 0
  #   xlabel = "Point estimate (BMD % change)"
  # }
  # 
  # p.funnel = suppressWarnings( significance_funnel( yi = d$yi,
  #                                                   vi = d$vi,
  #                                                   xmin = xmin,
  #                                                   xlab = xlabel,
  #                                                   favor.positive = TRUE,
  #                                                   est.N = estWorst,  # on analysis scale, not transformed to original scale
  #                                                   est.all = estNaive,
  #                                                   alpha = 0.05,
  #                                                   plot.pooled = TRUE ) )
  # 
  # p.funnel = p.funnel + 
  #   theme(axis.title = element_text(size=axis.label.size),
  #         axis.text = element_text(size=axis.font.size) ) 
  # 
  # 
  # #################### ~~ SAVE ALL THE THINGS #################### 
  # 
  # # plots
  # if ( redo.line.plot == TRUE ) ggsave( filename = paste( meta.name, "_robu_plot.png", sep = "" ),
  #                                       device = "png",
  #                                       p.rob,
  #                                       width = 9,
  #                                       height = 6)
  # 
  # ggsave( filename = paste( meta.name, "_funnel.png", sep = "" ),
  #         device = "png",
  #         p.funnel,
  #         width = 6,
  #         height = 4)
  # 
  # # also save straight to Overleaf
  # setwd(overleaf.dir)
  # if ( redo.line.plot == TRUE ) ggsave( filename = paste( meta.name, "_robu_plot.png", sep = "" ),
  #                                       device = "png",
  #                                       p.rob,
  #                                       width = 9,
  #                                       height = 6)
  # 
  # ggsave( filename = paste( meta.name, "_funnel.png", sep = "" ),
  #         device = "png",
  #         p.funnel,
  #         width = 6,
  #         height = 4)
  # 
  
  
}  # end huge for-loop over meta-analyses




