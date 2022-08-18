
# Corrected estimate, Mhat_adj, for specified publication bias and internal bias
# - dat: dataset, which must contain column of two-tailed p-values called "pval", point estimates "yi", variances "vi"
#   point estimates must be on log-RR scale if using g-transformation
# - cluster: indicator for clusters of point estimates, if needed
# - Ci: indicator for whether each study is internally biased
# - EB.affirm.obs: muB|A=1 as defined in paper
# - EB.nonaffirm.obs: muB|A=0 as defined in paper
# - eta: selection ratio
# - transformation: g-transformation to go from MULTIPLICATIVE bias factor to your E-value of choice
#      depending on type of bias (defaults to no transformation); see report_meta for more information
corrected_meta_mbma = function(dat,
                               cluster = 1:nrow(dat),
                               Ci,
                               
                               # sens params (on log scale)
                               EB.affirm.obs,
                               EB.nonaffirm.obs,
                               eta,
                               
                               # effect-size transformations for reporting
                               # see report_meta for info
                               transformation = function(x) x,
                               transformed.scale.name = "",
                               
                               # for pasting the analysis label
                               suffix = NA
) {
  
  dat$Ci = Ci
  dat$cluster = cluster
  
  # assumes favor.positive = TRUE or that estimates are already direction-flipped
  dat$affirm = dat$pval < 0.05 & dat$yi > 0
  
  # get P(A^*_i = 1 | C^*_i = 1) from P(A^*_i = 1 | C^*_i = 1, D^*_i = 1) 
  P.affirm.pub = mean( dat$affirm[ dat$Ci == 1 ] )
  # and P(A^*_i = 0 | C^*_i = 1):
  P.nonaffirm.pub = 1 - P.affirm.pub
  
  
  denom = P.affirm.pub + eta * P.nonaffirm.pub
  
  # a sample estimate of muB*
  # but note that if you set EB.affirm.obs = EB.nonaffirm.obs, MhatB = muB
  ( MhatB = (1/denom) * ( P.nonaffirm.pub * eta * EB.nonaffirm.obs +
                            P.affirm.pub * EB.affirm.obs ) )
  
  # c.f. bias in *published* studies (will be larger than MhatB):
  #  P.nonaffirm.pub * EB.nonaffirm.obs + P.affirm.pub * EB.affirm.obs
  
  # adjusted yi's using the estimated MhatB
  dat$yi.adj.est = dat$yi - dat$Ci * MhatB
  
  
  weights = rep( 1, length(dat$yi.adj.est) )
  # weight based on the affirm indicator of the *confounded* estimates
  weights[ dat$affirm == FALSE ] = eta
  
  # initialize a dumb (unclustered and uncorrected) version of tau^2
  # which is only used for constructing weights
  meta.re = rma.uni( yi = dat$yi.adj.est,
                     vi = dat$vi)
  t2hat.naive = meta.re$tau2 
  
  
  # fit weighted robust model
  meta.mbma = robu( yi.adj.est ~ 1,
                    studynum = cluster,
                    data = dat,
                    userweights = weights / (vi + t2hat.naive),
                    var.eff.size = vi,
                    small = TRUE )
  
  # sanity check: compare to direct calculation of user-weighted estimator
  expect_equal( as.numeric( meta.mbma$b.r ),
                sum( dat$yi.adj.est * ( weights / (dat$vi + t2hat.naive) ) ) / sum( weights / (dat$vi + t2hat.naive) ) )
  
  label = ifelse( is.na(suffix), "mbma", paste("mbma", suffix, sep = "-") )
  .res = report_meta(meta.mbma,
                     .mod.type = "robu",
                     .analysis.label = label,
                     .transformation = transformation,
                     .transformed.scale.name = transformed.scale.name)
  
  # add more info
  .res$eta_assumed = eta
  .res$EB.affirm.obs_assumed = EB.affirm.obs
  .res$EB.nonaffirm.obs = EB.nonaffirm.obs
  
  return(.res)
  
}


# E-value (i.e., g(muB*)) for a fixed choice of eta
# args: see corrected_meta above
# EB_grid_hi: largest value of EB to consider in grid search
# q: value to which to shift Mhat_adj
evalue_mbma = function(dat,
                       Ci,
                       eta,
                       EB_grid_hi = log(10),
                       q = 0,
                       
                       evalue.transf = function(B) B + sqrt(B^2 - B)  
) {
  
  #### Additive Bias Factor Required to Explain Away the CI ####
  func = function(.EB.shared) {
    
    # note that below does NOT transform Mhat in order to have it on same scale as q
    corrected = corrected_meta_mbma(dat = dat,
                                    cluster = dat$cluster,
                                    Ci = Ci,
                                    
                                    # sens params
                                    EB.affirm.obs = .EB.shared,
                                    EB.nonaffirm.obs = .EB.shared,
                                    
                                    eta = eta )
    
    return( abs( corrected$Mhat - q ))
  }
  
  opt = optimize( f = func,
                  interval = c(0, EB_grid_hi),
                  maximum = FALSE )
  
  
  # bias on ADDITIVE scale (e.g., log if yi's are log-RRs) required to shift Mhat to q
  ( EB_est = opt$minimum )
  
  # discrepancy between the corrected estimate and the s-value
  diff = opt$objective
  
  # if the optimal value is very close to the upper range of grid search
  #  AND we're still not very close to the target q,
  #  that means the optimal value was above eta_grid_hi
  if ( abs(EB_est - EB_grid_hi) < 0.0001 & diff > 0.0001 ) EB_est = paste(">", EB_grid_hi)
  
  
  
  #### Additive Bias Factor Required to Explain Away the CI ####
  func = function(.EB.shared) {
    
    # note that below does NOT transform Mhat in order to have it on same scale as q
    corrected = corrected_meta_mbma(dat = dat,
                                    cluster = dat$cluster,
                                    Ci = Ci,
                                    
                                    # sens params
                                    EB.affirm.obs = .EB.shared,
                                    EB.nonaffirm.obs = .EB.shared,
                                    
                                    eta = eta )
    
    return( abs( corrected$MLo - q ))
  }
  
  opt = optimize( f = func,
                  interval = c(0, EB_grid_hi),
                  maximum = FALSE )
  
  # bias on ADDITIVE scale (e.g., log if yi's are log-RRs) required to shift CI to q
  ( EB_ci = opt$minimum )
  
  # discrepancy between the corrected estimate and the s-value
  diff = opt$objective
  
  # if the optimal value is very close to the upper range of grid search
  #  AND we're still not very close to the target q,
  #  that means the optimal value was above eta_grid_hi
  if ( abs(EB_ci - EB_grid_hi) < 0.0001 & diff > 0.0001 ) EB_ci = paste(">", EB_grid_hi)
  
  
  #### Transform to E-value Scale ####
  eval_est = evalue.transf( exp(EB_est) )
  eval_ci = evalue.transf( exp(EB_ci) )
  
  #### Return ####
  return( data.frame(EB_est = EB_est,
                     EB_ci = EB_ci,
                     eval_est = eval_est,
                     eval_ci = eval_ci,
                     eta_assumed = eta) )
}












# INPUT/OUTPUT FNS ----------------------------------------------

# one or both dirs can be NA
my_ggsave = function(name,
                     .plot = last_plot(),
                     .width,
                     .height,
                     .results.dir = results.dir,
                     .overleaf.dir = overleaf.dir) {
  
  dirs = c(.results.dir, .overleaf.dir)
  dirIsNA = sapply(dirs, is.na)
  validDirs = dirs[ !dirIsNA ]
  
  
  for ( dir in validDirs ) {
    setwd(dir)
    ggsave( name,
            plot = .plot,
            width = .width,
            height = .height,
            device = "pdf" )
  }
}

# for reproducible manuscript-writing
# adds a row to the file "stats_for_paper" with a new statistic or value for the manuscript
# optionally, "section" describes the section of code producing a given result
update_result_csv = function( name,
                              section = NA,
                              value = NA,
                              print = FALSE ) {
  setwd(results.dir)
  
  new.rows = data.frame( name,
                         value = as.character(value),
                         section = as.character(section) )
  
  # to avoid issues with variable types when overwriting
  new.rows$name = as.character(new.rows$name)
  new.rows$value = as.character(new.rows$value)
  new.rows$section = as.character(new.rows$section)
  
  
  if ( "stats_for_paper.csv" %in% list.files() ) {
    .res = read.csv( "stats_for_paper.csv",
                     stringsAsFactors = FALSE,
                     colClasses = rep("character", 3 ) )
    
    # if this entry is already in the results file, overwrite the
    #  old one
    if ( all(name %in% .res$name) ) .res[ .res$name %in% name, ] = new.rows
    else .res = rbind(.res, new.rows)
  }
  
  if ( !"stats_for_paper.csv" %in% list.files() ) {
    .res = new.rows
  }
  
  write.csv( .res, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  # also write to Overleaf
  setwd(overleaf.dir)
  write.csv( .res, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  if ( print == TRUE ) {
    View(.res)
  }
}




# stands for "wipe results"
wr = function(){
  setwd(results.dir)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
  setwd(overleaf.dir)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
}

# stands for "view results"
vr = function(){
  setwd(results.dir)
  View( read.csv("stats_for_paper.csv") )
}


# GENERIC SMALL HELPERS ----------------------------------------------

# quick mean with NAs removed
meanNA = function(x){
  mean(x, na.rm = TRUE)
}

# quick median with NAs removed
medNA = function(x){
  median(x, na.rm = TRUE)
}

# quick length(unique) equivalent
uni = function(x){
  length(unique(x))
}


# ~ Handling strings -----------------

# return strings containing anything in pattern vector
stringsWith = function(pattern, x){
  # make regex expression 
  patterns = paste(pattern, collapse="|")
  x[ grepl(pattern = patterns, x = x)]
}
# stringsWith( pattern = c("dog", "cat"),
#  x = c("dogcat", "horse", "cat", "lion") )


# return indices of strings containing anything in pattern vector
whichStrings = function(pattern, x){
  patterns = paste(pattern, collapse="|")
  grepl(pattern = pattern, x = x)
}

names_with = function(.dat, .pattern) {
  names(.dat)[ grepl(pattern = .pattern, x = names(.dat) ) ]
}



# ~ Calculate simple stats -----------------

quick_ci = function( est, var ) {
  c( est - qnorm(.975) * sqrt(var),
     est + qnorm(.975) * sqrt(var) )
}

quick_pval = function( est, var ) {
  2 * ( 1 - pnorm( abs( est / sqrt(var) ) ) )
}


# ~ Formatting stats as strings -----------------

# round while keeping trailing zeroes
my_round = function(x, digits) {
  formatC( round( x, digits ), format='f', digits=digits )
}

format_CI = function( lo, hi, digits ) {
  paste( "[", my_round( lo, digits ), ", ", my_round( hi, digits ), "]", sep="" )
}

# round down to nearest integer
format_sval = function( sval, digits ) {
  if ( as.character(sval) == "--" ) return("Already NS")
  else if ( as.character(sval) == "Not possible" ) return("Not possible")
  else return( as.character( round(sval, digits) ) )
  #else return( floor(sval) )
}

format_pval = function(p) {
  if (p >= 0.01) return( my_round( p, 2 ) )
  if (p < 0.01 & p > 10^-5 ) return( formatC( p, format = "e", digits = 0 ) )
  if ( p < 10^-5 ) return("< 1e-05")
}

# ~ Meta-analysis fns -----------------

# nicely report a metafor or robumeta object with optional suffix to denote which model
# in case yi's have been transformed, applies function .transformation to put Mhat and its CI on interpretable scale (e.g., set .transformation = function(x) exp(x) if yi's are log-RRs)
# **Shat and its CI will NOT be transformed
report_meta = function(.mod,
                       .mod.type = "rma",  # "rma" or "robu"
                       .transformation = function(x) x,  
                       .transformed.scale.name = "",  # e.g., "RR" if you set .transformation = exp
                       .analysis.label = "" ) {
  
  if ( !is.null(.mod) ) {
    
    
    if ( .mod.type == "rma" ) {
      tau.CI = tau_CI(.mod)
      .res = data.frame( .transformation(.mod$b),
                         .transformation(.mod$ci.lb),
                         .transformation(.mod$ci.ub),
                         
                         .mod$pval,
                         
                         sqrt(.mod$tau2),
                         tau.CI[1],
                         tau.CI[2] )
    } 
    
    
    if ( .mod.type == "robu" ) {
      
      # catch possibility that tau.sq is NULL, which should only happen if
      #  you passed userweights to robu
      tau = ifelse( is.null(.mod$mod_info$tau.sq),
                    NA,
                    sqrt(.mod$mod_info$tau.sq) )
      
      .res = data.frame( .transformation(.mod$b.r),
                         .transformation(.mod$reg_table$CI.L),
                         .transformation(.mod$reg_table$CI.U),
                         
              
                         .mod$reg_table$prob,
                         
                         tau,
                         NA,
                         NA )
    } 
    
  } else {
    .res = data.frame( rep(NA, 6) )
  }
  
  # rename columns
  names(.res) = c("Mhat", "MLo", "MHi", "MPval",
                  "Shat", "SLo", "SHi")
  
  .res = .res %>% add_column(.before = 1, 
                             Analysis = .analysis.label)
  
  .res$MhatScale = .transformed.scale.name
  
  
  # names(.res) = paste( c("Mhat", "MLo", "MHi", "MPval",
  #                        "Shat", "SLo", "SHi"), .suffix, sep = "" )
  row.names(.res) = NULL
  
  #return( list(stats = .res) )
  
  return( .res )
}


