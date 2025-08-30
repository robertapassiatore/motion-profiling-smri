
# Convert t-statistic to partial r
t2r = function(t, df) {
  sign(t) * sqrt(t^2 / (t^2 + df))
}

# Fisher’s r-to-z transformation
r2z = function(r) {
  0.5 * log((1 + r) / (1 - r))
}

# z to Fisher’s r transformation
z2r <- function(z) {
  (exp(2 * z) - 1) / (exp(2 * z) + 1)
}



  
res.dir = "/Users/roberta.passiatore/Documents/Motion-Profiling/analyses/NCvsSCZ_GroupDifferences/results/"
setwd(res.dir)
datasets <- list.files(pattern=glob2rx("*regression_results_permutations_corticalROI*csv*"))

all.data = NULL
for (dataset in datasets) {
  #dataset = datasets[1]
  data <- read_csv(dataset)
  #data = data[-which(data$Analysis=='standard'),]
  dataset_name <- gsub(paste0("regression_results_permutations_corticalROI_"), "", basename(dataset))
  dataset_name <- gsub(".csv", "", basename(dataset_name))
  data$R = t2r(data$T_values, data$N-5)
  data$zscore = r2z(data$R)
  data$Site = strsplit(dataset_name,'_')[[1]][1]
  data$Session = strsplit(dataset_name,'_')[[1]][2]
  all.data = rbind(all.data, data)
}

# ---------------------------
# 1) Create and initialize the meta.data data frame
# ---------------------------
meta.data <- as.data.frame(unique(all.data$ROI))
colnames(meta.data)[1] <- 'ROI'
meta.data$zval  <- NA
meta.data$meanR  <- NA
meta.data$ci.lb <- NA
meta.data$ci.ub <- NA
meta.data$pval <- NA
meta.data$I2   <- NA
meta.data$H2   <- NA
meta.data$QE   <- NA
meta.data$QEp  <- NA

meta.data$session <- NA  # to store the session info

# Initialize an empty data frame to store the results
meta.data.res <- data.frame()

# ---------------------------
# 2) Define sites and sessions
# ---------------------------
sites <- unique(all.data$Site)
all.data.org = all.data
#all.data = all.data[-which(all.data$Site=='UCLA'),]
#sites = sites[c(1:4,6,7)]
session <- unique(all.data$Session)
session <- session[1:3]  # using first three sessions - no FE (only one site)

# ---------------------------
# 3) Loop through sessions and ROIs
# ---------------------------
for (i in 1:length(session)) {
  # Subset the data for the current session
  data.session <- all.data[all.data$Session == session[i], ]
  
  # Iterate over each ROI in meta.data
  for (r in 1:nrow(meta.data)) {
    # Subset the data for the current ROI
    data.sub <- data.session[data.session$ROI == meta.data$ROI[r], ]
    
    # If there's no data for this ROI in the current session, skip this iteration
    if (nrow(data.sub) == 0) next
    
    # Build a data frame for meta-analysis
    data.test <- data.frame(
      study    = data.sub$Site,
      cohens_d = data.sub$R,
      std.err  = abs(1-data.sub$SE)
    )
    
    # Perform meta-analysis using the metafor package
    res <- rma(yi = cohens_d, sei = std.err, data = data.test, method = "ML")
    # Two‐level: random intercept for each study
    res2 <- rma.mv(yi, vi,
                   random = ~ 1 | study,    # one random intercept per study
                   data   = dat,
                   method = "REML")
    
    
    # Populate the meta.data data frame with the meta-analysis results
    meta.data$zval[r]  <- res$zval
    meta.data$ci.lb[r] <- res$ci.lb
    meta.data$ci.ub[r] <- res$ci.ub
    meta.data$pval[r]  <- res$pval
    meta.data$I2[r]    <- res$I2
    meta.data$H2[r]    <- res$H2
    meta.data$QE[r]    <- res$QE
    meta.data$QEp[r]   <- res$QEp
  }
  
  # Assign the current session value to all rows in meta.data for this iteration
  meta.data$session <- session[i]
  
  # Append the updated meta.data for this session to meta.data.res
  meta.data.res <- rbind(meta.data.res, meta.data)
}

# View the combined results
#print(meta.data.res)

  
write.csv(meta.data.res, file=paste0('res.LinearEffect_allPC_corticalROI_metanalysis.csv'))


# MIXED EFFECT METANALYSIS
# ---------------------------
# 1) Create and initialize the meta.data data frame
# ---------------------------
meta.data <- as.data.frame(unique(all.data$ROI))
colnames(meta.data)[1] <- 'ROI'
meta.data$zval  <- NA
meta.data$meanR  <- NA
meta.data$ci.lb <- NA
meta.data$ci.ub <- NA
meta.data$pval <- NA
meta.data$I2   <- NA
meta.data$H2   <- NA
meta.data$QE   <- NA
meta.data$QEp  <- NA

meta.data$session <- NA  # to store the session info

# Initialize an empty data frame to store the results
meta.data.res <- data.frame()

# ---------------------------
# 2) Define sites and sessions
# ---------------------------
sites <- unique(all.data$Site)
all.data.org = all.data
#all.data = all.data[-which(all.data$Site=='UCLA'),]
#sites = sites[c(1:4,6,7)]
session <- unique(all.data$Session)
session <- session[1:3]  # using first three sessions - no FE (only one site)

# ---------------------------
# 3) Loop through sessions and ROIs
# ---------------------------
data.session = all.data[-which(all.data$Session=='FE'),]
# Iterate over each ROI in meta.data
  for (r in 1:nrow(meta.data)) {
    # Subset the data for the current ROI
    data.sub <- data.session[data.session$ROI == meta.data$ROI[r], ]
    
    # If there's no data for this ROI in the current session, skip this iteration
    if (nrow(data.sub) == 0) next
    
    # Build a data frame for meta-analysis
    data.test <- data.frame(
      study    = data.sub$Site,
      r = data.sub$R,
      session = data.sub$Session,
      std.err  = abs(1-data.sub$SE)
    )
    
    # Perform meta-analysis using the metafor package
    #res <- rma(yi = cohens_d, sei = std.err, data = data.test, method = "ML")
    # Two‐level: random intercept for each study
    res <- rma.mv(yi = r, vi = std.err, 
                  data = data.test,
                  random   = ~ study + session,
                  #random = ~ 1 | study,    # one random intercept per study
                  method = "ML")
    summary(res)
    
    # Populate the meta.data data frame with the meta-analysis results
    meta.data$zval[r]  <- res$zval
    meta.data$ci.lb[r] <- res$ci.lb
    meta.data$ci.ub[r] <- res$ci.ub
    meta.data$pval[r]  <- res$pval
    meta.data$I2[r]    <- res$I2
    meta.data$H2[r]    <- res$H2
    meta.data$QE[r]    <- res$QE
    meta.data$QEp[r]   <- res$QEp
    
    # Append the updated meta.data for this session to meta.data.res
    meta.data.res <- rbind(meta.data.res, meta.data)
}
  

  



# Three‐level: e.g. effects nested within studies, studies within labs
res3 <- rma.mv(yi, vi,
               random = list(~ 1 | lab,
                             ~ 1 | study),
               data   = dat,
               method = "REML")
summary(res3)


