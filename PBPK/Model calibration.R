# Loading R packages
library(mrgsolve)
library(dplyr)
library(zoo)  # the package was loaded for the function na.locf 
library(ggplot2)
library(patchwork)
library(FME)
library(Metrics)
library(tidyr)
library(minpack.lm)
library(GenSA)

## Build mrgsolve-based PBPK Model
mod <- mcode("NanoPBPK.code", NanoPBPK.code,compile=TRUE)

# read the predicted parameters
data <- read.csv(file = 'Data_PBPK parameters.csv',fileEncoding = 'UTF-8-BOM')
#########################
Physic_pars<-c(
  #Cardiac Output and Blood flow
  QCC            = 16.5,   #L/h/kg^0.75, Cardio output,               (Brown, 1997)
  QLC            = 0.02,   #unitless, Fraction blood flow to liver,   (Brown, 1997, Table 23)
  QLuC           = 1,      #unitless, Fraction blood flow to lung,    (Brown, 1997, Table 23)
  QKC            = 0.091,  #unitless, Fraction blood flow to kidney,  (Brown, 1997, Table 23)
  QBrC           = 0.033,  #unitless, Fraction blood flow to brain,   (Brown, 1997, Table 23)
  QSC            = 0.011,  #unitless, Fraction blood flow to spleen,  (Lin, 2008; Davies and Morries, 1993)
  QMC            = 0.159,  #unitless, Fraction blood flow to muscle, (Brown, 1997, Table 23)

  #Tissue Volume
  VLC            = 0.055,  #unitless, Fraction liver tissue,          (Brown, 1997, Table 21)
  VLuC           = 0.007,  #unitless, Fraction lung tissue,           (Brown, 1997, Table 21)
  VKC            = 0.017,  #unitless, Fraction kidney tissue,         (Brown, 1997, Table 21)
  VBrC           = 0.017,  #unitless, Fraction brain tissue,          (Brown, 1997, Table 21)
  VSC            = 0.005,  #unitless, Fraction spleen tissue,         (Lin, 2008; Davies and Morries, 1993)
  VBldC          = 0.06,   #unitless, Fraction blood,                 (Chen, 2015)
  VPlasC         = 0.0355, #unitless, Fraction plasma,                (Davies and Morris, 1993
  VMC            = 0.384,  #unitless, Fraction muscle tissue,         (Brown, 1997, Table 21)

   #Blood volume fraction in organs and tissues
  BVL            = 0.31,   #unitless, Liver,                          (Brown, 1997, Table 30)
  BVBr           = 0.03,   #unitless, Brain,                          (Brown, 1997, Table 30)
  BVK            = 0.24,   #unitless, Kidney,                         (Brown, 1997, Table 30)
  BVS            = 0.17,   #unitless, Spleen,                         (Brown, 1997, Table 30)
  BVLu           = 0.5,    #unitless, lungs,                          (Brown, 1997, Table 30)
  BVM            = 0.04,   #unitless, muscle,                         (Brown, 1997, Table 30)
  BVR            = 0.04,   #unitless, Rest of body (assumed to equal to muscle), (Brown, 1997, Table 30)

  #//Partition coefficients(PC, tissue:plasma)
  PL             = 0.08,   #unitless, liver,                          (Lin, 2016)
  PK             = 0.15,   #unitless, kidney,                         (Lin, 2016)
  PBr            = 0.15,   #unitless, brain,                          (Lin, 2016)
  PS             = 0.15,   #unitless, spleen,                         (Lin, 2016)
  PLu            = 0.15,   #unitless, lungs,                          (Lin, 2016)
  PH             = 0.15,   #unitless, heart,                          (Lin, 2016)
  PM             = 0.15,   #unitless, muscle,                         (Lin, 2016)
  PR             = 0.15,   #unitless, rest of body,                   (Lin, 2016)
  
  #Membrane-limited permeability coefficient constants
  PALC           = 0.001,  #unitless, liver,                          (Lin, 2016)
  PABrC          = 0.000001, #unitless, brain,                        (Lin, 2016)
  PAKC           = 0.01,   #unitless, kidney,                         (Lin, 2016)
  PASC           = 0.15,   #unitless, spleen,                         (Lin, 2016)
  PALuC          = 0.001,  #unitless, lung,                           (Lin, 2016)
  PAMC           = 0.00005,#unitless, muscle,                         (Lin, 2016)
  PARC           = 0.00005,#unitless, rest of body,                   (Lin, 2016)

  #Endocytic parameters; RES represent endocytic/phagocytic cells
  KLRESrelease   = 0.0015,   #1/h,                liver,                          Release rate constant of phyagocytic cells
  KLRESmax       = 0.3,      #1/h,                liver,                          Maximmum uptake rate constant of phyagocytic cells
  KLRES50        = 48,       #h,                  liver,                          Time reaching half maximum uptake rate
  KLRESn         = 5,        #unitless,           liver,                          Hill coefficient

  KSRESrelease   = 0.001,    #1/h,                spleen,                         Release rate constant of phyagocytic cells
  KSRESmax       = 5,        #1/h,                spleen,                         Maximmum uptake rate constant of phyagocytic cells
  KSRES50        = 36,       #h,                  spleen,                         Time reaching half maximum uptake rate
  KSRESn         = 5,        #unitless,           spleen,                         Hill coefficient

  KKRESrelease   = 0.001,    #1/h,                kidney,                         Release rate constant of phyagocytic cells
  KKRESmax       = 0.12,     #1/h,                kidney,                         Maximmum uptake rate constant of phyagocytic cells
  KKRES50        = 48,       #h,                  kidney,                         Time reaching half maximum uptake rate
  KKRESn         = 5,        #unitless,           kidney,                         Hill coefficient

  KLuRESrelease  = 0.003,    #1/h,               lung,                           Release rate constant of phyagocytic cells
  KLuRESmax      = 0.085,    #1/h,               lung,                           Maximmum uptake rate constant of phyagocytic cells
  KLuRES50       = 48,       #h,                 lung,                           Time reaching half maximum uptake rate
  KLuRESn        = 5,        #unitless,          lung,                           Hill coefficient

  KMRESrelease   = 0.005,    #1/h,                muscle,                         Release rate constant of phyagocytic cells
  KMRESmax       = 0.4,      #1/h,                muscle,                         Maximmum uptake rate constant of phyagocytic cells
  KMRES50        = 48,       #h,                  muscle,                         Time reaching half maximum uptake rate
  KMRESn         = 5,        #unitless,           muscle,                         Hill coefficient
  
  KRRESrelease   = 0.005,    #1/h,                rest of body,                   Release rate constant of phyagocytic cells
  KRRESmax       = 0.4,      #1/h,                rest of body,                   Maximmum uptake rate constant of phyagocytic cells
  KRRES50        = 48,       #h,                  rest of body,                   Time reaching half maximum uptake rate
  KRRESn         = 5,        #unitless,           rest of body,                   Hill coefficient

  #Excretion parameters
  KbileC         =  0.00003,  #L/hr/kg^0.75,       Bile clearance
  KurineC        =  0.000003  #L/hr/kg^0.75,       Urine clearance
  )
##################
vpars <-c(
  QTC=0.01,
  VTC=0.006,
  BVT=0.02,
  PT             = 0.265,     #unitless, tumor,                          fitted
  PATC           = 0.01,      #unitless, tumor,                          fitted
  KTRESrelease   = 0.001,     #1/h,                tumor, Release rate constant of phyagocytic cells 
  KTRESmax       = 0.055,     #1/h,                tumor, Maximmum uptake rate constant of phyagocytic cells
  KTRES50        = 0.5,       #h,                  tumor, Time reaching half maximum uptake rate
  KTRESn         = 0.5        #unitless,           tumor, Hill coefficient
)
##################

#Read dataset
TvsC_Dat <- read.csv(file='EvaDat_1.csv',fileEncoding = 'UTF-8-BOM')
TvsC_Dat$ID<-na.locf(TvsC_Dat$ID) # replace na values with last no-na value
colnames(TvsC_Dat) <- c("ID","Time", "DETumor")

# Define the prediction function
# Predict concentration profile
pred.nano <- function(fixpars, Vpars, BW, Dose) {
  Vpars_exp <- exp(Vpars)
  DOSEiv <- Dose * BW
  
  ex <- ev(
    ID = 1, amt = DOSEiv, ii = 24, tinf = 0.005,
    addl = 0, cmt = "AV", replicate = FALSE
  ) + ev(
    ID = 1, amt = DOSEiv, ii = 24, tinf = 0.005,
    addl = 0, cmt = "ADOSE", replicate = FALSE
  )
  
  sim_end_time <-24
  tsamp <- tgrid(0, sim_end_time, 0.1)
  
  out <- mod %>%
    param(fixpars) %>%
    param(Vpars_exp) %>%
    update(atol = 1e-4, rtol = 1e-3, maxsteps = 50000) %>%
    mrgsim_d(data = ex, tgrid = tsamp)
  
  outdf <- data.frame(
    Time = out$time,
    DETumor = (out$Tumor / DOSEiv) * 100
  )
  
  return(outdf)
}

# ------------Define the input parameters-----------------

#Fault Parameter value
theta <- c(
  PT             = 0.04152393,
  PATC           = 0.00350176,
  KTRESrelease   = 0.5133574,
  KTRESmax       = 0.879352,
  KTRES50        = 17.80145031,
  KTRESn         = 3.1818864
  )

# Individual parameter bounds
param_bounds <- list(
  PT            = c(0.000202, 8.68),
  PATC          = c(0.0000049, 0.9),
  KTRESrelease  = c(0.000102, 12),
  KTRESmax      = c(0.001, 45.60),
  KTRES50       = c(0.00001, 180),
  KTRESn        = c(0.01, 10)
)

# Perturbation range
perturbation_range <- list(
  PT            = c(0.5, 2.0),   
  PATC          = c(0.3, 3.0),   
  KTRESrelease  = c(0.5, 2.0),    
  KTRESmax      = c(0.5, 2.0),   
  KTRES50       = c(0.3, 3.0),   
  KTRESn        = c(0.7, 1.4) 
)

# ---------Generate fixed parameters for each individual--------
Fix_func <- function(Physic_pars, indiv) {
  mod_pars <- Physic_pars
  mod_pars["BW"]  <- indiv$BW
  mod_pars["QTC"] <- 0.01
  mod_pars["VTC"] <- indiv$TW / indiv$BW
  mod_pars["BVT"] <- 0.02
  return(mod_pars)
}

# -------Objective function: cost-----------
Cost_nano <- function(pars, indiv, Physic_pars, w = "mean") {
  full_pars <- theta
  full_pars[names(pars)] <- exp(pars)
  fixed_pars <- Fix_func(Physic_pars, indiv)
  
  prediction <- pred.nano(
    fixpars = fixed_pars,
    Vpars = log(full_pars),
    BW = indiv$BW,
    Dose = indiv$Dose
  )
  
  modCost(model = prediction, obs = indiv$obs, x = "Time", weight = w)
}

#---------------------Parameter fitting---------------------------------
fit_single_id <- function(id, theta, Physic_pars, data, TvsC_Dat,
                          r2_threshold = 0.85, max_try =1, min_attempts = 1) {
  
  # ----------------- Step 1: Prepare Individual Data -----------------
  A_i <- data %>% filter(ID == id) %>%
    mutate(across(c(BW, TW, Dose), as.numeric)) %>%
    select(BW, TW, Dose)
  Obs_i <- TvsC_Dat %>% filter(ID == id) %>% select(Time, DETumor)
  indiv <- list(BW = A_i$BW / 1000, TW = A_i$TW / 1000, Dose = A_i$Dose, obs = Obs_i)
  
  # ----------------- Step 2: Bounds -----------------
  fit_param_names <- names(theta)
  lower_bound <- log(sapply(fit_param_names, function(p) param_bounds[[p]][1]))
  upper_bound <- log(sapply(fit_param_names, function(p) param_bounds[[p]][2]))
  
  # ----------------- Step 3: Gradient-Guided Perturbation -----------------
  generate_gradient_guided_init <- function(base_theta, prev_theta = NULL, prev_score = NULL, attempt = 1, max_try = 500) {
    step_scale <- min(1.0, 0.2 + log1p(attempt) / log1p(max_try))
    
    sapply(fit_param_names, function(p) {
      lower <- log(param_bounds[[p]][1])
      upper <- log(param_bounds[[p]][2])
      center <- log(base_theta[[p]])
      
      direction <- if (!is.null(prev_theta) && !is.null(prev_score)) {
        delta <- center - log(prev_theta[[p]])
        if (abs(delta) < 1e-6) delta <- runif(1, -1, 1)
        if (prev_score > 0) -sign(delta) else sign(delta)
      } else {
        sample(c(-1, 1), 1)
      }
      
      for (j in 1:200) {
        val <- center + direction * runif(1, 0.2, 1.5) * step_scale * (upper - lower) / 8
        if (val >= lower && val <= upper) return(val)
      }
      return(center)
    })
  }
  
  # ----------------- Step 4: Iterative Fit Loop -----------------
  best_score <- -Inf
  best_result <- NULL
  all_results <- list()
  attempt <- 1
  prev_theta <- NULL
  prev_score <- NULL
  p_init <- log(theta)
  
  while (attempt <= max_try) {
    cat("\U0001F501 Fitting ID =", id, "(Attempt", attempt, ") ...\n")
    
    fit <- tryCatch({
      modFit(
        f = function(pars) Cost_nano(pars, indiv = indiv, Physic_pars = Physic_pars),
        p = p_init,
        method = "Marq",
        lower = lower_bound,
        upper = upper_bound,
        control = nls.lm.control(nprint = 0)
      )
    }, error = function(e) {
      cat("Optimization error:", e$message, "\n")
      NULL
    })
    
    if (is.null(fit)) {
      attempt <- attempt + 1
      p_init <- generate_gradient_guided_init(theta, prev_theta, prev_score, attempt, max_try)
      next
    }
    #---Evaluation paraeter fitting
    final_theta <- theta
    final_theta[names(fit$par)] <- exp(fit$par)
    res <- Cost_nano(log(final_theta), indiv = indiv, Physic_pars = Physic_pars)
    
    PDat <- res$residuals %>%
      mutate(ID = id,
             Log.OBS = log(obs, 10),
             Log.PRE = log(mod, 10)) %>%
      filter(is.finite(Log.OBS) & is.finite(Log.PRE))
    
    R2   <- summary(lm(Log.OBS ~ Log.PRE, data = PDat))$r.squared
    MAPE <- mean(abs((PDat$obs - PDat$mod) / PDat$obs))
    
    if (is.na(R2) || is.na(MAPE)) {
      cat("R2 or MAPE is NA — skipping this attempt.\n")
      attempt <- attempt + 1
      p_init <- generate_gradient_guided_init(theta, prev_theta, prev_score, attempt, max_try)
      next
    }
    #Evaluated metrics
    score <- R2 * (1 - MAPE)
    is_success <- (R2 >= r2_threshold && MAPE <= 0.3)
    Status <- ifelse(is_success, "Success", "Tried")
    
    cat(sprintf("   → R2 = %.4f | MAPE = %.4f | Score = %.4f | Status = %s\n", 
                R2, MAPE, score, Status))
    cat("   → Parameters used:\n")
    for (param in names(final_theta)) {
      cat(sprintf("      %s = %.6f\n", param, final_theta[[param]]))
    }
    
    
    PDat <- PDat %>%
      mutate(Attempt = attempt, R2 = R2, MAPE = MAPE, Score = score)
    #Output dataframe
    current_result <- list(
      theta_df = data.frame(ID = id, Parameter = names(final_theta),
                            Value = final_theta, R2 = R2, MAPE = MAPE, Score = score, Status = Status),
      mod_obs = PDat,
      R2 = R2,
      Score = score
    )
    
    all_results[[length(all_results) + 1]] <- current_result
    
    if (is_success && attempt >= min_attempts) {
      cat("success criteria met, exiting early.\n")
      best_result <- current_result
      break
    }
    
    if (!is.na(score) && score > best_score) {
      best_score <- score
      best_result <- current_result
    }
    
    prev_theta <- final_theta
    prev_score <- -sum((PDat$obs - PDat$mod)^2, na.rm = TRUE)
    p_init <- generate_gradient_guided_init(final_theta, prev_theta, prev_score, attempt, max_try)
    attempt <- attempt + 1
  }
  
  if (is.null(best_result)) {
    theta_df <- data.frame(ID = id, Parameter = NA, Value = NA, R2 = NA, MAPE = NA, Score = NA, Status = "Fail")
    return(list(theta_df = theta_df, mod_obs = NULL, all_results = all_results))
  } else {
    return(list(
      theta_df = best_result$theta_df,
      mod_obs = best_result$mod_obs,
      all_results = all_results
    ))
  }
}
# -----------Run fitting for all IDs------------

ID_list <- c(532)
theta_df_list <- list()
mod_obs_list <- list()
mod_obs_all_attempts <- list()

for (id in ID_list) {
  result <- fit_single_id(id, theta, Physic_pars, data, TvsC_Dat)
  
  
  theta_df_list[[as.character(id)]] <- result$theta_df
  mod_obs_list[[as.character(id)]] <- result$mod_obs
  
  for (i in seq_along(result$all_results)) {
    mod_i <- result$all_results[[i]]$mod_obs
    if (!is.null(mod_i)) {
      mod_obs_all_attempts[[paste0("ID", id, "_Attempt", i)]] <- mod_i
    }
  }
}

# Combine all outputs
theta_all_df <- do.call(rbind, theta_df_list)
mod_obs_all <- do.call(rbind, mod_obs_list)
mod_obs_all_attempt_df <- do.call(rbind, mod_obs_all_attempts)

# View results
print(theta_all_df)
print(mod_obs_all)

