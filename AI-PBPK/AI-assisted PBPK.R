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
library(ggplot2)
library(openxlsx)

## Build mrgsolve-based PBPK Model
mod <- mcode("NanoPBPK.code", NanoPBPK.code,compile=TRUE)

# read the predicted parameters
#NP properties
data <- read.csv(file = 'data_total.csv',fileEncoding = 'UTF-8-BOM')

# Time concentration
TvsC_Dat <- read.csv(file='EvaDat_1.csv',fileEncoding = 'UTF-8-BOM') 
TvsC_Dat$ID<-na.locf(TvsC_Dat$ID) # replace na values with last no-na value
colnames(TvsC_Dat) <- c("ID","Time", "DETumor")

#PBPK NP parameters
Theta<- read.csv(file='data_total.csv',fileEncoding = 'UTF-8-BOM')
Thata_parameter<- as.data.frame(Theta)

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
# Define the prediction function
# Predict concentration profile and delivery efficiency
pred_DE <- function(fixpars, Vpars, BW, Dose, timepoint,pred=FALSE) {
  Vpars_exp <- exp(Vpars)
  DOSEiv <- Dose * BW
  
  ex <- ev(
    ID = 1, amt = DOSEiv, ii = 24, tinf = 0.005,
    addl = 0, cmt = "AV", replicate = FALSE
  ) + ev(
    ID = 1, amt = DOSEiv, ii = 24, tinf = 0.005,
    addl = 0, cmt = "ADOSE", replicate = FALSE
  )
  
  sim_end_time <-336
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
  
  if(pred) return(outdf)
  
  outdf<-cbind.data.frame(
    DE24  = outdf %>% filter(Time == 24) %>% select(DE24 = DETumor),
    DE168 = outdf %>% filter(Time == 168) %>% select(DE168 = DETumor),
    DEmax = outdf %>% filter(DETumor == max(DETumor)) %>% select(DEmax = DETumor)
    # DETlast = outdf %>%filter(Time == timepoint) %>%select(DETlast = DETumor)
  )

  return(outdf)
}

# ---------Generate fixed parameters for each individual--------
Fix_func <- function(Physic_pars, indiv) {
  mod_pars <- Physic_pars
  mod_pars["BW"]  <- indiv$BW
  mod_pars["QTC"] <- indiv$QTC
  mod_pars["VTC"] <- indiv$TW / indiv$BW
  mod_pars["BVT"] <- indiv$BVT
  return(mod_pars)
}

#Fault PBPK parameter value
theta <- c(
  PT             = 0.008859185,
  PATC           = 0.0107586996 ,
  KTRESrelease   = 0.7735230240,
  KTRESmax       = 0.3479386042,
  KTRES50        = 0.0007910517,
  KTRESn         = 0.102063846
)
#
TvsC_Dat <- read.csv(file='EvaDat_1.csv',fileEncoding = 'UTF-8-BOM')
TvsC_Dat$ID<-na.locf(TvsC_Dat$ID) # replace na values with last no-na value
colnames(TvsC_Dat) <- c("ID","Time", "DETumor")

# -------Objective function: cost-----------
Cost_nano <- function(ID, Physic_pars, w = "mean") {
  A_i <- data %>% filter(.data$ID == !!ID) %>%
    mutate(across(c(BW, TW, Dose,QTC, BVT), as.numeric)) %>%
    select(BW, TW, Dose,QTC, BVT)
  
  Obs_i <- TvsC_Dat %>% filter(.data$ID == !!ID) %>% select(Time, DETumor)
  indiv <- list(BW = A_i$BW / 1000, TW = A_i$TW / 1000, Dose = A_i$Dose, 
                QTC= A_i$QTC,BVT=A_i$BVT, obs = Obs_i)
  #Combine as the fixed parameters
  fixed_pars <- Fix_func(Physic_pars, indiv)
  
  #Get the tumor related parameters
  pars <- Thata_parameter %>%
    filter(.data$ID == !!ID) %>%
    select(PT, PATC, KTRESrelease, KTRESmax, KTRES50, KTRESn) %>%
    as.list() %>%                
    unlist(use.names = TRUE)     
  
  full_pars <- theta
  full_pars[names(pars)] <- pars
  
  prediction <- pred_DE(
    fixpars = fixed_pars,
    Vpars = log(full_pars),
    BW = indiv$BW,
    Dose = indiv$Dose,
    timepoint= timepoint,
    pred=TRUE
  )

  cc<-modCost(model = prediction, obs = indiv$obs, x = "Time", weight = w)
  
  DE<-pred_DE(fixpars = fixed_pars,
              Vpars = log(full_pars),
              BW = indiv$BW,
              Dose = indiv$Dose,
              timepoint= timepoint,
              pred=FALSE)
  
  return(list(ID=ID, DE = DE, residuals = cc$residuals))
}

#List set
ID_list <- 1:530
all_DE_list <- list()
all_residuals_list <- list()

# Loop for prediction of DE 
for (id in ID_list) {
  result <- tryCatch({
    Cost_nano(ID = id, Physic_pars = Physic_pars, w = "mean")
  }, error = function(e) {
    message(sprintf("Error at ID = %s: %s", id, e$message))
    return(NULL)
  })
  
  if (!is.null(result)) {
    DE_tmp <- result$DE
    if (!all(is.na(DE_tmp))) {
      DE_tmp$ID <- id
      all_DE_list[[as.character(id)]] <- DE_tmp
    }
    
    residuals_tmp <- result$residuals
    if (!all(is.na(residuals_tmp))) {
      residuals_tmp$ID <- id
      all_residuals_list[[as.character(id)]] <- residuals_tmp
    }
  }
}

# Combine the results
all_DE_df <- bind_rows(all_DE_list)
all_residuals_df <- bind_rows(all_residuals_list)

all_residuals_df<-all_residuals_df%>% mutate(OPR=mod/obs)

# 2-fold error analysis
n   <- all_residuals_df %>%summarise (count = n())
n_2 <- all_residuals_df %>%
  filter(OPR>=0.5 & OPR<=2) %>% summarise (count = n())
n_3 <- all_residuals_df %>%
  filter(OPR>=0.33 & OPR<=3) %>% summarise (count = n())

N2 <- (n_2$count [1]/n$count [1])*100
N3 <- (n_3$count [1]/n$count [1])*100

#Linear regression analysis

log_all_residuals_df<- all_residuals_df %>% 
  mutate(Log_obs = log10(obs),
         Log_mod = log10(mod))

log_all_residuals_df_clean <- log_all_residuals_df %>%
  filter(is.finite(Log_obs), is.finite(Log_mod))


lm_fit_all <- lm(Log_obs ~ Log_mod, data = log_all_residuals_df_clean )
summary(lm_fit_all)


p3 <- ggplot(log_all_residuals_df_clean, aes(x = Log_obs, y = Log_mod)) +
  geom_point(color = "black", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  annotate("text",
           x = min(log_all_residuals_df_clean$Log_mod, na.rm = TRUE),
           y = max(log_all_residuals_df_clean$Log_mod, na.rm = TRUE),
           label = paste0("R² = ", round(summary(lm(Log_mod ~ Log_obs, data = log_all_residuals_df_clean))$r.squared, 3)),
           hjust = -0.1, vjust = 1.2, size = 5, fontface = "bold") +
  labs(x = "Observed DE (log10)", y = "Predicted DE(log10)") +
  theme_minimal(base_size = 14)









# Export Excel files
write.xlsx(list(DE = all_DE_df, Residuals = all_residuals_df),
           file = "DE_results.xlsx")



