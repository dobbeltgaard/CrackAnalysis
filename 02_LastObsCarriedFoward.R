
rm(list = ls())

### LOCF per defect level
### Objective: To analyse crack growth in terms of track behavior.
### Assumption: assuming last observation is valid until new measurement, creating regular temporal grid (1month resolution)

d = read.csv("crack_data_v2.csv")

library(dplyr)
library(tidyr)
library(lubridate)
library(data.table)
library(zoo)
d$Date <- as.Date(d$Date)
d <- d %>% mutate(Date_month = floor_date(Date, "month"))
constant_cols = c("ID", "in_straight","in_transition_curve","in_curve","covariates_complete","rail_age", "R200", "R260", "R350HT", "rail_1m_weight","traffic_direction", 
                  "Curve", "Overheight", "MGT_min", "MGT_max", "EMGT_min", "EMGT_max", "Line_speed", "Turnout_indicator","rail_left","BTR", "Track", "From", "To", "Profile","Fastening","detection_method","crack_counts", "BTR_first2digits")
measurement_cols = c("wear_h_mu", "wear_h_sigma", "wear_h_rms", "wear_h_p01","wear_h_p99", "wear_w_mu", "wear_w_sigma","wear_w_rms", "wear_w_p01", "wear_w_p99",
                  "corrugation_mu", "corrugation_sigma","corrugation_rms", "corrugation_p01", "corrugation_p99", "inclination_mu","inclination_sigma", "inclination_rms",
                  "inclination_p01", "inclination_p99","twist2m_mu", "twist2m_sigma","twist2m_rms", "twist2m_p01","twist2m_p99",  "twist3m_mu","twist3m_sigma","twist3m_rms",
                  "twist3m_p01","twist3m_p99","gauge_mu","gauge_sigma","gauge_rms","gauge_p01","gauge_p99","cant_mu","cant_sigma","cant_rms","cant_p01","cant_p99","align_horiz_d0_mu",
                  "align_horiz_d0_sigma","align_horiz_d0_rms","align_horiz_d0_p01","align_horiz_d0_p99","align_horiz_d1_mu","align_horiz_d1_sigma","align_horiz_d1_rms","align_horiz_d1_p01",
                  "align_horiz_d1_p99","align_vert_d0_mu","align_vert_d0_sigma","align_vert_d0_rms","align_vert_d0_p01","align_vert_d0_p99","align_vert_d1_mu","align_vert_d1_sigma","align_vert_d1_rms","align_vert_d1_p01","align_vert_d1_p99")

first_non_na <- function(x) {
  out <- na.omit(x)
  if (length(out) == 0) {
    if (is.character(x)) return(NA_character_)
    if (is.numeric(x))  return(NA_real_)
    if (is.integer(x))  return(NA_integer_)
    if (is.logical(x))  return(NA)
    return(NA)
  } else {return(out[1])}
}



interp_data <- function(data, ID, constant_cols, measurement_cols, order = 0,cutoff_date = NULL){
  data=data[data$ID == ID, ] #look at specific ID
  original_cols <- names(data) #for later reordering
  data <- data %>% mutate(Date_month = floor_date(Date, "month")) #regridding
  data <- data %>%group_by(Date_month) %>%summarise(across(all_of(measurement_cols), ~ mean(.x, na.rm = TRUE)),across(setdiff(names(data), c("Date_month", measurement_cols)),first_non_na),.groups = "drop") #summarizing over repeated months
  data = data  %>% complete(Date_month = seq.Date(min(Date_month), max(Date_month), by = "month"))  %>% arrange(Date_month) #make time grid
  if(!is.null(cutoff_date)) data=data %>% filter(Date_month >= cutoff_date) #filter out old observations
  data <- data %>% fill(all_of(constant_cols), .direction = "updown") # LOCF/BOCF for constant columns
  if(order == 0){ data = data %>% fill(all_of(measurement_cols), .direction = "down") }  #LOCF
  else if(order == 1){for (col in measurement_cols) { data[[col]] <- zoo::na.approx(data[[col]], x = data$Date_month, na.rm = FALSE) }} # Linear interpolation for measurement columns}
  data <- data %>% filter(if_all(all_of(measurement_cols), ~ !is.na(.))) #trunacte by available measurement observations
  removed_idx = which(data$State %in% c("Removed","No_defect"))+1
  removed_idx = removed_idx[removed_idx <= NROW(data) & is.na(data$crack_size[removed_idx])] #if end obs, cutoff so we're not out of bounds
  data$crack_size[removed_idx] = 0 #month after removed, set crack size to zero
  data = data %>% fill(crack_size, .direction = "down")  #LOCF crack size
  data$crack_size[is.na(data$crack_size)] = 0 #overwrite NaNs
  data = data %>% mutate(across(all_of(c("milling", "grinding", "planing")), ~replace_na(., 0))) %>% as.data.frame() #fill maintenance data as one-hot
  if (any(data$crack_size > 0)) return(data[,intersect(original_cols, names(data))]) else return(data[0, intersect(original_cols, names(data))])
}

## THIS EXTRACTION HAS BEEN PARRALLELIZED ON THE HPC CLUSTER.
# IDs = unique(d$ID)
# interpolated_data <- list()
# for (i in seq_along(IDs)) {
#   foo = interp_data(data = d, ID = IDs[i],constant_cols = constant_cols, measurement_cols = measurement_cols[!grepl("d0",measurement_cols)], order = 1, cutoff_date = as.Date("2011-12-01"))
#   interpolated_data[[i]] = foo
# }
# interpolated = bind_rows(interpolated_data)


### QUESTIONS ###
# - Correlation between cracks and meas (also lagged versions)
# - When do cracks initiate? Logistic regression.
# - How cracks propagate? diff crack_growth ~ time + meas

rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
d = fread('crack_data_interpolated_1order.csv', header = T, sep = ',')
d = as.data.frame(d)

measurement_cols = c("wear_h_mu", "wear_h_sigma", "wear_h_rms", "wear_h_p01","wear_h_p99", "wear_w_mu", "wear_w_sigma","wear_w_rms", "wear_w_p01", "wear_w_p99",
                     "corrugation_mu", "corrugation_sigma","corrugation_rms", "corrugation_p01", "corrugation_p99", "inclination_mu","inclination_sigma", "inclination_rms",
                     "inclination_p01", "inclination_p99","twist2m_mu", "twist2m_sigma","twist2m_rms", "twist2m_p01","twist2m_p99",  "twist3m_mu","twist3m_sigma","twist3m_rms",
                     "twist3m_p01","twist3m_p99","gauge_mu","gauge_sigma","gauge_rms","gauge_p01","gauge_p99","cant_mu","cant_sigma","cant_rms","cant_p01","cant_p99","align_horiz_d0_mu",
                     "align_horiz_d0_sigma","align_horiz_d0_rms","align_horiz_d0_p01","align_horiz_d0_p99","align_horiz_d1_mu","align_horiz_d1_sigma","align_horiz_d1_rms","align_horiz_d1_p01",
                     "align_horiz_d1_p99","align_vert_d0_mu","align_vert_d0_sigma","align_vert_d0_rms","align_vert_d0_p01","align_vert_d0_p99","align_vert_d1_mu","align_vert_d1_sigma","align_vert_d1_rms","align_vert_d1_p01","align_vert_d1_p99")



idx_class5 = d$ID %in% d$ID[which(d$crack_class == 5)]
idx_class4 = d$ID %in% d$ID[which(d$crack_class == 4)]
idx_class1 = d$ID %in% d$ID[which(d$crack_class == 1)]

#CORRELATION STRUCTURE
cormat_class5 = cor(as.matrix(d[idx_class5, colnames(d) %in% c("crack_size", measurement_cols)]), use = "complete.obs")
cormat_class4 = cor(as.matrix(d[idx_class4, colnames(d) %in% c("crack_size", measurement_cols)]), use = "complete.obs")
cormat_class1 = cor(as.matrix(d[idx_class1, colnames(d) %in% c("crack_size", measurement_cols)]), use = "complete.obs")
cormat_crack = cbind(cormat_class5[,1],cormat_class4[,1],cormat_class1[,1])[-1, ]
colnames(cormat_crack) = c("Head_checks", "Weld_cracks", "Squats")
round(cormat_crack,4)

# In general:
#   - align_horiz_* variables show consistent and strong positive correlations across all defect classes, 
#     particularly in their central tendency (mu), spread (sigma), and overall magnitude (rms). This suggests that
#     horizontal alignment irregularities are a robust indicator of crack presence and severity, regardless of class.

#   - twist2m_* and twist3m_* variables (especially sigma, rms, and p99) are positively correlated with crack size 
#     for all classes, though with slightly higher values for weld cracks and squats. These results support the 
#     understanding that twist (variation in cant over short base lengths) reflects local instability that may 
#     promote defect growth.

#   - cant_* variables (mu, p99, rms) show a **divergent pattern**:
#       - For head checks (class 5), they correlate *positively* with crack size.
#       - For squats (class 1), they correlate *negatively*.
#     This might reflect differing failure mechanisms: head checks tend to form under higher cant and wheel load 
#     concentration in curves, while squats may appear more under conditions of insufficient lateral support.

#   - gauge_* variables (mu, sigma, p99, etc.) have moderate positive correlations across classes, especially class 5.
#     This suggests that widened gauge may allow lateral wheel motion or stress redistribution, potentially contributing 
#     to surface cracks like head checks.

#   - corrugation_* and wear_* variables show mostly weak or negative correlations across all classes:
#       - Slightly more negative for squats (class 1), indicating they may form under smoother track but higher impact.
#       - Head checks show near-zero or weakly negative correlation, implying corrugation is not a strong driver.
#     These patterns suggest that while wear and surface roughness may degrade ride quality, they are not strongly
#     predictive of crack severity, especially in the presence of dynamic geometry variables (alignment, twist).

#   - inclination_* and alignment_vert_* variables are weakly correlated overall. They contribute marginally to 
#     the correlation structure and may act as secondary factors or covariates rather than direct drivers.

#In general, weaker correlation with squats than head checks.

#Testing linear model:
summary(lm(formula = crack_size ~ align_horiz_d1_mu + I(grinding + milling + planing) + cant_mu +twist2m_sigma, d[idx_class5, ]))
summary(lm(formula = crack_size ~ align_horiz_d1_mu + I(grinding + milling + planing) + cant_mu +twist2m_sigma, d[idx_class4, ]))
summary(lm(formula = crack_size ~ align_horiz_d1_mu + I(grinding + milling + planing) + cant_mu +twist2m_sigma, d[idx_class1, ]))
#estimates in accordance with physical interpretation of variable influence on cracks

lag.max = 36

lags_df_class1  <- d[idx_class1, ] %>% 
  group_by(ID) %>%
  reframe(
    max_lag_twist = {
      ts1 <- crack_size; ts2 <- twist2m_sigma
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$lag[which.max(abs(cc$acf))]
    },
    corr_twist = {
      ts1 <- crack_size; ts2 <- twist2m_sigma
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$acf[which.max(abs(cc$acf))]
    },
    max_lag_align = {
      ts1 <- crack_size; ts2 <- align_horiz_d1_mu
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$lag[which.max(abs(cc$acf))]
    },
    corr_align = {
      ts1 <- crack_size; ts2 <- align_horiz_d1_mu
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$acf[which.max(abs(cc$acf))]
    }
  ) %>% as.data.frame()
lags_df_class4  <- d[idx_class4, ] %>% 
  group_by(ID) %>%
  reframe(
    max_lag_twist = {
      ts1 <- crack_size; ts2 <- twist2m_sigma
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$lag[which.max(abs(cc$acf))]
    },
    corr_twist = {
      ts1 <- crack_size; ts2 <- twist2m_sigma
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$acf[which.max(abs(cc$acf))]
    },
    max_lag_align = {
      ts1 <- crack_size; ts2 <- align_horiz_d1_mu
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$lag[which.max(abs(cc$acf))]
    },
    corr_align = {
      ts1 <- crack_size; ts2 <- align_horiz_d1_mu
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$acf[which.max(abs(cc$acf))]
    }
  ) %>% as.data.frame()
lags_df_class5  <- d[idx_class5, ] %>% 
  group_by(ID) %>%
  reframe(
    max_lag_twist = {
      ts1 <- crack_size; ts2 <- twist2m_sigma
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$lag[which.max(abs(cc$acf))]
    },
    corr_twist = {
      ts1 <- crack_size; ts2 <- twist2m_sigma
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$acf[which.max(abs(cc$acf))]
    },
    max_lag_align = {
      ts1 <- crack_size; ts2 <- align_horiz_d1_mu
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$lag[which.max(abs(cc$acf))]
    },
    corr_align = {
      ts1 <- crack_size; ts2 <- align_horiz_d1_mu
      cc <- ccf(ts2, ts1, lag.max = lag.max, plot = FALSE)
      cc$acf[which.max(abs(cc$acf))]
    }
  ) %>% as.data.frame()


table(lags_df_class5$max_lag_twist)
table(lags_df_class4$max_lag_twist)
table(lags_df_class1$max_lag_twist)

table(lags_df_class5$max_lag_align)
table(lags_df_class4$max_lag_align)
table(lags_df_class1$max_lag_align)

# Class-wise interpretation of lag patterns:
#
# Class 5 (head checks):
#   - Peak at lag = 0, but long tail toward negative lags (e.g., -1 to -36).
#   - Interpretation: Track geometry (twist/alignment) often changes *before* crack size increases.
#     => Indicates a *predictive* relationship: geometric irregularities may lead to head check formation.
#
# Class 4 (weld cracks):
#   - Dominant peak at lag = 0, with a fairly symmetric distribution.
#   - Interpretation: Crack size and geometry tend to change simultaneously.
#     => Suggests a more *static* or *localized* cause (e.g., welding defects), not dynamically evolving geometry.
#
# Class 1 (squats):
#   - Peak at lag = 0, but strong skew toward *positive lags* (1 to 15+).
#   - Interpretation: Crack size often increases *before* geometric deviations occur.

##############################
### CRACK INITIATION MODEL ###
##############################
d_initiation <- d[idx_class5 & d$ID > 50000, ] %>% arrange(ID, Date_month) %>% group_by(ID) %>% filter(!first(crack_size > 0)) %>% mutate(crack_initiation = crack_size > 0,post_crack = cumsum(crack_initiation) > 1) %>%filter(!post_crack) %>%mutate(crack_initiation = as.integer(crack_initiation)) %>% ungroup()
d_initiation <- d_initiation %>% select(crack_initiation, State, everything()) %>% as.data.frame()

summary(glm(crack_initiation ~ align_horiz_d1_mu + twist2m_sigma, data = d_initiation, family = "binomial")) #fixed logistic regression

d_initiation$ID2 <- match(d_initiation$ID, unique(d_initiation$ID)) #relabel IDs



library(TMB)
nam = "logistic"
compile(paste0(nam, ".cpp"))
dyn.load(dynlib(nam))  

X <- model.matrix(~ align_horiz_d1_mu + twist2m_sigma, data = d_initiation)  # add covariates here
Y <- d_initiation$crack_initiation  # binary response

data <- list(Y = Y, X = X)
parameters <- list(beta = runif(NCOL(X),-1,1))

obj <- MakeADFun(data, parameters, DLL = nam,hessian = T)
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt #model estimates
sqrt(diag(solve(obj$he()))) #standard errors


nam = "logistic_mixed"
compile(paste0(nam, ".cpp"))
dyn.load(dynlib(nam))  

X <- model.matrix(~ align_horiz_d1_mu + twist2m_sigma, data = d_initiation)  # add covariates here
Y <- d_initiation$crack_initiation  # binary response
id <- as.integer(factor(d_initiation$ID2)) - 1  # zero-based factor

data <- list(Y = Y, X = X, ID = id)
parameters <- list(beta = runif(NCOL(X),-1,1), u = rep(0.1, length(unique(id))),log_sigma_u = -2)

obj <- MakeADFun(data, parameters, random = "u", DLL = nam,hessian = T)
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt #model estimates
sqrt(diag(solve(obj$he()))) #standard errors

#Trouble shoot: 
# - find glmm examle from TMB package and try on this windows: https://github.com/admb-project/tmb-examples/blob/master/orange-trees/ora.cpp
# - try this code on the HPC cluster

library(glmmTMB)
glmmTMB(crack_initiation ~ 1 + (1 | ID2), data = d_initiation, family = binomial)








#create target for logistic regression.(need to filter out IDs that dont have zero crack at first obs***)

d <- d %>% group_by(ID) %>% mutate(crack_initiated = crack_size > 0,crack_start = crack_initiated & !lag(crack_initiated, default = FALSE))



p1 = d %>% ggplot(aes(x = Date_month, y = crack_size, group = ID)) + geom_line(alpha = 0.2)



fit = lm(crack_size ~ wear_h_mu + wear_w_mu + corrugation_mu + inclination_sigma + twist2m_mu + twist3m_sigma + gauge_mu + cant_mu, data = d)
summary(fit)




### PRESENTS RESULTS FOR ###
# - descriptive data analysis
# - correlation structures
# - crack initiation
## - logistic regression
## - logistic regression with random effects
## - assessment of predictive cabability (locf and linear interp. data comparison as well)

# - Do above for different crack types (class 5, 4, and 1).











