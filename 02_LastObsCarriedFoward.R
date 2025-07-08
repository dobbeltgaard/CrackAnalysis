
rm(list = ls())

### LOCF per defect level
### Objective: To analyse crack growth in terms of track behavior.
### Assumption: assuming last observation is valid until new measurement, creating regular temporal grid (1month resolution)

d = read.csv("crack_data_v2.csv")
#d_copy = d

library(dplyr)
library(tidyr)
library(lubridate)
library(data.table)
library(zoo)


d$Date <- as.Date(d$Date)
d <- d %>% mutate(Date_month = floor_date(Date, "month"))

# colnames(foo)
# [1] "Date_month"              "ID"                      "Date"                    "dTime"                  
# [5] "crack_size"              "crack_class"             "crack_class_description" "in_straight"            
# [9] "in_transition_curve"     "in_curve"                "covariates_complete"     "rail_age"               
# [13] "R200"                    "R260"                    "R350HT"                  "grinding"               
# [17] "milling"                 "planing"                 "rail_1m_weight"          "traffic_direction"      
# [21] "Curve"                   "Overheight"              "MGT_min"                 "MGT_max"                
# [25] "EMGT_min"                "EMGT_max"                "Line_speed"              "Turnout_indicator"      
# [29] "wear_h_mu"               "wear_h_sigma"            "wear_h_rms"              "wear_h_p01"             
# [33] "wear_h_p99"              "wear_w_mu"               "wear_w_sigma"            "wear_w_rms"             
# [37] "wear_w_p01"              "wear_w_p99"              "corrugation_mu"          "corrugation_sigma"      
# [41] "corrugation_rms"         "corrugation_p01"         "corrugation_p99"         "inclination_mu"         
# [45] "inclination_sigma"       "inclination_rms"         "inclination_p01"         "inclination_p99"        
# [49] "twist2m_mu"              "twist2m_sigma"           "twist2m_rms"             "twist2m_p01"            
# [53] "twist2m_p99"             "twist3m_mu"              "twist3m_sigma"           "twist3m_rms"            
# [57] "twist3m_p01"             "twist3m_p99"             "gauge_mu"                "gauge_sigma"            
# [61] "gauge_rms"               "gauge_p01"               "gauge_p99"               "cant_mu"                
# [65] "cant_sigma"              "cant_rms"                "cant_p01"                "cant_p99"               
# [69] "align_horiz_d0_mu"       "align_horiz_d0_sigma"    "align_horiz_d0_rms"      "align_horiz_d0_p01"     
# [73] "align_horiz_d0_p99"      "align_horiz_d1_mu"       "align_horiz_d1_sigma"    "align_horiz_d1_rms"     
# [77] "align_horiz_d1_p01"      "align_horiz_d1_p99"      "align_vert_d0_mu"        "align_vert_d0_sigma"    
# [81] "align_vert_d0_rms"       "align_vert_d0_p01"       "align_vert_d0_p99"       "align_vert_d1_mu"       
# [85] "align_vert_d1_sigma"     "align_vert_d1_rms"       "align_vert_d1_p01"       "align_vert_d1_p99"      
# [89] "passages_maintenance"    "stone_maintenance"       "removed_material"        "maintenance_source"     
# [93] "defect_source"           "track_source"            "rail_left"               "visible"                
# [97] "combined_defect"         "under_limit"             "BTR"                     "Track"                  
# [101] "From"                    "To"                      "Defect_group"            "Profile"                
# [105] "Fastening"               "State"                   "detection_method"        "cluster_method"         
# [109] "crack_counts"           
# 


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
  #first_non_na <- function(x) {out <- na.omit(x); if (length(out) == 0) NA else out[1]} #function to summarize all but measurement_cols
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

IDs = unique(d$ID)
interpolated_data <- list()
for (i in seq_along(IDs)) {
  foo = interp_data(data = d, ID = IDs[i],constant_cols = constant_cols, measurement_cols = measurement_cols[!grepl("d0",measurement_cols)], order = 1, cutoff_date = as.Date("2011-12-01"))
  interpolated_data[[i]] = foo
}



interpolated = bind_rows(interpolated_data)







write.csv(interpolated_data, file = "crack_data_interpolated_list.csv", row.names = FALSE)





as.data.frame(interpolated_data[[1]])






foo1 = interp_data(data = d, ID = unique(d$ID)[2001],constant_cols = constant_cols, measurement_cols = measurement_cols[!grepl("d0",measurement_cols)], order = 1, cutoff_date = as.Date("2013-12-01"))
View(foo1)


cor_per_id <- d %>%
  group_by(ID) %>%
  summarise(cor_wear_incl = cor(wear_h_mu, inclination_mu, use = "pairwise.complete.obs"),
            cor_crack_gauge = cor(crack_size, gauge_mu, use = "pairwise.complete.obs"),
            .groups = "drop")

ccf_result <- ccf(interp_data$wear_interp, interp_data$crack_size_interp, na.action = na.pass)

### QUESTIONS ###
# - Collinearity in measurement variables
# - Correlation between cracks and meas (also lagged versions)
# - When do cracks initiate? Logistic regression.
# - How cracks propagate? diff crack_growth ~ time + meas


# create_locf_data <- function(data,ID,constant_cols, measurement_cols) {
#   data = data[data$ID == ID,]
#   
#   data <- data %>% mutate(Date_month = floor_date(Date, "month")) %>% arrange(Date_month) %>% 
#     complete(Date_month = seq.Date(min(Date_month), max(Date_month), by = "month")) %>%
#     fill(all_of(constant_cols), .direction = "updown") %>% fill(all_of(measurement_cols), .direction = "down") %>% 
#     filter(if_all(all_of(measurement_cols), ~ !is.na(.)))
#   
#   if (any(!is.na(data$crack_size))) return(data) else return(data[0, ])
# }
# 
# 
# create_linear_interp_data <- function(data, ID, constant_cols, measurement_cols) {
#   data <- data[data$ID == ID, ]
#   
#   data <- data %>% mutate(Date_month = floor_date(Date, "month")) %>% arrange(Date_month) %>%complete(Date_month = seq.Date(min(Date_month), max(Date_month), by = "month"))
#   data <- data %>% fill(all_of(constant_cols), .direction = "updown") # LOCF/BOCF for constant columns
#   for (col in measurement_cols) { data[[col]] <- na.approx(data[[col]], x = data$Date_month, na.rm = FALSE) } # Linear interpolation for measurement columns
#   data <- data %>% filter(if_all(all_of(measurement_cols), ~ !is.na(.)))
#   
#   if (any(!is.na(data$crack_size))) return(data) else return(data[0, ])
# }
