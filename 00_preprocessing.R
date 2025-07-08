
####################
### EDDY CURRENT ###
####################
rm(list = ls())
library(readxl)
library(dplyr)
library(igraph)
source("C:/Users/askbi/OneDrive - Danmarks Tekniske Universitet/SpatioTemporal/99_translations_and_naming_convetions.R")
path = "C:/Users/askbi/OneDrive - COWI/Documents - A235142 - Predictive maintanance ph.d project (1502)/Data/BDK/"
D = NULL

E = read_excel(list.files(path, full.names = T)[grepl("Eddy",list.files(path))])
E$Skinnestreng = translate_strings(E$Skinnestreng)
E = as.data.frame(E)[, c("BTR", "Spor", "Event (Date)", "From", "To", "Skinnestreng", "Dybde")]
names(E) <- c("BTR", "Track", "Date", "From", "To", "Rail_string", "Depth")
E$Date <- as.Date(E$Date)
E <- E[!duplicated(E[, c("BTR", "Track", "Date", "From", "To", "Rail_string", "Depth")]), ]
E$Year =  as.numeric(format(as.Date(E$Date), "%Y"))
E = E[order(E$Date, decreasing = T),] 
E$IDfoo = seq_len(NROW(E))
E$ID = NA

average_overlap <- function(a_from, a_to, b_from, b_to) {
  len_a <- a_to - a_from
  len_b <- b_to - b_from
  overlap <- max(0, min(a_to, b_to) - max(a_from, b_from))  # Calculate overlap length
  if (len_a <= 0 || len_b <= 0) return(0)
  # Compute average of overlap ratios
  overlap_a <- overlap / len_a
  overlap_b <- overlap / len_b
  avg_overlap <- (overlap_a + overlap_b) / 2
  return(avg_overlap)
}

start.time <- Sys.time()
used_elements = c()
overlap_parameter = 0.5
for(i in 1:NROW(E)){ 
  idx1 = !(E$IDfoo %in% used_elements) & E[i, c("BTR")] == E[, c("BTR")] & E[i, c("Track")] == E[, c("Track")] & E[i, c("Rail_string")] == E[, c("Rail_string")] 
  if(NROW(E[idx1, ]) < 1) next
  foo = E[idx1, ]#[idx2, ]
  foo = foo[order(foo$Date, decreasing = T),] 
  
  #calculate similarities
  similarity_matrix = matrix(NA,  NROW(foo),  NROW(foo))
  for(j in 1:NROW(foo)){
    for(k in 1:NROW(foo)){
      similarity_matrix[j,k]=average_overlap(a_from = foo$From[j], a_to = foo$To[j], b_from = foo$From[k], b_to = foo$To[k]) * (abs(foo$Date[j] - foo$Date[k]) > 60)
    }
  }
  bin_mat = similarity_matrix > overlap_parameter; diag(bin_mat) <- FALSE 
  
  #clustering algo. 
  match_check = T
  while(match_check){
    g = graph_from_adjacency_matrix(bin_mat, mode = "undirected") #make graph
    cliques = max_cliques(g) #find cliques
    
    #calculate clique scores
    rewards = rep(0, length(cliques)) 
    for(j in 1:length(rewards)){
      indices = cliques[[j]]
      res = 0
      if (length(indices) >= 2) {
        mat = similarity_matrix[indices, indices]
        res = sum(mat[lower.tri(mat)])
      }
      rewards[j] = res
    }
    if(any(rewards > overlap_parameter)) match_check = TRUE else { match_check = FALSE; break }
    match_idx = cliques[[which.max(rewards)]]
    if (length(match_idx) >= nrow(bin_mat)) break
    matched = foo[match_idx,]; matched = matched[order(matched$Date, decreasing = F),] #find best clique
    matched$ID = min(matched$IDfoo)
    if(is.null(D)) { D = matched} else {D = rbind(D, matched)}
    used_elements = c(used_elements, matched$IDfoo)
    
    bin_mat = bin_mat[-match_idx,-match_idx] #remove matched elements
    similarity_matrix = similarity_matrix[-match_idx,-match_idx] #remove matched elements
    foo = foo[-match_idx, ] #update candidates
  }
  end.time <- Sys.time()
  if(i %% 1000 == 0) print(end.time - start.time)
}
Esingles = E[!E$IDfoo %in% unique(D$IDfoo),] # find singles in E
Esingles$ID = Esingles$IDfoo #Assign singles their unique ID
DD = rbind(D, Esingles) #bind singles and matches
#write.csv(x = DD[, c("ID","BTR", "Track", "Rail_string", "From", "To", "Date", "Depth")], row.names = F,  file = "defects_eddy_current.csv")

##################
### UT DEFECTS ###
##################
rm(list = ls())
library(readxl)
library(stringr)
library(plyr); library(dplyr)
library(tibble)
source("C:/Users/askbi/OneDrive - Danmarks Tekniske Universitet/SpatioTemporal/99_translations_and_naming_convetions.R")
path = "C:/Users/askbi/OneDrive - COWI/Documents - A235142 - Predictive maintanance ph.d project (1502)/Data/BDK/"

UT_files = list.files(paste0(path, "BDK_UT"), full.names = T)
UT_files = UT_files[grepl("Defects", UT_files)]
mat = matrix(NA, nrow = length(UT_files), ncol = 60)
for(i in 1:length(UT_files)){
  f = read_excel(UT_files[i])
  mat[i, 1:length(names(f))] = names(f) 
  if(i == 1){ ff = f} else { ff = rbind(ff,f) }
}

cols.remove = 
  c("Event (Time)",  "UT-operatør", "GPS", "Svejsernr.","Svejsning lbnr.","Journal nr.",
    "B-scan nummers","Flere over 3 sveller","Hovedbredde forøget","Fejlbeskrivelse",
    "Fejlens status","Målt håndholdt","Lasket fejl","Valseværk","Måned","Højdeslid",
    "Sideslid","Indpasser","Svejsningens nummer","Fejlen kan fjernes ved","Sikret med lasker",
    "Sikret ved LA","Defect class","Isolerklæbestød","Fejl udbedret","FejlenKanSikresMed",
    "clDate","closedBy", "Gammel fejlIDs")
ff = ff[,!names(ff) %in% cols.remove] #remove irrelevant columns
ff = ff[!is.na(ff$FejlID), ] #remove NaNs
ff = ff[!grepl("^5000", as.character(ff$FejlID)), ] #remove defects starting with 5000
ff = ff[!ff$`Under grænsen` == "Yes", ]
ff = ff[!ff$`Fejl fundet i` == "Ingen fejl", ]

#TRANSLATION
names(ff) = translate_strings(names(ff))
fff = ff
fff[, c("Rail_string","Defect_found_in","Track_type","Visible","Combined_defect","Under_limit","Curvature","State")] = 
  apply(ff[,c("Rail_string","Defect_found_in","Track_type","Visible","Combined_defect","Under_limit","Curvature","State")], 2, translate_strings)
fff = fff[!apply(X = fff[ ,c("BTR", "Track", "From", "To", "Rail_string")],1, function(x) any(is.na(x))), ] #remove NaNs in important columns
d <- as.data.frame(fff)
d$Method = "UT"
d[d$Defect_ID == -1,"Defect_ID"] = 1012 #overwrite only negative ID 
d$Date <- as.Date(d$Date)

#Matching within UTs
eps = 0.001
d$Defect_ID2 = NA
for(i in 1:NROW(d)){
  #check = (d$BTR == d$BTR[i] & d$Track == d$Track[i] & d$From == d$From[i] &d$To == d$To[i] & d$Rail_string == d$Rail_string[i] & d$Defect_ID != d$Defect_ID[i])
  check = (d$BTR == d$BTR[i] & d$Track == d$Track[i] & abs(d$From - d$From[i]) < eps & abs(d$To - d$To[i]) < eps & d$Rail_string == d$Rail_string[i] & d$Defect_ID != d$Defect_ID[i])
  if(sum(check) > 0){ d$Defect_ID2[i] = paste0(d[c(i,which(check)), "Defect_ID"],collapse=" ") }
}
d$Defect_ID3 = NA
for(i in 1:NROW(d)){ if(!is.na(d$Defect_ID2[i])) d$Defect_ID3[i] = min(as.numeric(unlist(strsplit(d$Defect_ID2[i], " ")))) } #use minimum defect number for match cycles
d$Defect_ID3[is.na(d$Defect_ID3)] = d$Defect_ID[is.na(d$Defect_ID3)]

IDs = unique(d$Defect_ID3) #collect defect IDs
IDs = IDs[!is.na(IDs)] #remove nans
defect_set = NULL
for(i in 1:length(IDs)){
  idx = d$Defect_ID3 == IDs[i]
  foo = d[idx,]
  x = split(foo,list(foo$Rail_string,foo$Defect_found_in), drop=TRUE) #split to unique defects, so rail string and defect_found_in surely matches
  
  for(j in 1:NROW(x)){
    if(NROW(x) < 2){ nam = toString(IDs[i])} else {nam = paste0(sprintf("%.0f", 100000000),toString(IDs[i]),j) }# defect_set[[toString(IDs[i])]] = x[[1]][order(x[[1]]$Date),]; next} #x[[1]]; next} #if data correct, then no splitting is needed
    foo2 = x[[j]][order(x[[j]]$Date),]
    
    if(NROW(foo2) > 1){
      idxfoo2 = rep(T, NROW(foo2));  
      if(any(head(foo2$State == "Removed",-1)) & !all(foo2$State == "Removed")){ idxfoo2[foo2$State == "Removed"] = F; idxfoo2[length(idxfoo2)]=T;} #filter out observations with "Removed" under state, excect from endpoint
      foo2 = foo2[idxfoo2, ]
      if(sum(idxfoo2) > 1){ #if more than 1 obs is left, then check for time differences
        idx2foo2 = rep(T, NROW(foo2));
        idx2foo2 = outer(foo2$Date, foo2$Date, FUN = difftime, units = "days")[,1] > 30 #only obs with bigger distance than 30 days
        idx2foo2[!idx2foo2][which.max(foo2$Date[!idx2foo2])] = T
        foo2 = foo2[idx2foo2, ]
      }
    }
    foo2$Maintenance_date = NA
    foo2$Maintenance_type = NA
    foo2$Amount = NA
    foo2$Passages = NA
    foo2$ID = as.numeric(nam)
    if(is.null(defect_set)) { defect_set = foo2} else {defect_set = rbind(defect_set, foo2)}
    }
  if(i %% 1000 == 0){print(i)}
  }
#write.csv(x = defect_set[,!colnames(defect_set) %in% c("Defect_ID", "Defect_ID2", "Defect_ID3")], row.names = F,  file = "ultrasonic_testing_defects.csv")
#save(defect_set, file = "ultasonic_testing_defects.RData")

#############################################
### UT DEFECTS COMBINED WITH EDDY CURRENT ###
#############################################
rm(list = ls())
library(dplyr)
U = read.csv("defects_ultrasonic_testing.csv")
E = read.csv("defects_eddy_current.csv")

#format EC like UT
U$Date <- as.Date(U$Date)
E$Date <- as.Date(E$Date)
E$`Switch_nr.` = NA; E$Defect_length = NA; E$Defect_width = NA; E$Depth_from = 0; E$Depth_to = E$Depth; 
E$Defect_group = NA; E$Defect_found_in = NA; E$Track_type = NA; E$Visible = NA; 
E$Combined_defect = NA; E$Under_limit = NA; E$Profile = NA; E$Curvature=NA; E$Year = NA; E$Steel = NA; 
E$Fastening = NA; E$Switch_type = NA; E$State = "Open"; E$Defect_group_calc = NA; E$UIC = NA; E$Method = "EC"
E$Maintenance_date = NA; E$Maintenance_type = NA; E$Amount = NA; E$Passages = NA; 
cols = c("BTR","Track","Date","From","To","Rail_string", "Switch_nr.",
         "Defect_length","Defect_width","Depth_from","Depth_to","Defect_group","Defect_found_in",
         "Track_type","Visible","Combined_defect","Under_limit","Profile","Curvature","Year","Steel",
         "Fastening","Switch_type" ,"State","Defect_group_calc","UIC","Method", 
         "Maintenance_date", "Maintenance_type", "Amount", "Passages", "ID")
E = E[,cols]
E$ID = -(E$ID) #reverse sign of EC IDs
E$ID_UT = NA 


#Requirements for match: 
# latest EC is earlier than earliest UT. But time difference must be 60 days.
# Earliest UT is contained by all EC in a clustrer. No slack allowed, due relative large span of EC.
# No splitting of UT or EC, into entire new Defect clusters.

#Process is: 
# Loop through UT IDs. 
# Choose the earliest.
# Check for matches with EC clusters, based on the requirements for a match. 
# If one match, assign UT ID to ECs in the matched cluster.
# If multiple matches, choose EC cluster with latest obs. closest in time to the earliest UT, still not exceeding 60 days limit.
# If no matches, proceed to next UT ID.


start.time <- Sys.time()
iters = 0
for(ut_id in unique(U$ID) ){
  iters = iters + 1
  ut_group <- U[U$ID == ut_id, ]
  ut_date_min <- min(ut_group$Date) # Get earliest UT observation
  ut_from <- min(ut_group$From)
  ut_to   <- max(ut_group$To)
  
  idx = ut_group$BTR[1] == E$BTR & ut_group$Track[1] == E$Track & ut_group$Rail_string[1] == E$Rail_string 
  ec_candidates <- list()
  for (ec_id in unique(E$ID[idx])) { # Loop through EC clusters at same BTR, Track and rail string as UT group
    ec_group <- E[E$ID == ec_id, ]
    ec_date_max <- max(ec_group$Date)
    if (as.numeric(ut_date_min - ec_date_max) < 30) next #latest EC must be at least 30 days before earliest UT
    
    ec_from_max <- max(ec_group$From)
    ec_to_min   <- min(ec_group$To)
    ec_ok <- ec_from_max <= ut_from & ec_to_min >= ut_to
    #ec_ok <- all(ec_group$From <= ut_from & ec_group$To >= ut_to) #all ECs in cluster must fully contain UT
    if (!ec_ok) next
    ec_candidates[[as.character(ec_id)]] <- ec_group # Candidate passes all checks
  }
  
  if (length(ec_candidates) == 1) { # One match: assign it
    ec_match_id <- as.numeric(names(ec_candidates)[1])
    E$ID_UT[E$ID == ec_match_id] <- ut_id
  } else if (length(ec_candidates) > 1) { # Multiple matches: pick the one with EC date closest to UT, still outside day span
    candidate_diffs <- sapply(ec_candidates, function(g) abs(ut_date_min - max(g$Date)))
    ec_match_id <- as.numeric(names(ec_candidates)[which.min(candidate_diffs)])
    E$ID_UT[E$ID == ec_match_id] <- ut_id
  }
  if(iters %% 1000 == 0){end.time <- Sys.time();print(end.time - start.time)}
}
U$ID_UT = NA
D = rbind(E, U)
D$ID_UT[is.na(D$ID_UT)] = D$ID[is.na(D$ID_UT)]

D_sorted <- D %>%
  group_by(ID_UT) %>%
  arrange(Date, .by_group = TRUE) %>%
  ungroup()
D_sorted = as.data.frame(D_sorted)
D_sorted$ID = D_sorted$ID_UT
#write.csv(x = D_sorted[,!colnames(D_sorted) %in% c("ID_UT")], row.names = F,  file = "defects_UC_EC.csv")

#####################
### CLEANING DATA ###
#####################
rm(list = ls())
library(dplyr)
d = read.csv("defects_UC_EC.csv")

#track info, cleaning, etc.
#What about Wear?

#Description of cluster method
foo <- d %>%
  group_by(ID) %>%
  filter(n_distinct(Method) > 1) %>%
  arrange(Date) %>%
  slice(1) %>%
  ungroup()
foo = as.data.frame(foo)
idx = d$ID %in% foo$ID
d$Cluster_method = d$Method
d$Cluster_method[idx] = "EC+UT"

#Summary of sequence size per cluster method
mat = matrix(NA, 12,4)
mat[1:length(table(table(d$ID[d$Cluster_method == "EC"]))), 1] = table(table(d$ID[d$Cluster_method == "EC"]))
mat[1:length(table(table(d$ID[d$Cluster_method == "UT"]))), 2] = table(table(d$ID[d$Cluster_method == "UT"]))
mat[2:(1+length(table(table(d$ID[d$Cluster_method == "EC+UT"])))), 3] = table(table(d$ID[d$Cluster_method == "EC+UT"]))
mat[1:length(table(table(d$ID))), 4] = rowSums(mat, na.rm = T)
rownames(mat) = 1:12
colnames(mat) = c("EC", "UT", "EC+UT", "Total")

#################################################
### COMBINE DEFECT DATA WITH MAINTENANCE DATA ###
#################################################
main = read.csv(file = file.path("rail_maintenance.csv"))
iters = 0; start.time <- Sys.time(); 
for(i in unique(d$ID) ){
  candidates = d[d$ID == i,]
  if(NROW(candidates)>1){ #if there are more than one observation
    idx = candidates$BTR[1] == main$BTR & candidates$Track[1] == main$Track
    main_sub <- main[idx, ]
    overlaps <- pmin(max(candidates$To), main_sub$To) > pmax(min(candidates$From), main_sub$From)
    foo <- main_sub[overlaps, ]
    if(NROW(foo)>0){
      for(j in 1:NROW(foo)){
        cond3 = foo$Time[j] > candidates$Date  
        cond4 = foo$Time[j] < candidates$Date
        if(any(cond3) & any(cond4)){ #if maintenance is between
          candidates[which(cond4)[1],c("Maintenance_date","Maintenance_type","Amount", "Passages")] =
            foo[j,c("Time", "action","removed_mm", "Passages") ]
        }
      }
      d[match(rownames(candidates), rownames(d)), ] <- candidates
    }
  }
  iters = iters + 1
  if(iters %% 1000 == 0){end.time <- Sys.time();print(end.time - start.time); start.time <- Sys.time()}
}
#write.csv(x = d, row.names = F,  file = "defects_UC_EC_maintenance.csv")

#############################################################
#### COMBINE with static data. Tonnage. Speed. Curveture. ###
#############################################################
rm(list = ls())
library(readxl)
source("C:/Users/askbi/OneDrive - Danmarks Tekniske Universitet/SpatioTemporal/99_translations_and_naming_convetions.R") #source("99_translations_and_naming_convetions.R")
d = read.csv(file = "defects_UC_EC_maintenance.csv")

path = "C:/Users/askbi/OneDrive - COWI/Documents - A235142 - Predictive maintanance ph.d project (1502)/Data/BDK/"
asset.files = list.files(paste0(path, "Asset Data"), full.names = T)

#Read superstructure
super.file = asset.files[which(grepl("Superstructure",asset.files))]
R <- read_excel(super.file)
R$BTRn = as.numeric(R$BTR)
R$Sportype = translate_strings(R$Sportype) 
R$`Skinne-streng` = translate_strings(R$`Skinne-streng`) 
R$`Skinne: Type` = translate_strings(R$`Skinne: Type`)

#Read curvature
curve.file = asset.files[which(grepl("urvature",asset.files))]
C <- read_excel(curve.file)
C$DK_type_element = translate_strings(C$DK_type_element) 
curvature <- gsub(".*\\] ", "", C$Curvature)
curvature <- gsub("[{}]", "", curvature) #as.data.frame(do.call(rbind, str_match_all(curvature, "(-?\\d+\\.?\\d*)")))
curvature <- matrix(as.numeric(unlist(strsplit(curvature, ";"))), ncol=2,byrow = 1)
curvature <- apply(curvature, 1, mean)
C$curve = curvature
overheight <- gsub(".*\\] ", "", C$`Overhøjde [mm]`)
overheight <- gsub("[{}]", "", overheight) #as.data.frame(do.call(rbind, str_match_all(curvature, "(-?\\d+\\.?\\d*)")))
overheight <- matrix(as.numeric(unlist(strsplit(overheight, ";"))), ncol=2,byrow = 1)
overheight <- apply(overheight, 1, mean)
C$overheight = overheight

#Read loading
load.file = asset.files[which(grepl("load",asset.files))]
L <- read_excel(load.file)

#Read Speeds
speed.file = asset.files[which(grepl("speed",asset.files))]
S <- read_excel(speed.file)

#Read turnover segments
turnover.file = asset.files[which(grepl("segments_turnouts",asset.files))]
SC <- read_excel(turnover.file)
split_vec <- strsplit(gsub(" ", "", SC$`BTR-Spor`), "-")
left_part <- sapply(split_vec, `[`, 1)
right_part <- sapply(split_vec, `[`, 2)
SC$BTR = left_part
SC$Spor = right_part
SC$BTRn = as.numeric(SC$BTR)

### GET ASSET INFORMATION ###
d$Rail_type2=NA;d$Track_type2=NA;d$Steel_type2=NA
for(i in unique(d$ID)) {
  idx_base <- d$ID == i
  candidates <- d[idx_base, ]
  idx <- candidates$BTR[1] == R$BTR & candidates$Track[1] == R$Spor & candidates$Rail_string[1] == R$`Skinne-streng`
  sub <- R[idx, ]
  defect_from <- min(candidates$From)
  defect_to   <- max(candidates$To)
  overlaps <- pmin(defect_to, sub$Til) > pmax(defect_from, sub$Fra)
  foo <- sub[overlaps, ]
  if (NROW(foo) > 1) {
    p1 <- c(defect_from, defect_to)
    p2s <- as.matrix(foo[, c("Fra", "Til")])
    idx2 <- which.min(colSums((t(p2s) - p1)^2))
    foo <- foo[idx2, ]
  } else if (NROW(foo) == 0) {
    if (NROW(sub) == 0) {
      d$Rail_type2[idx_base]  <- NA
      d$Track_type2[idx_base] <- NA
      d$Steel_type2[idx_base] <- NA
      next
    }
    p1 <- c(defect_from, defect_to)
    p2s <- as.matrix(sub[, c("Fra", "Til")])
    idx2 <- which.min(colSums((t(p2s) - p1)^2))
    foo <- sub[idx2, ]
  }
  d$Rail_type2[idx_base]  <- foo$`Skinne: Type`
  d$Track_type2[idx_base] <- foo$Sportype
  d$Steel_type2[idx_base] <- foo$`Skinne: Stålkvalitet`
}

### GET CURVATURE INFORMAION ###
d$Track_type3=NA;d$Curve=NA;d$Overheight=NA
for(i in unique(d$ID) ){
  idx_base <- d$ID == i
  candidates <- d[idx_base, ]
  idx <- candidates$BTR[1] == C$BTR & candidates$Track[1] == C$Spor
  
  sub <- C[idx, ]
  defect_from <- min(candidates$From)
  defect_to   <- max(candidates$To)
  overlaps <- pmin(defect_to, sub$Til) > pmax(defect_from, sub$Fra)
  foo <- sub[overlaps, ]

  if(NROW(foo) > 1){ #if there are more than one match, then find the minimum distance
    p1 <- c(defect_from, defect_to)
    p2s <- as.matrix(foo[, c("Fra", "Til")])
    idx2 = which.min(apply(p2s, MARGIN = 1, FUN = function(x) sum((x - p1)^2) )) #find minimum distance
    foo <- foo[idx2, ]
  } else if(NROW(foo)  == 0){ #if there are no matches, then it is straight track
    d$Track_type3[idx_base]  = "Straight_track"
    d$Curve[idx_base] = 0
    d$Overheight[idx_base] = 0
    next
  }
  d$Track_type3[idx_base] = foo$DK_type_element
  d$Curve[idx_base] = foo$curve
  d$Overheight[idx_base] = foo$overheight
}

### GET TONNAGE INFORMATION ###
d$MGT_min=NA;d$MGT_max=NA;d$EMGT_min=NA; d$EMGT_max = NA; 
for(i in unique(d$ID) ){
  idx_base <- d$ID == i
  candidates <- d[idx_base, ]
  idx <- candidates$BTR[1] == L$BTR & candidates$Track[1] == L$Spor
  
  sub <- L[idx, ]
  defect_from <- min(candidates$From)
  defect_to   <- max(candidates$To)
  overlaps <- pmin(defect_to, sub$To) > pmax(defect_from, sub$From)
  foo <- sub[overlaps, ]
  
  if(NROW(foo) > 1){ #if there are more than one match, then find the minimum distance
    p1 <- c(defect_from, defect_to)
    p2s <- as.matrix(foo[, c("From", "To")])
    idx2 = which.min(apply(p2s, MARGIN = 1, FUN = function(x) sum((x - p1)^2) )) #find minimum distance
    foo <- foo[idx2, ]
  } else if(NROW(foo) == 0){ #if there are no matches, then it is straight track
    d[idx_base, c("MGT_min", "MGT_max", "EMGT_min", "EMGT_max")] <- NA
    next
  }
  d$MGT_min[idx_base] = foo$MGT_min
  d$MGT_max[idx_base] = foo$MGT_max
  d$EMGT_min[idx_base] = foo$EMGT_min
  d$EMGT_max[idx_base] = foo$EMGT_max
}

### GET LINE SPEED ###
d$Line_speed=NA;
for(i in unique(d$ID) ){
  idx_base <- d$ID == i
  candidates <- d[idx_base, ]
  idx <- candidates$BTR[1] == S$BTR & candidates$Track[1] == S$Spor
  
  sub <- S[idx, ]
  defect_from <- min(candidates$From)
  defect_to   <- max(candidates$To)
  overlaps <- pmin(defect_to, sub$Til) > pmax(defect_from, sub$Fra)
  foo <- sub[overlaps, ]
  
  if(NROW(foo) > 1){ #if there are more than one match, then find the minimum distance
    p1 <- c(defect_from, defect_to)
    p2s <- as.matrix(foo[, c("Fra", "Til")])
    idx2 = which.min(apply(p2s, MARGIN = 1, FUN = function(x) sum((x - p1)^2) )) #find minimum distance
    foo <- foo[idx2, ]
  } else if(NROW(foo)  == 0){ #if there are no matches, then it is straight track
    d$Line_speed[idx_base]  = NA
    next
  }
  d$Line_speed[idx_base] = foo$`Line speed [km/h]`
}

### GET TURNOUT INDICATOR ###
d$Turnout_indicator=NA
for(i in unique(d$ID) ){
  idx_base <- d$ID == i
  candidates <- d[idx_base, ]
  idx <- candidates$BTR[1] == SC$BTR & candidates$Track[1] == SC$Spor
  sub <- SC[idx, ]
  defect_from <- min(candidates$From)
  defect_to   <- max(candidates$To)
  overlaps <- pmin(defect_to, sub$To) > pmax(defect_from, sub$From)
  foo <- sub[overlaps, ]
  
  if(NROW(foo) >= 1){ 
    d$Turnout_indicator[idx_base] = 1
  } else if(NROW(foo) == 0){ #if there are no matches, then it is straight track
    d[idx_base, c("Turnout_indicator")] = 0
  }
}
#write.csv(x = d, row.names = F,  file = "defects_UC_EC_maintenance_covariates.csv")

################
### CLEANING ###
################
rm(list = ls())
library(stringr)
d = read.csv(file = "defects_UC_EC_maintenance_covariates.csv")


#how many with same ID, dont share BRT, Track, Rail_string
save = c()
for(i in 1:length(unique(d$ID))){
  if(any(1!=apply(X = d[d$ID == unique(d$ID)[i], c("BTR", "Track", "Rail_string") ], MARGIN = 2, FUN =function(x) length(unique(x))))) save = c(save,unique(d$ID)[i])
}


newids = max(d$ID, na.rm = TRUE) + (1:length(save))

for (i in 1:length(save)) {
  tryCatch({
    check_uni = TRUE
    
    while (check_uni) {
      foo = d[d$ID == save[i], ]
      
      if (NROW(foo) < 3) { 
        d[foo$ID[1] == d$ID & foo$BTR[1] == d$BTR & foo$Track[1] == d$Track & foo$Rail_string[1] == d$Rail_string,"ID"] = newids[i]
        break
      }
      
      uni = apply(foo[, c("BTR", "Track", "Rail_string")], 2, table)
      
      if (length(uni$BTR) > 1) {
        idx_btr = which.min(uni$BTR)
        odd_btr = names(uni$BTR[idx_btr])
        btr_out = odd_btr == foo$BTR
      } else {
        btr_out = rep(FALSE, NROW(foo))
      }
      
      if (length(uni$Track) > 1) {
        idx_trk = which.min(uni$Track)
        odd_trk = names(uni$Track[idx_trk])
        trk_out = odd_trk == foo$Track
      } else {
        trk_out = rep(FALSE, NROW(foo))
      }
      
      if (length(uni$Rail_string) > 1) {
        idx_rai = which.min(uni$Rail_string)
        odd_rai = names(uni$Rail_string[idx_rai])
        rai_out = odd_rai == foo$Rail_string
      } else {
        rai_out = rep(FALSE, NROW(foo))
      }
      
      res = foo[apply(cbind(btr_out, trk_out, rai_out), 1, any), ]
      
      if (NROW(res) > 0) {
        for (j in 1:NROW(res)) {
          d[d$ID == res$ID[j] & d$BTR == res$BTR[j] & 
              d$Track == res$Track[j] & d$Rail_string == res$Rail_string[j], 
            "ID"] = newids[i]
        }
      }
      
      foo = d[d$ID == save[i], ]
      uni = apply(foo[, c("BTR", "Track", "Rail_string")], 2, table)
      check_uni = any(sapply(uni, length) > 1)
    }
    
  }, warning = function(w) {
    cat("Warning in iteration", i, ":", conditionMessage(w), "\n")
  })
}


save2 = c()
for(i in 1:length(unique(d$ID))){
  if(any(1!=apply(X = d[d$ID == unique(d$ID)[i], c("BTR", "Track", "Rail_string") ], MARGIN = 2, FUN =function(x) length(unique(x))))) save2 = c(save2,unique(d$ID)[i])
  if(i %% 1000 == 0){print(i)}
  }


dnew = d[!d$ID %in% save2, ]
#write.csv(x = dnew, row.names = F,  file = "defects_UC_EC_maintenance_covariates_v2.csv")



#############################
### CLEANING AND ENCODING ###
#############################
# rm(list = ls())
# D = read.csv(file = "defects_UC_EC_maintenance_covariates_v2.csv")
# 
# 
# #MAINTENANCE ENCODING
# D$Removed_status = D$State == "Removed" #Encoding of binary Remove control action
# D$Grinding = D$Maintenance_type %in% c("Grinding", "Grinding (HS)") 
# D$Milling = D$Maintenance_type %in% c("Milling", "Milling (HS)")
# D$Planing = D$Maintenance_type %in% c("Planing")
# D$Maintenance_indicator = 0
# D$Maintenance_indicator[!is.na(D$Maintenance_date)] = 1
# D$Passages[is.na(D$Maintenance_date)] = 0
# 
# sum(is.na(D$Maintenance_date) == (D$Grinding | D$Milling | D$Planing)) #CHECK IF MAINTENANCE DATES EXISTS WHEN MAINTENANCE IS CONDUCTED:
# 
# #TURNOUT ENCODING
# idx = D$Turnout_indicator %in% c(1) & D$Track_type %in% c("Standard_track"); D$Track_type[idx] = "Switch" #Find turnout_indications not repored in track_type
# idx = !D$Turnout_indicator %in% c(1) & !D$Track_type %in% c("Standard_track", "Unknown"); D$Turnout_indicator[idx] = 1 #Overwrite turnout indications, where track type reports turnout
# idx = D$Rail_type %in% c("Switch") & D$Turnout_indicator %in% c(0); D$Turnout_indicator[idx] = 1 #Overwrite turnout indications where rail_type is reported as switch
# 
# #WELD ENCODINGS
# D$Aluminothermic_weld = D$Defect_found_in %in% c("Aluminothermic_welding")
# D$Flash_butt_weld = D$Defect_found_in %in% c("Flash_butt_welding")
# D$Other_weld = !D$Defect_found_in %in% c("Aluminothermic_welding","Flash_butt_welding","Clean_rail")
# D$Weld = !D$Defect_found_in %in% c("Clean_rail")
# 
# #CURVE ENCODINGS
# D$In_curve = D$Track_type3 %in% c("Curve")
# D$In_trans_curve = D$Track_type3 %in% c("Transition_curve")
# D$In_straight_track = D$Track_type3 %in% c("Straigt_track")
# 
# #PROFILE ENCODINGS
# idx = (is.na(D$Profile) & is.na(D$Rail_type)) | 
#   (is.na(D$Profile) & D$Rail_type %in% c("Switch")); 
# 
# #D = D[!idx,]  #remove rows with no useful information on profile (*= 8 rows)
# D$Rail_type[is.na(D$Rail_type)] = D$Profile[is.na(D$Rail_type)] #overwrite NANs
# D$Profile[is.na(D$Profile)] = D$Rail_type[is.na(D$Profile)] #Overwrite nans
# library(stringr)
# extract_two_digits <- function(strings) {str_extract(strings, "\\d{2}")} #convert strings to digits
# D$Rail_weight_1m = extract_two_digits(D$Profile)
# 
# #STEEL ENCODINGS
# extract_three_digits <- function(strings) {str_extract(strings, "\\d{3}")}
# D$Steel_hardness = extract_three_digits(D$Steel)
# 
# #REMOVE OBS MISSING SUBSTANTIAL INFORMATION
# idx = is.na(D$UIC) | is.na(D$defect_size) | is.na(D$Date); D = D[!idx, ]
# 
# 
# D$Age = (as.numeric( as.Date(D$Date) - as.Date(D$Year)) / 365)
# 
# 
# cols1 = c("ID", "BTR", "Track", "From", "To", "Rail_string")
# cols2 = c("Date", "defect_size")
# cols3 = c("Defect_group", "Visible", "UIC","Combined_defect")
# cols4 = c("Maintenance_date","Maintenance_indicator", "Grinding", "Milling", "Planing", "Passages","Removed_status")
# cols5 = c("In_straight_track", "In_curve", "In_trans_curve","Turnout_indicator","Weld","Aluminothermic_weld","Flash_butt_weld","Other_weld")
# cols6 = c("Age", "Curve", "Rail_weight_1m", "Steel_hardness", "Line_speed", "MGT_max")
# 
# idx = rowSums(apply(D[, cols6], 2, FUN = function(x) is.na(x))) == 0; D = D[idx, ] #only keep complete covariates
# 
# 
# #Count how many times ID appears
# library(dplyr)
# id_counts <- D %>%
#   count(ID)
# 
# freq_of_freqs <- id_counts %>%
#   count(n, name = "num_ids") %>%
#   rename(times = n)
# 
# D2 <- D %>%
#   group_by(ID) %>%
#   filter(n() > 1) %>%
#   ungroup()
# 
# D2 = as.data.frame(D2)

#write.csv(x =D[, c(cols1, cols2, cols3, cols4, cols5, cols6)],row.names = F,  file = "C:/Users/askbi/OneDrive - Danmarks Tekniske Universitet/SpatioTemporal/defect_trajectories_long_format.csv")



#####################################################
### COMBINE WEAR, CORRUGATION, TWIST, INCLINATION ###
#####################################################
rm(list = ls())
library(dplyr)
file_info <- list(
  align_horiz_d0 = list(file = "Processed_meas_at_defect_locs/Align_horiz_D0_at_defect_locs.csv", prefix = "align_horiz_d0_"),
  align_horiz_d1 = list(file = "Processed_meas_at_defect_locs/Align_horiz_D1_at_defect_locs.csv", prefix = "align_horiz_d1_"),
  align_vert_d0  = list(file = "Processed_meas_at_defect_locs/Align_vert_D0_at_defect_locs.csv", prefix = "align_vert_d0_"),
  align_vert_d1  = list(file = "Processed_meas_at_defect_locs/Align_vert_D1_at_defect_locs.csv", prefix = "align_vert_d1_"),
  cant           = list(file = "Processed_meas_at_defect_locs/Cant_at_defect_locs.csv", prefix = "cant_"),
  corrugation    = list(file = "Processed_meas_at_defect_locs/Corrugation_at_defect_locs.csv", prefix = "corrugation_"),
  gauge          = list(file = "Processed_meas_at_defect_locs/Gauge_at_defect_locs.csv", prefix = "gauge_"),
  height_wear    = list(file = "Processed_meas_at_defect_locs/Height_Wear_at_defect_locs.csv", prefix = "wear_h_"),
  inclination    = list(file = "Processed_meas_at_defect_locs/Inclination_at_defect_locs.csv", prefix = "inclination_"),
  twist2m        = list(file = "Processed_meas_at_defect_locs/Twist2m_at_defect_locs.csv", prefix = "twist2m_"),
  twist3m        = list(file = "Processed_meas_at_defect_locs/Twist3m_at_defect_locs.csv", prefix = "twist3m_"),
  width_wear     = list(file = "Processed_meas_at_defect_locs/Width_Wear_at_defect_locs.csv", prefix = "wear_w_")
)

dfs <- lapply(file_info, function(info) {
  df <- read.csv(info$file, stringsAsFactors = FALSE)
  colnames(df)[8:12] <- paste0(info$prefix, colnames(df)[8:12])
  df <- df[!duplicated(df[c("ID", "Date")]), ]
  df
})
merge_order <- c("height_wear", "width_wear", "corrugation", "inclination", "twist2m", "twist3m", "gauge", "cant", "align_horiz_d0", "align_horiz_d1", "align_vert_d0", "align_vert_d1")
smart_merge <- function(df1, df2) {
  redundant_cols <- c("BTR", "Track", "Rail_string", "From", "To")
  df <- full_join(df1, df2, by = c("ID", "Date"))
  df <- df %>%
    mutate(across(ends_with(".x"), ~ ifelse(is.na(.), get(gsub("\\.x$", ".y", cur_column())), .))) %>%
    select(-all_of(paste0(redundant_cols, ".y"))) %>%
    rename_with(~ gsub("\\.x$", "", .x), ends_with(".x")) %>%
    arrange(ID, Date)
  return(df)
}

df_combined <- dfs[[merge_order[1]]]
for (i in merge_order[-1]) {
  df_combined <- smart_merge(df_combined, dfs[[i]])
  message("Rows after merging ", i, ": ", NROW(df_combined))
}

write.csv(df_combined, file = "track_measurements_at_defects.csv", row.names = FALSE)


#################################
### DEFECTS: FURTHER CLEANING ###
#################################
rm(list = ls())
library(dplyr)
de = read.csv("defects_UC_EC_maintenance_covariates_v2.csv")

#propagates obs with UIC = 2223 to all rows under an ID if all UT observations under that ID have UIC == 2223
mixedIDs = unique(de$ID[de$Cluster_method=="EC+UT"])
for(i in mixedIDs){
  idx = de$ID == i
  idx2 = de$Method[idx] == "UT"
  if( all(de[idx,][idx2,"UIC" ]%in% c(2223) ) ){
   de$UIC[idx] = 2223 
  }
}

#UIC to defect type encoding:
de <- de %>%
  mutate(
    UIC_group = case_when(
      UIC %in% c(211, 227)                  ~ "Squats",
      UIC %in% c(113, 213)                  ~ "Vertical longitudinal",
      UIC %in% c(112, 212, 232, 412)        ~ "Horizontal longitudinal",
      UIC %in% c(411, 421, 422, 471, 472)   ~ "Weld defects",
      UIC %in% c(2223)                      ~ "Head checks",
      UIC %in% c(135, 235)                  ~ "Fish plate chamber",
      TRUE                                  ~ "Other/Unknown"
    ),
    UIC_group_id = case_when(
      UIC_group == "Squats"                    ~ 1,
      UIC_group == "Vertical longitudinal"     ~ 2,
      UIC_group == "Horizontal longitudinal"   ~ 3,
      UIC_group == "Weld defects"              ~ 4,
      UIC_group == "Head checks"               ~ 5,
      UIC_group == "Fish plate chamber"        ~ 6,
      TRUE                                     ~ NA_integer_
    )
  )

#classify EC only measurements to headchecks if they appear on outer rail in curve
de <- de %>%
  group_by(ID) %>%
  mutate(
    UIC_group_id = case_when(
      all(Cluster_method == "EC", na.rm = TRUE) &
        all(Rail_string == "Left", na.rm = TRUE) &
        all(Curve > 0, na.rm = TRUE) ~ 5,
      
      all(Cluster_method == "EC", na.rm = TRUE) &
        all(Rail_string == "Right", na.rm = TRUE) &
        all(Curve < 0, na.rm = TRUE) ~ 5,
      
      TRUE ~ UIC_group_id  # keep existing if not matched
    )
  ) %>%
  ungroup()

de$size = de$Depth_to - de$Depth_from

#remove all rows with nan uic and all defect ID where defect classification is not matching
de_clean <- de %>%
  filter(!is.na(UIC_group_id)) %>%         
  filter(!is.na(size)) %>%
  filter(size > 0) %>%
  group_by(ID) %>%
  filter(n_distinct(UIC_group_id) == 1) %>%  #keep only consistent IDs
  ungroup()


d <- de_clean %>%
  select(ID, Date, size, UIC_group_id, everything())
d = as.data.frame(d)
#head(d)
#colSums(is.na(de_clean))
#head(as.data.frame(de_clean))

write.csv(d, file = "defects_UC_EC_maintenance_covariates_cleaned.csv", row.names = FALSE)

#####################################
### MAINTENANCE DATA ID ASSIGNING ### 
#####################################
rm(list = ls())
library(dplyr)
m = read.csv("rail_maintenance.csv")
d = read.csv("defects_UC_EC_maintenance_covariates_cleaned.csv")
m <- m %>%
  rename(Date = Time) %>%
  mutate(source = "maintenance")
M <- list() 
for (i in 1:nrow(m)) {
  idx <- d$BTR == m$BTR[i] & d$Track == m$Track[i]
  if (!any(idx)) next
  d_sub <- d[idx, ]
  overlap <- d_sub$From <= m$To[i] & d_sub$To >= m$From[i]
  if (!any(overlap)) next
  matched_IDs <- unique(d_sub$ID[overlap])
  m_rep <- m[rep(i, length(matched_IDs)), ]
  m_rep$ID <- matched_IDs
  M[[length(M) + 1]] <- m_rep
}
M_df <- bind_rows(M)
M_df$Grinding = M_df$action %in% c("Grinding", "Grinding (HS)") 
M_df$Milling = M_df$action %in% c("Milling", "Milling (HS)")
M_df$Planing = M_df$action %in% c("Planing")

write.csv(M_df, file = "rail_maintenance_with_defect_IDs.csv", row.names = FALSE)


###############################################
### Combine track measurements with defects ###
###############################################
rm(list = ls())
library(dplyr)
t = read.csv("track_measurements_at_defects.csv")
de = read.csv("defects_UC_EC_maintenance_covariates_cleaned.csv")
m = read.csv("rail_maintenance_with_defect_IDs.csv")

de = de[, !colnames(de)%in%c("Maintenance_date","Maintenance_type","Amount","Passages")]
all_cols <- union(names(de), names(t))

de_aligned <- de %>%
  mutate(source = "defect") %>%
  tibble::add_column(!!!setNames(
    rep(list(NA), length(setdiff(all_cols, names(.)))),
    setdiff(all_cols, names(.))
  )) %>%
  select(all_of(all_cols), source)
t_aligned <- t %>%
  filter(ID %in% de$ID) %>%               # restrict to relevant IDs
  mutate(source = "track") %>%
  tibble::add_column(!!!setNames(
    rep(list(NA), length(setdiff(all_cols, names(.)))),
    setdiff(all_cols, names(.))
  )) %>%
  select(all_of(all_cols), source)
rm(de, t)

combined_df <- bind_rows(de_aligned, t_aligned) %>%
  arrange(ID, Date)

# combined_df <- combined_df %>%
#   mutate(Date = as.Date(Date))
# m_aligned <- m %>%
#   mutate(Date = as.Date(Date),  # if not already Date type
#          source = "maintenance") 


all_cols <- union(names(combined_df), names(m))  # recompute to include new columns from m

m_aligned <- m %>%
  tibble::add_column(!!!setNames(
    rep(list(NA), length(setdiff(all_cols, names(.)))),
    setdiff(all_cols, names(.))
  )) %>%
  select(all_of(all_cols))
rm(m)

combined_df <- bind_rows(combined_df, m_aligned) %>%
  arrange(ID, Date)


write.csv(combined_df, file = "defects_UC_EC_maintenance_covariates_cleaned_track_measurements_v2.csv", row.names = FALSE)

################
### CLEANING ###
################
rm(list = ls())
library(dplyr)
d = read.csv("defects_UC_EC_maintenance_covariates_cleaned_track_measurements_v2.csv")
#d = d_copy
#d_copy = d

#fill out missing static information on ID level
cols_to_fill <- c("Rail_string", "Switch_nr.", "Track_type", "Profile", "Curvature", "Year", "Steel", 
                  "Fastening", "Switch_type", "Rail_type2", "Track_type2", "Steel_type2", "Track_type3", "Curve",
                  "Overheight", "MGT_min","MGT_max", "EMGT_min", "EMGT_max", "Line_speed", "Turnout_indicator")
d <- d %>%
  group_by(ID) %>%
  mutate(across(all_of(cols_to_fill), ~ {
    non_na <- .[!is.na(.)]
    if (length(non_na) == 0) return(.)            # all NA → leave unchanged
    replace(., is.na(.), non_na[1])               # replace NA with first non-NA
  })) %>%
  ungroup()
d = as.data.frame(d)
get_rail_profile <- function(x) {
  as.numeric(sub(".*?(\\d+).*", "\\1", x))
}
#overwrite untrustworthy year information
d <- d %>% mutate(Year = ifelse(Year < 1930, NA, Year))
d$rail_age = as.numeric(format(as.Date(d$Date), "%Y")) - d$Year
d <- d %>% mutate(rail_age = ifelse(rail_age < 0, NA, rail_age))
#Encoding´
d <- d %>%
  mutate(
    in_curve            = as.integer(Track_type3 == "Curve"),
    in_straight         = as.integer(Track_type3 == "Straight_track"),
    in_transition_curve = as.integer(Track_type3 == "Transition_curve"),
    maintenance_source = as.integer(source == "maintenance"),
    defect_source         = as.integer(source == "defect"),
    track_source = as.integer(source == "track"), 
    rail_left = as.integer(Rail_string == "Left"), 
    visible = as.integer(Visible == "Yes"), 
    combined_defect = as.integer(Combined_defect == "Yes"), 
    under_limit = as.integer(Under_limit == "Yes"), 
    R200 = as.integer(Steel == "R200"), 
    R260 = as.integer(Steel == "R260"), 
    R350HT = as.integer(Steel == "R350HT"), 
    grinding = as.integer(Grinding), 
    milling = as.integer(Milling), 
    planing = as.integer(Planing), 
    weight = get_rail_profile(Profile)
  )
#Cleaning
remove.cols = c("Switch_nr.", "Defect_length", "Defect_width", "Defect_found_in", "Track_type", 
                "Track_type2", "Switch_type", "Defect_group_calc", "UIC", "Rail_type2","Reason", "action", 
                "Track_type3", "source", "Rail_string", "Visible", "Combined_defect", "Under_limit", "Steel", "Steel_type2",
                "Curvature", "Grinding", "Milling", "Planing", "Year", "Depth_from", "Depth_to", "amount")
d = d[, !colnames(d) %in% remove.cols]
d$UIC_group[d$UIC_group=="Other/Unknown"] = "Head checks"
#compute dTime
d <- d %>%
  arrange(ID, Date) %>%
  group_by(ID) %>%
  mutate(
    dTime = as.numeric(c(0, diff(as.Date(Date))))
  ) %>%
  ungroup()
d = as.data.frame(d)
#indicate complete exo. information
cols = c("weight", "R200","R260", "R350HT","MGT_max", "Curve", "Line_speed", "rail_age")
d$covariates_complete = as.integer(!apply(d[, cols], 1 , function(x) (any(is.na(x)))))
#reorder
d = (d %>% select(ID, Date, dTime, size, UIC_group_id,UIC_group, in_straight,in_transition_curve,in_curve, covariates_complete,rail_age,R200,R260,R350HT,grinding,milling,planing,weight,Curve,Overheight,
                  everything(),
                  -c(BTR,Track,From,To,Defect_group,Profile,Fastening,State,Method,Cluster_method), 
                  c(BTR,Track,From,To,Defect_group,Profile,Fastening,State,Method,Cluster_method),) ) 
#rename
d <- d %>%
  rename(
    crack_size = size,
    crack_class = UIC_group_id,
    crack_class_description = UIC_group,
    rail_1m_weight = weight,
    passages_maintenance = Passages,
    stone_maintenance = stone,
    removed_material = removed_mm,
    detection_method = Method, 
    cluster_method = Cluster_method
  )
d <- d %>%
  group_by(ID) %>%
  mutate(crack_counts = sum(!is.na(crack_size))) %>%
  ungroup()
d = as.data.frame(d)


write.csv(d, file = "crack_data.csv", row.names = FALSE)



#################################################
### Incorporate traffic direction in data set ###
#################################################
rm(list = ls())
library(readxl)
library(dplyr)
library(purrr)

d = read.csv("crack_data.csv")
d$traffic_direction = NA


f = read_excel("direction_train.xlsx")
f = as.data.frame(f)
colnames(f)[5] = "Text1"

f <- f %>% arrange(BTR, Spor, From, To)
resolve_group <- function(group) {
  #H: ascending km
  #V: descending km
  #E: both
  #O: Overlap of ascending and descending category, meaning either both or unknown.
  resolved <- list()
  for (i in seq_len(nrow(group))) {
    row <- group[i, ]
    start <- row$From
    end <- row$To
    label <- row$Text1
    inserted <- FALSE
    # Look for overlaps in the current resolved list
    new_resolved <- list()
    for (seg in resolved) {
      r_start <- seg$From
      r_end   <- seg$To
      r_label <- seg$Text1
      if (end <= r_start || start >= r_end) {
        new_resolved <- append(new_resolved, list(seg))
      } else {
        # Overlapping case
        overlap_start <- max(start, r_start)
        overlap_end   <- min(end, r_end)
        # Before overlap
        if (r_start < overlap_start)
          new_resolved <- append(new_resolved, list(data.frame(BTR=row$BTR, Spor=row$Spor, From=r_start, To=overlap_start, Text1=r_label)))
        if (start < overlap_start)
          new_resolved <- append(new_resolved, list(data.frame(BTR=row$BTR, Spor=row$Spor, From=start, To=overlap_start, Text1=label)))
        # Overlap
        new_label <- if (r_label == label) r_label else "O"
        new_resolved <- append(new_resolved, list(data.frame(BTR=row$BTR, Spor=row$Spor, From=overlap_start, To=overlap_end, Text1=new_label)))
        # After overlap
        if (r_end > overlap_end)
          new_resolved <- append(new_resolved, list(data.frame(BTR=row$BTR, Spor=row$Spor, From=overlap_end, To=r_end, Text1=r_label)))
        if (end > overlap_end)
          new_resolved <- append(new_resolved, list(data.frame(BTR=row$BTR, Spor=row$Spor, From=overlap_end, To=end, Text1=label)))
        inserted <- TRUE
      }
    }
    if (!inserted) {
      new_resolved <- append(new_resolved, list(row))
    }
    resolved <- new_resolved
  }
  bind_rows(resolved) %>% arrange(From, To)
}
f <- f %>% group_split(BTR, Spor) %>% map_dfr(resolve_group)


average_overlap <- function(a_from, a_to, b_from, b_to) {
  len_a <- a_to - a_from
  len_b <- b_to - b_from
  overlap <- max(0, min(a_to, b_to) - max(a_from, b_from))
  if (len_a <= 0 || len_b <= 0) return(0)
  return(overlap / len_a)  # Use only the fraction of overlap relative to A
}

ids <- unique(d$ID)

for (i in seq_along(ids)) {
  level <- d[d$ID == ids[i], c("BTR", "Track", "From", "To")]
  candidates <- f[f$BTR == level$BTR[1] & f$Spor == level$Track[1], ]
  if (nrow(candidates) < 1) next
  overlaps <- numeric(nrow(candidates))
  for (j in seq_len(nrow(candidates))) overlaps[j] <- average_overlap(min(level$From), max(level$To),candidates$From[j], candidates$To[j])
  if (all(overlaps == 0)) next
  best_match <- which.max(overlaps)
  d[d$ID == ids[i], "traffic_direction"] <- candidates$Text1[best_match]
}

#ENCODING directions
d$traffic_direction <- dplyr::case_when(
  d$traffic_direction == "H" ~ 1L,
  d$traffic_direction == "V" ~ -1L,
  d$traffic_direction == "E" ~ 0L,
  d$traffic_direction == "O" ~ NA_integer_, #overlaps are set to NaN
  is.na(d$traffic_direction) ~ NA_integer_,
  TRUE ~ NA_integer_
)

#update complete covariate information
d = as.data.frame(d)
cols = c("rail_1m_weight", "R200","R260", "R350HT","MGT_max", "Curve", "Line_speed", "rail_age","traffic_direction")
d$covariates_complete = as.integer(!apply(d[, cols], 1 , function(x) (any(is.na(x)))))
#reorder
allcols <- colnames(d)
allcols <- setdiff(allcols, "traffic_direction")
insert_after <- match("rail_1m_weight", allcols)
new_order <- append(allcols, "traffic_direction", after = insert_after)
d <- d[, new_order]

write.csv(d, file = "crack_data_v2.csv", row.names = FALSE)

#####################
### ANONYMIZATION ###
#####################

rm(list = ls())
library(dplyr)

d = read.csv("crack_data_v2.csv")

d$BTR_first2digits <- as.numeric(substr(sprintf("%06d", d$BTR), 1, 2)) # 1. Extract the first two digits of BTR
btr_levels <- unique(d$BTR)
btr_map <- setNames(sprintf("BTR%03d", seq_along(btr_levels)), btr_levels)
d$BTR <- btr_map[as.character(d$BTR)]
track_levels <- unique(d$Track)
track_map <- setNames(sprintf("TRK%03d", seq_along(track_levels)), track_levels)
d$Track <- track_map[d$Track]

write.csv(d, file = "crack_data_v2_anonymized.csv", row.names = FALSE)




# - remember to combine rows in df, that share the same date and ID
# - remove nans where possible (e.g. static information contained elsewhere in the df)
# - is tonnage time-varying? Then take this into account
# - remove unneeded columns
# - make sure all observation are complete in terms of needed information
# - streamline NA and <NA> accross variables
# - naming schemes (translation where needed, e.g. maintenance Reason)
# - make description of data
# - indicators on ID level, so data is easy to filter (e.g. on defect types, number of observation per defect, complete covariates)
# - check if dates are the same format accross all obs






