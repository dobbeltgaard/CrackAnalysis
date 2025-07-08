
rm(list = ls())
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

d = read.csv("crack_data_v2.csv")


### EC vs UT distribution ###
table(d[!is.na(d$detection_method),"detection_method" ])

### Counting trajectories per crack_class ###
counts = 
  d %>% group_by(ID)%>% filter(!is.na(crack_class)) %>%
  group_by(ID, crack_class) %>%
  summarise(n = n(), .groups = "drop") %>%         # how many times each cluster_method appears per ID
  group_by(crack_class, n) %>%
  summarise(n_IDs = n(), .groups = "drop") %>% as.data.frame()
counts_wide <- counts %>% pivot_wider(names_from = crack_class,values_from = n_IDs,values_fill = 0)
counts_wide$Total = rowSums(counts_wide[2:7])

counts_wide %>%
  kable(format = "latex", booktabs = TRUE, 
        col.names = c("n", "Squats", "Vertical longitudinal", "Horizontal longitudinal", "Weld defects", "Head checks", "Fish plate chamber", "Total"),
        caption = "Counts per severity class and number of observations") %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))


sum((counts_wide$`5`*counts_wide$n)[2:12])


### Covariates 


### Number of track measurements per defect
d %>% group_by(ID) %>% summarise(n = n(), .groups = "drop")



#how to analyse correlation of a dependent variable, say crack size, as a function of track measurements at the location of the crack? I have a dataset with m






#Check Violations functionality
check_id_violations <- function(data, var) {
  data %>%
    group_by(ID) %>%
    summarize(n_unique = n_distinct(.data[[var]][!is.na(.data[[var]])]), .groups = "drop") %>%
    filter(n_unique > 1)
}

check_id_violations(d, "crack_class")










sum(cluster_histogram$n*cluster_histogram$n_IDs)


sum(!is.na(d$cluster_method))

foo = (d[!is.na(d$cluster_method), ])


foo2 <- foo %>%
  group_by(ID) %>%
  filter(n_distinct(cluster_method) > 1) %>%
  arrange(Date) %>%
  slice(1) %>%
  ungroup()


foo = as.data.frame(foo)
idx = d$ID %in% foo$ID
d$Cluster_method = d$Method
d$Cluster_method[idx] = "EC+UT"

#Summary of sequence size per cluster method
mat = matrix(NA, 12,4)
mat[1:length(table(table(d$ID[d$cluster_method == "EC"]))), 1] = table(table(d$ID[d$cluster_method == "EC"]))
mat[1:length(table(table(d$ID[d$cluster_method == "UT"]))), 2] = table(table(d$ID[d$cluster_method == "UT"]))
mat[1:length(table(table(d$ID[d$cluster_method == "EC+UT"]))), 3] = table(table(d$ID[d$cluster_method == "EC+UT"]))

mat[2:(1+length(table(table(d$ID[d$cluster_method == "EC+UT"])))), 3] = table(table(d$ID[d$cluster_method == "EC+UT"]))
mat[1:length(table(table(d$ID))), 4] = rowSums(mat, na.rm = T)
rownames(mat) = 1:12
colnames(mat) = c("EC", "UT", "EC+UT", "Total")



# cluster_count_table <- d %>%
#   filter(!is.na(cluster_method)) %>%
#   group_by(ID, cluster_method, crack_class) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   count(n, name = "n_IDs")
# print(cluster_count_table)
# sum(cluster_count_table$n*cluster_count_table$n_IDs)



