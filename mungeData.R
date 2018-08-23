## Moving all the data munging into the same script
# M. Espe & M. Johnston
# 2017 March

########## Munging the full data set ##########

load("data/input_matrices_and_vectors_plusreaches.Rdata")

tagdata <- read.csv("data/tagdata_index.csv", stringsAsFactors = FALSE)

# need to collapse the extra columns in the matrices
sac2012_matrix <- cbind(sac2012_matrix, 0)
colnames(sac2012_matrix)[12] <- "39" # Receiver #39 = phantom

routeSac12 = cbind(c(sac2012_RouteAvector[1:5], 39,sac2012_RouteAvector[6:8]),
                   sac2012_RouteBvector,
                   c(sac2012_RouteCvector[1:5], 39,sac2012_RouteCvector[6:8]))
group12 = as.integer(as.factor(
  tagdata$TrueRoute[match(rownames(sac2012_matrix), tagdata$TagID)]))

Sac12 <- do.call(rbind, lapply(seq_along(group12), function(i){
  if(group12[i] != 4)
    j <- group12[i]
  else
    j <- 1
  sac2012_matrix[i , as.character(routeSac12[,j])]
}))

# Repeat - turn into function later

sac2013_matrix <- cbind(sac2013_matrix, 0)
colnames(sac2013_matrix)[18] <- "39"

routeSac13 = cbind(sac2013_RouteAvector,
                   c(sac2013_RouteBvector[1:11], 39, sac2013_RouteBvector[12:14]),
                   c(sac2013_RouteCvector[1:11], 39, sac2013_RouteCvector[12:14]))

group13 = as.integer(as.factor(
  tagdata$TrueRoute[match(rownames(sac2013_matrix), tagdata$TagID)]))

Sac13 <- do.call(rbind, lapply(seq_along(group13), function(i){
  if(group13[i] != 4)
    j <- group13[i]
  else
    j <- 1
  sac2013_matrix[i , as.character(routeSac13[,j])]
}))

# Source in the script to get the reaches by unique ID
source("makeRouteID.R")


### Route hyperparameters ####
# Hacking together by hand
hypers = list(a13=c(14, 15, 16, 17, 18, 19, 20,
                    21, 22, 23, 24, 25, 26, 27),
              b13= c(29, 30, 28),  
              c13= c(31, 32),  
              a12=c(1, 2, 3, 4, 5, 6, 7, 8),
              b12=c(9, 10, 11),
              c12=c(12, 13),
              y12 = reachYolo12,
              y13 = reachYolo13)
hyper_groups = integer(length(levs))
for(i in seq_along(hypers))
  hyper_groups[hypers[[i]]] = i


#  Create objects out of FL

sac12fls <- c(118L, 127L, 120L, 129L, 134L, 125L, 138L, 144L, 130L, 139L, 
              121L, 123L, 110L, 114L, 100L, 95L, 144L, 129L, 131L, 118L, 122L, 
              145L, 133L, 137L, 115L, 133L, 125L, 133L, 134L, 115L, 110L, 133L, 
              108L, 144L, 125L, 135L, 141L)

sac13fls <- c(111L, 156L, 142L, 161L, 130L, 135L, 131L, 117L, 136L, 130L, 
              127L, 126L, 125L, 115L, 126L, 130L, 118L, 113L, 145L, 121L, 181L, 
              135L, 127L, 127L, 120L, 138L, 139L, 120L, 137L, 123L, 128L, 121L, 
              121L, 159L, 147L, 118L, 120L, 137L, 128L, 116L, 112L, 120L, 118L, 
              122L, 123L, 126L, 152L, 131L, 123L, 126L, 137L, 142L, 110L, 129L, 
              128L, 106L, 130L, 120L, 157L, 128L, 112L, 148L, 122L, 121L, 125L, 
              131L, 111L, 152L, 158L, 132L, 125L, 136L, 147L, 166L, 122L, 139L, 
              116L, 124L, 130L, 144L, 133L, 148L, 115L, 128L, 125L, 140L, 138L, 
              160L, 108L, 135L, 153L, 123L, 135L, 135L, 118L, 140L, 158L, 128L, 
              142L, 130L)

yolo12fls <- c(127L, 108L, 136L, 113L, 119L, 117L, 128L, 118L, 116L, 107L, 
               139L, 121L, 141L, 143L, 137L, 145L, 133L, 118L, 122L, 134L, 135L, 
               132L, 108L, 138L, 127L)

yolo13fls <- c(114L, 110L, 119L, 108L, 119L, 125L, 123L, 111L, 113L, 114L, 
               113L, 120L, 115L, 120L, 114L, 123L, 117L, 107L, 107L, 125L, 119L, 
               116L, 112L, 121L, 115L)



# Make a big list with all the needed data
mod_data <- list(n_rec = length(unique(ParameterIndex$ParameterIndex)) + 1,
                 n_reach = length(levs),
                 nSac12 = nrow(Sac12),
                 nYolo12 = nrow(yolo2012_matrix),
                 nObsSac12 = ncol(Sac12),
                 nObsYolo12 = ncol(yolo2012_matrix),
                 Sac12 = Sac12,
                 Yolo12 = yolo2012_matrix,
                 routeSac12 = routeSac12,
                 reachSac12 = reachSac12,
                 routeYolo12 = yolo2012_vector,
                 reachYolo12 = reachYolo12,
                 group12 = group12,
                 nSac13 = nrow(Sac13),
                 nYolo13 = nrow(yolo2013_matrix),
                 nObsSac13 = ncol(Sac13),
                 nObsYolo13 = ncol(yolo2013_matrix),
                 Sac13 = Sac13,
                 Yolo13 = yolo2013_matrix, 
                 routeSac13 = routeSac13,
                 reachSac13 = reachSac13,
                 routeYolo13 = yolo2013_vector,
                 reachYolo13 = reachYolo13,
                 group13 = group13,
                 # Scale to make sampling easier
                 reachKm = reach_length / 10,
                 hyper_group = hyper_groups,
                 FL_Yolo12 = (yolo12fls - 130)/ 10,
                 FL_Yolo13 = (yolo13fls - 130)/ 10,
                 FL_Sac12 = (sac12fls - 130)/10,
                 FL_Sac13 = (sac13fls - 130)/10)
