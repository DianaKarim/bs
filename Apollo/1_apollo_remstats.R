library(readxl)
library(fastDummies)
library(abind)

#REMverse packages
library(remify)
library(remstats)
#To install use this:
#library(devtools)
#devtools::install_github(repo = "TilburgNetworkGroup/remify", build_vignettes = TRUE)
#devtools::install_github("TilburgNetworkGroup/remstats")

# Complete data:
# Load data
load("apollo13.RData")
apollo13 <- Merged_ALL_parts_Apollo
names(apollo13)[2:3] <- c("actor1", "actor2")

# Select the events after the incident
apollo13 <- apollo13[97:nrow(apollo13),]

# Actors
actors <- sort(unique(c(apollo13$actor1, apollo13$actor2)))

# Information on the actors
info <- data.frame(
  id = actors,
  time = 0, 
  team = ifelse(actors %in% c("CDR", "LMP", "CMP"), "air", "ground")
)

# Add dummy variables 
info <- dummy_cols(info, select_columns = "id")
info <- dummy_cols(info, select_columns = "team")

# Prepare the edgelist
reh_apollo13 <- reh(apollo13[,c(1:3)], model = "tie", directed = TRUE, actors = actors)

# Get the riskset
riskset <- remstats(reh_apollo13, tie_effects = ~ 1)$riskset

# Add structural zero's (modification of the risk set according to the protocols of the communication)
riskset$allowed <- 1
riskset$allowed[riskset$sender %in% info$id[info$team == "ground" & info$id != "CAPCOM"] &
                  riskset$receiver %in% info$id[info$team == "air"]] <- 0
riskset$allowed[riskset$receiver %in% info$id[info$team == "ground" & info$id != "CAPCOM"] &
                  riskset$sender %in% info$id[info$team == "air"]] <- 0

# Tie-effects
getTie <- function(var1, var2, actors) {
  x <- var1 %*% t(var2)
  rownames(x) <- colnames(x) <- actors
  x
}

#ground_to_ground <- getTie(info$team_ground, info$team_ground, info$id)
air_to_air <- getTie(info$team_air, info$team_air, info$id)
ground_to_CAPCOM <- getTie(info$team_ground, info$id_CAPCOM, info$id)
air_to_CAPCOM <- getTie(info$team_air, info$id_CAPCOM, info$id)
CAPCOM_to_ground <- getTie(info$id_CAPCOM, info$team_ground, info$id)
CAPCOM_to_air <- getTie(info$id_CAPCOM, info$team_air, info$id)
ground_to_FLIGHT <- getTie(info$team_ground, info$id_FLIGHT, info$id)
FLIGHT_to_ground <- getTie(info$id_FLIGHT, info$team_ground, info$id)

# Model 1 
effects <- ~ -1 + 
  (rrankSend() + rrankReceive() +
     psABBA() + psABBY() + psABXA() + psABAY() +
     inertia(scaling = "prop") + reciprocity(scaling = "prop") +
     outdegreeSender(scaling = "prop") + indegreeReceiver(scaling = "prop")) :
  (#tie(ground_to_ground, "ground_to_ground") + 
    tie(air_to_air, "air_to_air") +
      tie(ground_to_CAPCOM, "ground_to_CAPCOM") + tie(air_to_CAPCOM, "air_to_CAPCOM") +
      tie(CAPCOM_to_ground, "CAPCOM_to_ground") + tie(CAPCOM_to_air, "CAPCOM_to_air") +
      tie(ground_to_FLIGHT, "ground_to_FLIGHT") + tie(FLIGHT_to_ground, "FLIGHT_to_ground"))

stats1 <- remstats(reh_apollo13, tie_effects = effects)$statistics

# Prepare data.frame for estimation
dictionary <- attributes(reh_apollo13)$dictionary$actors
Events_apollo <- apollo13[,c(2,3)]
colnames(Events_apollo) <- c("sender", "receiver")
Events_apollo$index <- reh_apollo13$edgelist[,2] + 1
Events_apollo$sender <- match(Events_apollo$sender, dictionary$actorName)
Events_apollo$receiver <- match(Events_apollo$receiver, dictionary$actorName)

# Deal with structural zero's
head(riskset)
riskset$new_stat_column <- cumsum(riskset$allowed)
riskset$new_stat_column[duplicated(riskset$new_stat_column)] <- NA

# Remove these columns from the statistics
stats1 <- stats1[,-which(riskset$allowed==0),]

# Replace the Events_index
Events_apollo$index <- riskset$new_stat_column[match(Events_apollo$index, riskset$stat_column)]

# Function to standardize the statistics
standardizeStats <- function(stats) {
  std_stats <- lapply(1:dim(stats)[3], function(i) {
    (stats[,,i]-mean(stats[,,i]))/sd(stats[,,i])
  })
  std_stats <- do.call(abind, args = list(std_stats, along = 3))
  apply(std_stats, 3, sd) 
  dimnames(std_stats) <- dimnames(stats)
  std_stats
}
stats1 <- standardizeStats(stats1)

save(stats1, Events_apollo, file = "apollo_estimation_data1.RData")

# Model 2 (used in the reviewed paper, 103 statistics)
effects <- ~ -1 + 
  (rrankSend() + rrankReceive() +
     psABBA() + psABBY() + psABXA() + psABAY() +
     inertia(scaling = "prop") + reciprocity(scaling = "prop") +
     outdegreeSender(scaling = "prop") + indegreeReceiver(scaling = "prop") +
     otp() + itp()) :
  (tie(air_to_air, "air_to_air") +
      tie(ground_to_CAPCOM, "ground_to_CAPCOM") + tie(air_to_CAPCOM, "air_to_CAPCOM") +
      tie(CAPCOM_to_ground, "CAPCOM_to_ground") + tie(CAPCOM_to_air, "CAPCOM_to_air") +
      tie(ground_to_FLIGHT, "ground_to_FLIGHT") + tie(FLIGHT_to_ground, "FLIGHT_to_ground"))

stats2 <- remstats(reh_apollo13, tie_effects = effects)$statistics

# Remove appropriate columns from the statistics
stats2 <- stats2[,-which(riskset$allowed==0),]

# Standardize
stats2 <- standardizeStats(stats2)

# Save
save(stats2, Events_apollo, file = "apollo_estimation_data2.RData")


# Some descriptive statistics 

nrow(apollo13[apollo13$actor1 == "FLIGHT",])/nrow(apollo13)*100
nrow(apollo13[apollo13$actor2 == "FLIGHT",])/nrow(apollo13)*100

nrow(apollo13[apollo13$actor1 == "CAPCOM",])/nrow(apollo13)*100
nrow(apollo13[apollo13$actor2 == "CAPCOM",])/nrow(apollo13)*100



corr(stats2)









