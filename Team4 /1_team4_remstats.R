rm(list=ls())

library(remdata)
library(remify)
library(remstats)
library(abind)

setwd("E:/Users/mlbosman/Desktop/Revision Diana/run_models_team4")

# Get data
data(team4)
team4_edgelist <- team4$edgelist
team4_edgelist <- subset(team4_edgelist, sensor == "bluetooth")
team4_attributes  <- team4$attributes
advice <- team4$advice
social <- team4$social

# Time interval
max(team4_edgelist$time)-min(team4_edgelist$time)
hist(team4_edgelist$time, breaks = "hours")
table(lubridate::day(team4_edgelist$time))

# Subset: second day 
an.subset <- min(which(lubridate::day(team4_edgelist$time) == 18)):
	max(which(lubridate::day(team4_edgelist$time) == 18))

# Prepare data
reh_team4 <- reh(team4_edgelist[,c(1:3)], model = "tie", directed = FALSE, actors = team4_attributes$id) 

# Compute statistics
# Endogenous stats: inertia, shared partners, totaldegreeDyad, recencyContinue (4)
# Exogenous stats: same gender, age difference, tenure difference, advice seeking,
# friendship (5) 
# Interactions: 4 x 5 = 20? 

# Prepare attributes
team4_attributes$gender <- as.numeric(team4_attributes$gender)
team4_attributes$hqual <- as.numeric(team4_attributes$hqual)

# Transform advice and social to mean advice and social 
transform.to.mean <- function(x) {
	y <- x
	for(i in 1:nrow(x)) {
		for(j in 1:ncol(x)) {
			y[i,j] <- y[j,i] <- mean(c(x[i,j], x[j,i]))
		}
	}
	y
}
meanAdvice <- transform.to.mean(advice)
meanSocial <- transform.to.mean(social)

# Transform advice and social to difference in advice and social 
transform.to.diff <- function(x) {
	y <- x
	for(i in 1:nrow(x)) {
		for(j in 1:ncol(x)) {
			y[i,j] <- y[j,i] <- abs(x[i,j]-x[j,i])
		}
	}
	y
}
diffAdvice <- transform.to.diff(advice)
diffSocial <- transform.to.diff(social)

# Effects
effects <- ~ -1 + (inertia(scaling = "std") + sp(scaling = "std") +
									 	totaldegreeDyad(scaling = "std") + recencyContinue()):
	(same("gender") + difference("age", absolute = TRUE) + 
	 	difference("tenure", absolute = TRUE) + same("hqual") +
	 	tie(meanAdvice, "meanAdvice") + tie(meanSocial, "meanSocial") +
	 	tie(diffAdvice, "diffAdvice") + tie(diffSocial, "diffSocial"))

stats <- remstats(edgelist = reh_team4, tie_effects = effects, 
									attributes = team4_attributes, directed = FALSE,
									start = min(an.subset), stop = max(an.subset))$statistics

# Standardize the statistics (!!! ask whether this is the correct way to do it??)
std_stats <- lapply(1:dim(stats)[3], function(i) {
	(stats[,,i]-mean(stats[,,i]))/sd(stats[,,i])
})
std_stats <- do.call(abind, args = list(std_stats, along = 3))
apply(std_stats, 3, sd) 
dimnames(std_stats) <- dimnames(stats)

# Prepare data.frame for estimation
dictionary <- attributes(reh_team4)$dictionary$actors
Events_team4 <- team4_edgelist[an.subset,c(2,3)]
colnames(Events_team4) <- c("sender", "receiver")
Events_team4$index <- reh_team4$edgelist[an.subset,2] + 1
Events_team4$sender <- match(Events_team4$sender, dictionary$actorName)
Events_team4$receiver <- match(Events_team4$receiver, dictionary$actorName)

save(stats, std_stats, Events_team4, file = "team4_estimation_data.RData")

# Prepare data for out-of-sample-prediction
# Subset: second and third day 
an.subset <- min(which(lubridate::day(team4_edgelist$time) %in% c(18, 19))):
	max(which(lubridate::day(team4_edgelist$time) %in% c(18, 19)))

# Compute stats
stats <- remstats(edgelist = reh_team4, tie_effects = effects, 
									attributes = team4_attributes, directed = FALSE,
									start = min(an.subset), stop = max(an.subset))$statistics

# Standardize the statistics (!!! ask whether this is the correct way to do it??)
std_stats <- lapply(1:dim(stats)[3], function(i) {
	(stats[,,i]-mean(stats[,,i]))/sd(stats[,,i])
})
std_stats <- do.call(abind, args = list(std_stats, along = 3))
apply(std_stats, 3, sd) 
dimnames(std_stats) <- dimnames(stats)

# Prepare data.frame for estimation
Events_team4 <- team4_edgelist[an.subset,c(2,3)]
colnames(Events_team4) <- c("sender", "receiver")
Events_team4$index <- reh_team4$edgelist[an.subset,2] + 1
Events_team4$sender <- match(Events_team4$sender, dictionary$actorName)
Events_team4$receiver <- match(Events_team4$receiver, dictionary$actorName)

# Indices
secondDay <- min(which(lubridate::day(team4_edgelist$time) == 18)):
	max(which(lubridate::day(team4_edgelist$time) == 18))
thirdDay <- min(which(lubridate::day(team4_edgelist$time) == 19)):
	max(which(lubridate::day(team4_edgelist$time) == 19))

# Save
save(stats, std_stats, Events_team4, secondDay, thirdDay, file = "team4_out_of_sample_prediction_data.RData")