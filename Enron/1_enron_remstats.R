library(remstats)
library(remify)
library(dplyr)
library(abind)

# Load the enron data
load('enron.rda')
load('new_enron_messages.RData')

# Look at the times
hist(as.POSIXct(new.enron.messages$time), breaks = "days") 
# Little events happen in the first few years. 
times <- as.POSIXct(new.enron.messages$time)
table(lubridate::year(times))

# Propose to take a new subset:
an.subset <- seq(min(which(lubridate::year(times) == 2001)), 
								 by = 1, length.out = 10000)

# Prepare a data.frame to compute statistics 
edgelist_enron <- data.frame(
	time   = as.double.POSIXlt(new.enron.messages$time)[1:max(an.subset)], 
	actor1 = new.enron.messages$sender[1:max(an.subset)],
	actor2 = new.enron.messages$receiver[1:max(an.subset)])

reh_enron <- reh(edgelist_enron, model = "actor", actors = 1:156) 

enron.attributes <- data.frame(
	id = 1:156, time = 0, 
	department = as.numeric(factor(enron.actors$department)),
	gender = as.numeric(factor(enron.actors$gender)),
	seniority = as.numeric(factor(enron.actors$seniority)),
	fromLegal = as.numeric(enron.actors$department == "Legal"), 
	fromTrading = as.numeric(enron.actors$department == "Trading"), 
	fromOther = as.numeric(enron.actors$department == "Other"),
	isMale = as.numeric(enron.actors$gender == "Male"),
	isFemale = as.numeric(enron.actors$gender == "Female"), 
	isJunior = as.numeric(enron.actors$seniority == "Junior"),
	isSenior = as.numeric(enron.actors$seniority == "Senior"),
	title = as.numeric(factor(enron.actors$title)), 
	division = as.numeric(factor(enron.actors$long.department))
) # note: many unique titles, one unique division

# Compute interval inertia, reciprocity, indegreeReceiver, outdegreeReceiver,
# otp, itp, osp and isp
intervals <- c(1,2,7,14,31,91)*24*60*60
effects <- ~ inertia(scaling = "std") + reciprocity(scaling = "std") +
	indegreeReceiver(scaling = "std") + outdegreeReceiver(scaling = "std") +
	otp(scaling = "std") + itp(scaling = "std") + osp(scaling = "std") + 
	isp(scaling = "std")

wstats <- lapply(seq_along(intervals), function(i) {
	stats <- remstats(reh_enron, receiver_effects = effects, memory = "window", 
										memory_value = intervals[i], start = min(an.subset), 
										stop = max(an.subset))$statistics$receiver_stats
	dimnames(stats)[[3]] <- paste0(dimnames(stats)[[3]], intervals[i]/(60*60))
	stats
}) 

wstats <- do.call(abind, args = list(wstats, along = 3))
str(wstats)
save(wstats, file = "wstats.rda")

# Prepare matrices for tie effects
get.combo <- function(var1, var2, actors) {
	x <- var1 %*% t(var2)
	rownames(x) <- colnames(x) <- actors
	x
}

attach(enron.attributes)
LegalLegal <- get.combo(fromLegal, fromLegal, enron.attributes$id)
TradingLegal <- get.combo(fromTrading, fromLegal, enron.attributes$id)
JuniorLegal <- get.combo(isJunior, fromLegal, enron.attributes$id)
FemaleLegal <- get.combo(isFemale, fromLegal, enron.attributes$id)
LegalTrading <- get.combo(fromLegal, fromTrading, enron.attributes$id)
TradingTrading <- get.combo(fromTrading, fromTrading, enron.attributes$id)
JuniorTrading <- get.combo(isJunior, fromTrading, enron.attributes$id)
FemaleTrading <- get.combo(isFemale, fromTrading, enron.attributes$id)
LegalJunior <- get.combo(fromLegal, isJunior, enron.attributes$id)
TradingJunior <- get.combo(fromTrading, isJunior, enron.attributes$id)
JuniorJunior <- get.combo(isJunior, isJunior, enron.attributes$id)
FemaleJunior <- get.combo(isFemale, isJunior, enron.attributes$id)
LegalFemale <- get.combo(fromLegal, isFemale, enron.attributes$id)
TradingFemale <- get.combo(fromTrading, isFemale, enron.attributes$id)
JuniorFemale <- get.combo(isJunior, isFemale, enron.attributes$id)
FemaleFemale <- get.combo(isFemale, isFemale, enron.attributes$id)
detach(enron.attributes)

# Compute a bunch of exogenous statistics  #!!! scaling of the exo variables?
# Note: division has missing values, therefore I remove it for now
receiver_effects <- ~ 
	same("department") + same("gender") + same("seniority") + same("title") +
	receive("fromLegal") + receive("fromTrading") + receive("fromOther") +
	receive("isMale") + receive("isFemale") +
	receive("isJunior") + receive("isSenior") +
	tie(LegalLegal, "LegalLegal") +
	tie(TradingLegal, "TradingLegal") +
	tie(JuniorLegal, "JuniorLegal") +
	tie(FemaleLegal, "FemaleLegal") +
	tie(LegalTrading, "LegalTrading") +
	tie(TradingTrading, "TradingTrading") +
	tie(JuniorTrading, "JuniorTrading") +
	tie(FemaleTrading, "FemaleTrading") +
	tie(LegalJunior, "LegalJunior") +
	tie(TradingJunior, "TradingJunior") +
	tie(JuniorJunior, "JuniorJunior") +
	tie(FemaleJunior, "FemaleJunior") +
	tie(LegalFemale, "LegalFemale") +
	tie(TradingFemale, "TradingFemale") +
	tie(JuniorFemale, "JuniorFemale") +
	tie(FemaleFemale, "FemaleFemale") 

enron.remstats <- remstats(reh_enron, receiver_effects = receiver_effects,
													 attributes = enron.attributes, start = min(an.subset), stop = max(an.subset)) 

# Combine all the stats
stats <- abind(wstats, enron.remstats$statistics$receiver_stats, along = 3)

# Add interactions
interactions <- lapply(60:75, function(i) {
	x <- lapply(1:48, function(j) {
		stats[,,i]*stats[,,j]
	})
	y <- do.call(abind, args = list(x, along = 3))
	dimnames(y)[[3]] <- paste0(dimnames(stats)[[3]][i], ".x.", dimnames(stats)[[3]][1:48])
	y
})
interactions <- do.call(abind, args = list(interactions, along = 3))
# Interaction before or after standardization??

# Combine all the stats
stats <- abind(stats, interactions, along = 3)
dimnames(stats)

# Standardize the statistics (!!! ask whether this is the correct way to do it??)
std_stats <- lapply(1:dim(stats)[3], function(i) {
	(stats[,,i]-mean(stats[,,i]))/sd(stats[,,i])
})
std_stats <- do.call(abind, args = list(std_stats, along = 3))
apply(std_stats, 3, sd) 
dimnames(std_stats) <- dimnames(stats)

# Prepare data.frame for estimation
dictionary <- attributes(reh_enron)$dictionary$actors
Events_enron <- new.enron.messages[an.subset,c(4,5,5)] 
colnames(Events_enron) <- c("sender", "receiver", "index")
Events_enron$sender <- match(Events_enron$sender, dictionary$actorName)
Events_enron$receiver <- match(Events_enron$receiver, dictionary$actorName)
Events_enron$index <- Events_enron$receiver

save(stats, std_stats, Events_enron, file = "enron_estimation_data.RData")
