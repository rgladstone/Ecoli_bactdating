#!/usr/bin/env Rscript

library("ape")
library("BactDating")

#4 arguments the prefix, mcmc chain length, sampling, model
args = commandArgs(trailingOnly=TRUE)

#mcmc chain length
nbIts <- as.numeric(args[2])
cluster_tree <-loadGubbins(args[1])
dates <- read.csv(paste(args[1], "_dates.csv", sep=""), header = TRUE, sep =",", stringsAsFactors = FALSE)
rownames(dates) <- dates$ID
cluster_dates <- dates[cluster_tree$tip.label,2]

results_eq_dates <- bactdate(tree=cluster_tree, date=rep(2019,length(cluster_dates)), initMu = NA, initAlpha = NA, initSigma = NA,
                    updateMu = T, updateAlpha = T, updateSigma = T, updateRoot = T,
                    nbIts = as.numeric(args[2]), thin = ceiling(as.numeric(args[3])), useCoalPrior = T,
                    model = args[4], useRec = T, showProgress = F)

save(results_eq_dates, file = paste(args[1], args[2], args[3],args[4],"BD_eq_dates.RData", sep="_"))
