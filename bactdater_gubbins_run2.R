#!/usr/bin/env Rscript

library("ape")
library("BactDating")

#two arguments the prefix and the mcmc chain length
args = commandArgs(trailingOnly=TRUE)

#mcmc chain length
nbIts <- as.numeric(args[2])
cluster_tree <-loadGubbins(args[1])
dates <- read.csv(paste(args[1], "_dates.csv", sep=""), header = TRUE, sep =",", stringsAsFactors = FALSE)
rownames(dates) <- dates$ID
cluster_dates <- dates[cluster_tree$tip.label,2]
results <- bactdate(tree=cluster_tree, date=cluster_dates, initMu = NA, initAlpha = NA, initSigma = NA,
                            updateMu = T, updateAlpha = T, updateSigma = T, updateRoot = T,
                            nbIts = as.numeric(args[2]), thin = ceiling(as.numeric(args[3])), useCoalPrior = T,
                            model = args[4], useRec = T, showProgress = F)

jpeg(paste(args[1],"_bactdate_plot_r2.jpg", sep=""))
plot(results,"trace")
dev.off()

save(results, file = paste(args[1], args[2], args[3],args[4],"BD_r2.RData", sep="_"))
