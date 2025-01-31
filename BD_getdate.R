library(BactDating)
library(ape)
library(coda)
library(dispRity)

#setwd("")
#replace mytree with file prefix

tmp <- load(file ="mytree_100000000_1000_arc_BD_r1.RData")
results <- get(tmp)
CI <- results$CI
rownames(CI) <- c(results$tree$tip.label, results$tree$node.label)

node_ages <- tree.age(results$tree,order="present", fossil=TRUE)
node_ages$ages <- node_ages$ages+results$rootdate
rownames(node_ages) <- node_ages$elements
dates <- cbind(node_ages[,2],node_ages[,1],CI[,1:2])
colnames(dates) <- c("node","date", "CI lower", "CI upper")
write.csv(dates, file = "mytree_datesCI.csv", row.names = FALSE)

#write out tree
write.tree(results$tree, file =  "mytree_BD.tre")

