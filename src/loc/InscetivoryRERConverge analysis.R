correlData = read.csv("Data/InsectivoryCorrelationFile.csv")
head(correlData[order(correlData$p),])
hist(correlData$p.adj)
hist(correlData$P)

foregroundTree = readRDS("Data/insectivoryBinaryForegroundTree.rds")
