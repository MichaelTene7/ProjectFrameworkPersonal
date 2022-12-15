library(RERconverge)
source("Src/Reu/ZonomNameConvertMatrixCommon.R")
source("Src/Reu/ZonomNameConvertVector.R")
source("Src/Reu/categorizePaths.R")

allInsectivoryData = read.csv("Data/allInsectivoryCorrelationFile.csv")

head(allInsectivoryData[order(allInsectivoryData$P),])
head(allInsectivoryData[order(allInsectivoryData$p.adj),])
hist(allInsectivoryData$P)
hist(allInsectivoryData$p.adj)
#-----------------------------

RERs = readRDS("Data/allInsectivoryRERFile.rds")
Paths = readRDS("Data/allInsectivoryPathsFile.rds")
plotRers(RERs,"KIAA0825", Paths )

CNRers = ZonomNameConvertMatrixCommon(RERs)
CNNames = CNRers[1,]
names(CNNames)
bats = grep("bat", names(CNNames))

TricolorPaths = Paths
TricolorPaths[intersect(which(Paths==1), bats)] =2
plotRers(CNRers,"KIAA0825", TricolorPaths )


grep("hog", names(CNNames))
TricolorPaths[grep("hog", names(CNNames))]
Paths[grep("hog", names(CNNames))]
Paths[869]


# make edit via tree instead

binaryTree = readRDS("Data/allInsectivoryBinaryForegroundTree.rds")
plotTree(binaryTree)

newTipName = ZonomNameConvertVectorCommon(binaryTree$tip.label)
RenameTree = binaryTree
RenameTree$tip.label = newTipName
dev.new(height = 40, width = 40)

plotTree(RenameTree)

#manualTree = click_select_foreground_branches(RenameTree)
#saveRDS(manualTree, "Results/manualallInsColorTree.rds")
manualTree = readRDS("Results/manualallInsColorTree.rds")
manualZNameTree = manualTree
plotTree(manualTree)
manualZNameTree$tip.label = binaryTree$tip.label
zonomMaster = readRDS("../RunRER/data/RemadeTreesAllZoonomiaSpecies.rds")

manualPaths = tree2Paths(manualZNameTree, zonomMaster)

plotRers(CNRers,"KIAA0825", manualPaths, sortrers = T)

?foreground2Tree()
str(manualTree)
manualTree$edge.length

functionPaths = categorizePaths(binaryTree, zonomMaster, "functionPath", overwrite = T)
quadcolorPaths = categorizePaths(binaryTree, zonomMaster, "quadcolor", overwrite = T)
plotTree(readRDS("Results/quadcolorManualFGTree.rds"))
plotRers(CNRers,"KIAA0825", manualPaths, sortrers = F)
plotRers(CNRers,"KIAA0825", functionPaths, sortrers = T)
plotRers(CNRers,"KIAA0825", Paths)


par(mfrow = 2)
plotRers(CNRers,"ZNF292", manualPaths, sortrers = T)
plotRers(CNRers,"ZNF292", Paths, sortrers = T)
