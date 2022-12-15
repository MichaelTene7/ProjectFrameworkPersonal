datadir <- "/Volumes/NEWRSD/MammalDietAdaptation/ConvergentEvolutionAnalysis/SavedData/"
histdir <- "/Volumes/NEWRSD/MammalDietAdaptation/ConvergentEvolutionAnalysis/PHistograms/"
maindir <- "/Volumes/NEWRSD/MammalDietAdaptation/ConvergentEvolutionAnalysis/"
listdir <- "/Volumes/NEWRSD/MammalDietAdaptation/ConvergentEvolutionAnalysis/GeneLists/"
scriptdir <- "/Volumes/NEWRSD/MammalDietAdaptation/ConvergentEvolutionAnalysis/scripts/"
library("phytools")
load(paste(datadir,"mammal_62_aa2.RData",sep=""))
load(paste(datadir,"SpeciesNames.RData",sep=""))
source(paste(scriptdir,"ERCfunctions16012WM.R",sep=""))
#Figure out which genes to look at. Then:
load(paste(datadir,"trophtree.RData",sep=""))
load(paste(datadir,"quanttroph.RData",sep=""))
noneuth <- c("monDom5","sarHar1","macEug2","ornAna1")
euthmammals <- trees$masterTree$tip.label[trees$masterTree$tip.label %in% noneuth == FALSE]
source(paste(scriptdir,"PlotRatesCategorical.R",sep=""))
tmp <- correlateTreesProj(trophtree, "PNLIPRP1", trees, plot=F, usePaths = T, tree1Bin = T, residfun = residLN0, species.list = euthmammals)
tmp$e1 <- getEdgeLengths(trees,tmp$names,quanttroph)
pdf(paste("../GeneRERByBranch/NewEdgeRERForPNLIPRP1.pdf",sep=""))
plotCategorical(tmp,title="PNLIPRP1")
dev.off()
moregenes <- c("GALNT10","PON1", "DHRS4","ABCA13","FBXO25","ACADSB","SLC27A2","GCDH","PARP3")
for (i in c(1:length(moregenes))) {
	tmp <- correlateTreesProj(trophtree, moregenes[i], trees, plot=F, usePaths = T, tree1Bin = T, residfun = residLN0, species.list = euthmammals)
	tmp$e1 <- getEdgeLengths(trees,tmp$names,quanttroph)
	pdf(paste("../GeneRERByBranch/NewEdgeRERFor",moregenes[i],".pdf",sep=""))
	plotCategorical(tmp,title=moregenes[i])
	dev.off()
}