.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
m=readTrees("C:/Users/Michael/AppData/Local/R/win-library/4.2/RERconverge/extdata/subsetMammalGeneTrees.txt")
mainTree = m
mainTree$masterTree$tip.label
# -- Import a tree --
mainTree = read.tree("filename.txt")

# -- Import manual annotations CSv -- 

manualAnnots = read.csv("Data/manualAnnotationsSheet1.csv")

# -- start subsetting Annots to the foreground only --

# Version which is prepped for multiple contitions 
subsetManualAnnots = manualAnnots[intersect(which(manualAnnots$Ins_v_herbs %in% c(0,1))),]

#single condition version
subsetManualAnnots = manualAnnots[(manualAnnots$Ins_v_herbs %in% c(0,1)),]

# -- setup foreground --

foregroundManualAnnots = subsetManualAnnots[(subsetManualAnnots$Ins_v_herbs %in% 1),]


foregroundNames = foregroundManualAnnots$Tip_Label..Red.is.included.in.CMU.enhancer.dataset..but.missing.alignment.
foregroundNames = c("Walrus", "Seal", "Killer_whale", "Dolphin", "Manatee")

binaryForegroundTreeOutput = foreground2Tree(foregroundNames, mainTree)
