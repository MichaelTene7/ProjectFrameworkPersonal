# --- non permulated allINsectivory (AI) ----
nonPermInsectivoryPVal = readRDS("Data/allInsectivoryCorrelationFile.rds")
sortedNonPermInsectivoryPVal= nonPermInsectivoryPVal[order(nonPermInsectivoryPVal$p.adj),]
topNonPermAIgenes = rownames(sortedNonPermInsectivoryPVal)[1:30]

# --- Permualted old Insectivory (OAI) ---
PermulatedOldAllInsectivoryPVal = readRDS("Data/oldAllInsectivoryCombinedAppendedPermulationsPValue.rds")
sortedPermulatedOldAllInsectivoryPVal = PermulatedOldAllInsectivoryPVal[order(PermulatedOldAllInsectivoryPVal)]
topPermOAIgenes = names(sortedPermulatedOldAllInsectivoryPVal)[1:30]
topPermOAIgenes



# --- Permulated new Insectivory (NAI) ----
PermulatedNewAllInsectivoryPVals = readRDS("Data/allInsectivoryCombinedAppendedPermulationsPValue.rds")
sortedPermulatedNewAllInsectivoryPVals = PermulatedNewAllInsectivoryPVals[order(PermulatedNewAllInsectivoryPVals)]
topPermedNAIgenes = names(sortedPermulatedNewAllInsectivoryPVals)[1:30]
topPermedNAIgenes

#--- Insectivory Overlaps ---

#perm-no-perm OAI Overlap
overlapOAIPnp = topNonPermAIgenes[topNonPermAIgenes %in% topPermOAIgenes]
overlapOAIPnp

#Perm-no-perm NAI Overlap
overlapNAIPnp = topPermOAIgenes[topPermedNAIgenes %in% topNonPermAIgenes]
overlapNAIPnp

#Perm-NAI-Perm-OAI Overlap
overlapNvOPermAI = topPermOAIgenes[topPermOAIgenes %in% topPermedNAIgenes]
overlapNvOPermAI



# -------  Carn V Herbs --------

# ---- non permulated New CarnvHerbs (New CVH)
nonPermNewCVHpVal = readRDS("Data/carnvHerbsCorrelationFile.rds")
sortedNonPermNewCVHpVal = nonPermNewCVHpVal[order(nonPermNewCVHpVal$p.adj),]
topNonPermNewCVHgenes = rownames(sortedNonPermNewCVHpVal)[1:30]
topNonPermNewCVHgenes

# ---- permulated New CarnvHerbs (New CVH)
permNewCVHpVal = readRDS("Data/carnvHerbsCombinedPrunedFastAppendedPermulationsPValue.rds")
sortedPermNewCVHpVal = permNewCVHpVal[order(permNewCVHpVal)]
topPermNewCvHgenes = names(sortedPermNewCVHpVal)[6:36]                          #First six are either 0 in new or NA in old
topPermNewCvHgenes

# ---- non permulated Old CarnvHerbs (Old CVH) 
oldCarnvHerbspVal = read.csv("Data/CorrelationFishCarnHerbLaurasiatheria_wminsp_wnopermp_wcount_waddlperms0907.csv")
sortedNoPermOldCVHpVal = oldCarnvHerbspVal[order(oldCarnvHerbspVal$p.adj),1:10]
topNoPermOldCVHgenes = sortedNoPermOldCVHpVal$X[1:30]
topNoPermOldCVHgenes

# ---- permulated Old CarvHerbs (Old CVH)
oldCarnvHerbspVal = read.csv("Data/CorrelationFishCarnHerbLaurasiatheria_wminsp_wnopermp_wcount_waddlperms0907.csv")
sortedOldCVHpVal= oldCarnvHerbspVal[order(oldCarnvHerbspVal$permp.2sided),1:10]
topPermOldCVHgenes = sortedOldCVHpVal$X[1:30]
topPermOldCVHgenes

# --- CVH Overlaps ---

#Perm no perm new CVH
overlapNewCVHPnp = topNonPermNewCVHgenes[topNonPermNewCVHgenes %in% topPermNewCvHgenes]
overlapNewCVHPnp

#Perm no Perm old CVH
overlapOldCVHPnp = topNoPermOldCVHgenes[topNoPermOldCVHgenes %in% topPermOldCVHgenes]
overlapOldCVHPnp

# New Vs Old Permulated CVH
overlapOldvNewPermCVH = topPermNewCvHgenes[topPermNewCvHgenes %in% topPermOldCVHgenes]
overlapOldvNewPermCVH

# New Vs Old No-permulated CVH
overlapNewvOldNoPermCVH = topNoPermOldCVHgenes[topNoPermOldCVHgenes %in% topNonPermNewCVHgenes]
overlapNewvOldNoPermCVH

# --- Expand top to 100 
topNonPermNewCVHgenes = rownames(sortedNonPermNewCVHpVal)[1:100]
topPermNewCvHgenes = names(sortedPermNewCVHpVal)[6:106]
topNoPermOldCVHgenes = sortedNoPermOldCVHpVal$X[1:100]
topPermOldCVHgenes = sortedOldCVHpVal$X[1:100]

# --- CVH Expanded Overlaps ---

#Perm no perm new CVH
overlapNewCVHPnp = topNonPermNewCVHgenes[topNonPermNewCVHgenes %in% topPermNewCvHgenes]
overlapNewCVHPnp

#Perm no Perm old CVH
overlapOldCVHPnp = topNoPermOldCVHgenes[topNoPermOldCVHgenes %in% topPermOldCVHgenes]
overlapOldCVHPnp

# New Vs Old Permulated CVH
overlapOldvNewPermCVH = topPermNewCvHgenes[topPermNewCvHgenes %in% topPermOldCVHgenes]
overlapOldvNewPermCVH

# New Vs Old No-permulated CVH
overlapNewvOldNoPermCVH = topNoPermOldCVHgenes[topNoPermOldCVHgenes %in% topNonPermNewCVHgenes]
overlapNewvOldNoPermCVH


# merge dataframes 
names(nonPermNewCVHpVal)
reformatNonPermNewCVHpVal = nonPermNewCVHpVal
names(reformatNonPermNewCVHpVal) = c("NewRho", "NewN", "NewP", "Newp.adj")
CVHCombined =  cbind(oldCarnvHerbspVal, reformatNonPermNewCVHpVal)
CVHCombined = cbind(CVHCombined, permNewCVHpVal)
names(CVHCombined)[19] = "NewPermP"

CVHCombinedSortNewPerm = CVHCombined[order(CVHCombined$NewPermP),]
NewPermComparePs = CVHCombinedSortNewPerm[,c(1,5,6,17,18,19)]
CVHCombinedSortOldPerm = CVHCombined[order(CVHCombined$permp.2sided),]
OldPermComparePs = CVHCombinedSortOldPerm[,c(1,5,6,17,18,19)]
CVHCombinedSortOldPerm1Side = CVHCombined[order(CVHCombined$permp.1sided),]
OldPerm1SideComparePs = CVHCombinedSortOldPerm[,c(1,5,6,7,17,18,19)]


APCDD1Location = grep("APCDD1", sortedOldCarnvHerbsPVal$X)

sortedOldCarnvHerbsPVal[11436,]
