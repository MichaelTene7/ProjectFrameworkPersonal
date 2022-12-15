library(RERconverge)
source("Src/Reu/ZonomNameConvertMatrixCommon.R")
library(devtools)



allInsectivoryData = read.csv("Data/allInsectivoryCorrelationFile.csv")

head(allInsectivoryData[order(allInsectivoryData$P),])
head(allInsectivoryData[order(allInsectivoryData$p.adj),])
hist(allInsectivoryData$P)
hist(allInsectivoryData$p.adj)
#-----------------------------

RERs = readRDS("Data/allInsectivoryRERFile.rds")
Paths = readRDS("Data/allInsectivoryPathsFile.rds")
plotRers(RERs,"KIAA0825", Paths )
?plotRers()

names = names(Paths)


CNRers = ZonomNameConvertMatrixCommon(RERs)
#pdf("Output/RERImage.pdf", width = 12, height = 20)
plotRers(CNRers,"KIAA0825", Paths, sortrers = T)
dev.off()
?plotRers
str(Paths)

CNNames = CNRers[1,]
names(CNNames)
bats = grep("bat", names(CNNames))
boPath = Paths[bats]
fgbatsPath = boPath[boPath ==1]

coloredPath = Paths
coloredPath[]

Pathsshort = Paths[Paths == 1 %in% bats]
Paths[bats]
length(which(Paths == 1 %in% bats))

TricolorPaths = Paths
TricolorPaths[intersect(which(Paths == 1), bats)] = 2
plotRers(CNRers,"KIAA0825", TricolorPaths)
dev.new(width = 12, height = 40)

function (rermat = NULL, index = NULL, phenv = NULL, rers = NULL, 
          method = "k", xlims = NULL, plot = 1, xextend = 0.2, sortrers = F) 
{
  if (is.null(rers)) {
    e1 = rermat[index, ][!is.na(rermat[index, ])]
    colids = !is.na(rermat[index, ])
    e1plot <- e1
    if (exists("speciesNames")) {
      names(e1plot) <- speciesNames[names(e1), ]
    }
    if (is.numeric(index)) {
      gen = rownames(rermat)[index]
    }
    else {
      gen = index
    }
  }
  else {
    e1plot = rers
    gen = "rates"
  }
  names(e1plot)[is.na(names(e1plot))] = ""
  if (!is.null(phenv)) {
    phenvid = phenv[colids]
    fgdcor = getAllCor(rermat[index, , drop = F], phenv, 
                       method = method)
    plottitle = paste0(gen, ": rho = ", round(fgdcor$Rho, 
                                              4), ", p = ", round(fgdcor$P, 4))
    fgd = setdiff(names(e1plot)[phenvid == 1], "")
    df <- data.frame(species = names(e1plot), rer = e1plot, 
                     stringsAsFactors = FALSE) %>% mutate(mole = as.factor(ifelse(phenvid > 
                                                                                    0, 2, 1)))
  }
  else {
    plottitle = gen
    fgd = NULL
    df <- data.frame(species = names(e1plot), rer = e1plot, 
                     stringsAsFactors = FALSE) %>% mutate(mole = as.factor(ifelse(0, 
                                                                                  2, 1)))
  }
  if (sortrers) {
    df = filter(df, species != "") %>% arrange(desc(rer))
  }
  if (is.null(xlims)) {
    ll = c(min(df$rer) * 1.1, max(df$rer) + xextend)
  }
  else {
    ll = xlims
  }
  g <- ggplot(df, aes(x = rer, y = factor(species, levels = unique(ifelse(rep(sortrers, 
                                                                              nrow(df)), species[order(rer)], sort(unique(species))))), 
                      col = mole, label = species)) + scale_size_manual(values = c(3, 
                                                                                   3)) + geom_point(aes(size = mole)) + scale_color_manual(values = c("deepskyblue3", 
                                                                                                                                                      "brown1")) + scale_x_continuous(limits = ll) + geom_text(hjust = 1, 
                                                                                                                                                                                                               size = 5) + ylab("Branches") + xlab("relative rate") + 
    ggtitle(plottitle) + geom_vline(xintercept = 0, linetype = "dotted") + 
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
          legend.position = "none", panel.background = element_blank(), 
          axis.text = element_text(size = 18, face = "bold", 
                                   colour = "black"), axis.title = element_text(size = 24, 
                                                                                face = "bold"), plot.title = element_text(size = 24, 
                                                                                                                          face = "bold")) + theme(axis.line = element_line(colour = "black", 
                                                                                                                                                                           size = 1)) + theme(axis.line.y = element_blank())
  if (plot) {
    print(g)
  }
  else {
    g
  }
}

