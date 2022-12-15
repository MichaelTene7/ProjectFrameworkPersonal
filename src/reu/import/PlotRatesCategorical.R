#Which part of correlateTreesChar3 creates a similar "projection object"?
singleTreeChar3=function(treeListObj, globaltree, char, gene, char.ref=NULL, cutoff=0.000004*3, start=NULL, sqrt=F){
  if(!is.null(char.ref)){
    cm=intersect(cm, names(char.ref))
    char.ref=char.ref[cm]
  }
  cm=intersect(names(char), globaltree$tip)
  master.tree=pruneTree(globaltree, cm)
  if(!is.null(char.ref)){
    refTree=edgeVarsDiff(master.tree, char.ref)
  }
  charTree=edgeVarsDiff(master.tree, char)
	tree1=treeListObj$trees[[gene]]
      both=intersect(tree1$tip.label, cm)
      tree1=unroot(pruneTree(tree1,both))
      charTreeUse=unroot(pruneTree(charTree, both))
      if(!is.null(char.ref)){
        refTreeUse=unroot(pruneTree(refTree,both))
      }
      allreport=treeListObj$report[,both] 
      ss=rowSums(allreport)
      iiboth=which(ss==length(both))
      
      if(length(both)<10){
        next
      }

      torm=setdiff(treeListObj$masterTree$tip.label, both)
       
      ee=edgeIndexRelativeMaster(tree1, treeListObj$masterTree)
      ii= match(namePaths(ee,T), colnames(treeListObj$paths)) 
      
      allbranch=treeListObj$paths[iiboth,ii]
      nv=t(projection(t(allbranch), method="AVE", returnNV = T))
    
      iibad=which(allbranch<cutoff)
      allbranch=scaleMat_c(allbranch)
      if(sqrt){
        nv=sqrt(nv)
        allbranch=sqrt(allbranch)
      }
      allbranch[iibad]=NA
      
      nb=length(both)
      ai=which(maxn[iiboth, 1]==nb)
      if(cutoff>0){
        proj=naresid(allbranch[ai, ,drop=F], as.vector(nv))	### <------ Here is the problem, in the "resid" and "naresid" functions above.
      }
      else{
        proj=resid(allbranch[ai, ,drop=F], cbind(as.vector(nv)))
      }

      message("ai")
      show(length(ai))
      ee=charTreeUse$edge.length
      if(!is.null(char.ref)){
        ee=resid(rbind(ee), mod<-model.matrix(~1+refTreeUse$edge.length))
        proj=naresid(proj, refTreeUse$edge.length) #This is the "projection object," I think.
        
      }
      ee=as.vector(ee)
      
      for (j in 1:length(ai)){
     
        jj=ai[j]
        corout[iiboth[jj],4]=paste(corout[iiboth[jj],4], as.character((i)))
        if (sum(!is.na(proj[j,]))>0.5*ncol(proj)){
          # show(names(trees$trees)[iiboth[j]])
          ii2=which(!is.na(proj[j,]))
          
          cres=cor.test(proj[j,ii2], ee[ii2], method="s")
       
          corout[iiboth[jj],1:3]=c(cres$estimate, cres$statistic, cres$p.value)
          
        }
      }
}
#modified to handle binary trees
correlateTreesProj=function(treeIn1, treeIn2, treeObj, residfun=residLN, plot=F, cutoff=-1, usePaths=T, tree1Bin=F, useIndex=F, species.list=NULL){
  if(is.character(treeIn1)){
    tree1=treeObj$trees[[treeIn1]]
  }
  else{
    tree1=treeIn1
  }
  if(is.character(treeIn2)){
    tree2=treeObj$trees[[treeIn2]]
  }
  else{    
    tree2=treeIn2
  }
  
  both=intersect(tree1$tip.label, tree2$tip.label)
  if(!is.null(species.list)){
    both=intersect(both, species.list)
  }
  if(is.character(treeIn1)){
    bothIndex=which(colSums(treeObj$report[c(treeIn1, treeIn2),])==2)
  }
  
  torm=setdiff(treeObj$masterTree$tip.label, both) 
  tree1=pruneTree(tree1, both)
  tree1=unroot(tree1)
  if(tree1Bin){ #fix any edges that were created through pruning
    tree1$edge.length[tree1$edge.length>1]=1
  }
  tree2=pruneTree(tree2, both)
  tree2=unroot(tree2)
  allreport=treeObj$report[,both]
  ss=rowSums(allreport)
  iiboth=which(ss==length(both))
  if (! usePaths){
    allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
    for ( i in 1:length(iiboth)){            
      tmptree=rescaleTree(drop.tip(treeObj$trees[[iiboth[i]]], torm)) 
      allbranch[i, ]=tmptree$edge.length
    }
    
  }
  else{
    if(! useIndex){
      ee=edgeIndexRelativeMaster(tree1, treeObj$masterTree)
      ii= match(namePaths(ee,T), colnames(treeObj$paths)) 
      allbranch=treeObj$paths[iiboth,ii]
      allbranch=scaleMat_c(allbranch)
      nv=projection(t(allbranch), method="AVE", returnNV = T)
      mastertree=treeObj$master
      mastertree$edge.length=nv  
      res=correlateTrees(tree1, tree2, mastertree, residfun=residfun, plot=plot, cutoff=cutoff, Tree1Bin=tree1Bin)
      
      res$nv=nv 
      res$allbranch=allbranch
    }
    else{
      allbranch=getBranch(treeObj, bothIndex)
      show(rownames(allbranch)[1])
      allbranch=scaleMat_c(allbranch)
      nv=projection(t(allbranch), method="AVE",returnNV = T)
      rr=resid(allbranch, model.matrix(~0+nv))
      rownames(rr)=rownames(allbranch)
      show(dim(rr))
      show(rownames(allbranch)[1])
      plot(rr[name1,], rr[name2,])
      res=list()
    }
  }
  
  return(res)
}

#Modify the following to plot the three categories
#Uses output from correlateTreesProj.... Does this work with categorical data?
plotCategorical=function(corOutput, title="Species Plot", hjust=1.2){
  show(simpleAUCmat(lab=corOutput$l1, value=corOutput$e2))
  usenames=corOutput$names
  
  tmpn=speciesNames[usenames,]
  usenames[!is.na(tmpn)]=tmpn[!is.na(tmpn)]
  
  #usenames[usenames==""]="internal"
  theme_set(theme_bw( ))
  df=list()
  theme_set(theme_bw())
  df=list()
  iis=order(corOutput$e2)
  df$x=corOutput$e2[iis]
  df$y=usenames[iis]
  df$col=as.factor(corOutput$e1[iis])
  if(hjust>0){
    ll=c(min(df$x)-0.12, max(df$x)+0.01)
  }
  else{
    ll=c(min(df$x)-0.01, max(df$x)+0.12)
  }
  ncol=length(unique(corOutput$e1))
  posscol=c("green","blue","red","pink","brown","gray","cyan","black","orange")
  df=as.data.frame(df)
  ggplot(df, aes(x = x, y=y, col=col, label=y)) + scale_size_manual(values=c(rep(2,ncol)))+ geom_point(aes(size=col))+
    scale_color_manual(values = posscol[1:ncol])+
    scale_x_continuous(limits=ll)+
    geom_text(hjust=hjust, size=3)+
    ylab("Branches")+
    xlab("relative rate")+
    ggtitle(title)+
    geom_vline(xintercept=0, linetype="dotted")+
    theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.position="none")
}

#Figure out edge lengths to reassign after correlateTreesProj
getEdgeLengths=function(treeIn, keepTips, binVec) {
	mammaltree <- pruneTree(treeIn$masterTree, as.vector(keepTips))
	fitER <- rerootingMethod(mammaltree, binVec, model = "ER")
	tipvals <- binVec[match(mammaltree$tip.label,names(binVec))]
	nodevals <- as.numeric(colnames(fitER$marginal.anc))[max.col(fitER$marginal.anc)] #take maxlik value
	alllengths <- c(tipvals, nodevals)
	follengths <- alllengths[mammaltree$edge[,2]]
	names(follengths) <- NULL
	return(follengths) 
}

# Unadorned version of plotBinary
plotRates=function(corOutput, title="Species Plot", hjust=1.2){
  usenames=corOutput$names
  tmpn=speciesNames[usenames,]
  usenames[!is.na(tmpn)]=tmpn[!is.na(tmpn)]
  
  theme_set(theme_bw( ))

  df=list()
  df$x=corOutput$e2
  df$y=usenames
  df=as.data.frame(df)

  ggplot(df, aes(x = x, y=y,  label=y)) +  geom_point() +
    geom_text(hjust=hjust, size=5)+
    ylab("Branches")+
    xlab("relative rate")+
    ggtitle(title)+
    geom_vline(xintercept=0, linetype="dotted")+
    theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.position="none")
}