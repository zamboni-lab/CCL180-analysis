
TreeEnrichment_Stats<-function(tree,valCellLines,AllCellLines,OutputNameTree,LOG_10=FALSE,MinSizeCLuster=18){
  
  PositionSignificant={}
  differentCells=setdiff(AllCellLines,labels(tree))
  if(length(differentCells)>0){
    warning("These Cell lines are not on the tree")
    warning("differentCells")
  }
  library(ggplot2)
  library(dendextend)
  PVal={}
  diffMean={}
  posOnTree={}
  subtrees={}
  PValAdj={}
  
  if(LOG_10 ==TRUE){
    valCellLines <- log10(valCellLines)
  }
  
  treeorder=order(labels(tree))
  
  for(i in 2:length(labels(tree)))
  {
    cut_avg <- cutree(tree, k = i)
    # Made a lot of change here, could be problematic
   #dend <- as.dendrogram(tree)
    for(j in 1:i)
    { 
      Cells=names(cut_avg)
      orderCells=unique(cut_avg[treeorder])
      Cluster1=Cells[cut_avg==orderCells[j]]
      clusterName=paste(Cluster1,collapse="")
      
      # In order to consider subtree that hevent been considered before
      if(!clusterName %in% subtrees & length(Cluster1)>=MinSizeCLuster) # if not subtree already considered
      {
        subtrees=rbind(subtrees,clusterName)
        posInCluster = match(AllCellLines,Cluster1)
        posNotInCluster = which(is.na(posInCluster))
        posInCluster=which(posInCluster>0)
        InCluster=valCellLines[posInCluster]
        OutCluster=valCellLines[posNotInCluster]
        
        InCluster= na.omit(InCluster)
        OutCluster=na.omit(OutCluster)
        
        if(length(InCluster)>=MinSizeCLuster & length(OutCluster)>=MinSizeCLuster){
          p=t.test(InCluster,OutCluster,var.equal = FALSE) # Since Sample size so different, use unequal variance
          #p=wilcox.test(InCluster,OutCluster)
          difval=mean(InCluster,na.rm = T)-mean(OutCluster,na.rm = T)
          diffMean=cbind(diffMean,difval)
          # To keep sign
          if(mean(InCluster,na.rm = T)>mean(OutCluster,na.rm = T)){
            PVal=cbind(PVal,p$p.value) 
          }
          else
          {
            PVal=cbind(PVal,-(p$p.value))
          }
          posOnTree=rbind(posOnTree,c(i,j))
        }
        else{
          PVal=cbind(PVal,NA)
          posOnTree=rbind(posOnTree,c(i,j))
        }

        if(FALSE){
          boxplot(InCluster,OutCluster,names=c('In Subtree','Not In Subtree'))
        }
        
        pos=cut_avg[order.dendrogram(tree)]
        findPos=unique(pos)
      
      }
      
    }
    cut_avg_previous=cut_avg
  }
  
  PVal
  posOnTree
  pvalposlist <- list("PVal" = PVal, "posOnTree" = posOnTree)
  
  return(pvalposlist)  
}

