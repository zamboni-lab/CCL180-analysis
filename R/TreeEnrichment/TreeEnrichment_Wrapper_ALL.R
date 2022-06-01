
TreeEnrichment_Wrapper_ALL<-function(tree,traitVal,traitName,cellLines,LOG_10=FALSE,MinSizeCLuster=18,plotting=TRUE){
  library(dendextend)
  pvalSubTree={}
  posSubTree={}
  for(i in 1:length(traitName)){
    if(length(traitName)==1){
      pvalposlist=TreeEnrichment_Stats(tree,traitVal,cellLines,LOG_10=FALSE,MinSizeCLuster)
    }
    else{
      pvalposlist=TreeEnrichment_Stats(tree,traitVal[,i],cellLines,LOG_10=FALSE,MinSizeCLuster)
    }
    pvalSubTree=rbind(pvalSubTree,pvalposlist$PVal)
    posSubTree=cbind(posSubTree,pvalposlist$posOnTree)
  }
  valOut=list(pvalSubTree,posSubTree)
  return(valOut)
}
