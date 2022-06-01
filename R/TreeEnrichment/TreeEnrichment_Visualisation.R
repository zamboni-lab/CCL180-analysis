TreeEnrichment_Visualisation<-function(tree,pvals,subtree,traitNames,plotting=TRUE){
  library("gplots") 
  library(ggplot2)
  library(dendextend)
  traitToKeep={}
  finalmat={}
  size=dim(pvals)

    # Loop over all traits
    for(i in 1:size[2]){
      # Create empty vector
      val=replicate(length(tree$labels),NA)
      pvalselect=pvals[,i]
      pvalorder=sort(abs(pvalselect), na.last=NA,decreasing = TRUE,index.return=TRUE)
      pvalorder$ix
      position=c((2*i)-1,2*i)
      
      # If there is no pvalue (empty value)
      if(length(pvalorder$ix)!=0){
        traitToKeep=cbind(traitToKeep,i)
        # Loop on all subtree
        for(j in 1:length(pvalorder$ix)){
          possub=subtree[pvalorder$ix[j],position]
          
          cut_avg <- cutree(tree, k = possub[1])
          orderCells=unique(cut_avg[tree$order])
          coloursub=Cells[cut_avg==orderCells[possub[2]]]
          posHit=match(coloursub,labels(tree)) 
          val[posHit]=pvalselect[pvalorder$ix[j]]
        } 
      }
      # To -log10(val) and keeping the sign
      val2=(as.integer(val>=0)*2-1)*-log10(abs(val))
      finalmat=rbind(finalmat,val2)
    }
    rownames(finalmat)=traitNames
    colnames(finalmat)=labels(tree)
    if(plotting){
        if(size[2]>1){
          # traitToKeep only show the
          heatmap.2(finalmat[traitToKeep,],dendrogram='none', Rowv=TRUE, Colv=FALSE,trace='none',col=bluered,key=T,keysize = 1, density.info="none",cexRow=0.8,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray")
        }
        else{
          finalmat=rbind(finalmat,finalmat)
         heatmap.2(finalmat[traitToKeep,],dendrogram='none', Rowv=TRUE, Colv=FALSE,trace='none',col=bluered,key=T,keysize = 1, density.info="none",cexRow=0.8,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray")
        }
    }

 return(finalmat)
}

