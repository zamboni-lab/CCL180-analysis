findTraitType<-function(Heatamp_clust,TreeEnrichment_Mat){
  library(ComplexHeatmap)
  library(dendextend)

  Heatamp_clust_=column_dend(Heatamp_clust)
  cut_avg <- cutree(Heatamp_clust_, k = 3)
  
  Type1_cells=names(cut_avg[cut_avg==1])
  Type3_cells=names(cut_avg[cut_avg==3])
  Type2_cells=names(cut_avg[cut_avg==2])
  
  TraitToKeep=c()
  # Find which trait
  for(i in 1:dim(TreeEnrichment_Mat)[1]){
    signCellLines=colnames(TreeEnrichment_Mat)[which(abs(TreeEnrichment_Mat[i,])>1)]
    # If all cell lines are significant from each type
    res3=length(which(signCellLines %in% Type3_cells))==length(Type3_cells)
    res2=length(which(signCellLines %in% Type2_cells))==length(Type2_cells)
    res1=length(which(signCellLines %in% Type1_cells))==length(Type1_cells)
    # If at least one 
    sumres=res1+res2+res3
    if(sumres>=1){
      TraitToKeep=c(TraitToKeep,i)
    }
  }

  return(TraitToKeep)
}