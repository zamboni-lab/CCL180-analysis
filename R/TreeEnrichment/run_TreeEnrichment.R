
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute association of omics on metabolic clusters (Tree Enrichment) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd('MetabCCLs/R/TreeEnrichment')

library("R.matlab")
library("factoextra")
library("fpc")
library("NbClust")
library("gplots")
library("ggplot2")
library("pvclust")
library("ComplexHeatmap")
library("R.matlab")
library("circlize")
library("RColorBrewer")

source('TreeEnrichment_Wrapper_ALL.R')
source('TreeEnrichment_Visualisation.R')
source('TreeEnrichment_Stats.R')
source('BasicFunctions.R')


# Load Pathway Score
PCAHeatMap=read.csv("../../data/PathwayScore_180CCL.csv", header = TRUE,row.names = 'X')

# Move to other folder
###### Figure 2 ###### 
Batch1=fixlc(CellLinesMeta$V2)
ha = HeatmapAnnotation(Tissue = CellLinesMeta$V4,   Batch=Batch1, col = list(Batch= c("1" = "#BB7E98", "2" = "pink3", "3" = "pink2","4" = "pink","5" = "steelblue2","6" = "steelblue3","7" = "steelblue4")))
Heatamp_clust=Heatmap(t(PCAHeatMap2), name = "mat",na_col = "grey",col = colorRamp2(c(-1,-0.5, 0,0.5, 1), c("slateblue3","slateblue1", "white", 'springgreen',"springgreen3")),clustering_method_rows='ward.D',clustering_method_columns='ward.D2', row_names_gp = gpar(fontsize = 8),    column_names_gp = gpar(fontsize = 5),top_annotation = ha)
Heatamp_clust

###### Metadata Enrichment ###### 
# Load Cell line Metadata 
CellLinesMeta=read.csv("../../data/CellLinesMeta.csv", header = FALSE)

## Experiment METADATA
CellLinesMeta$V1=fixlc(CellLinesMeta$V1)
CellLinesMeta$V1[CellLinesMeta$V1=='786O']="7860"
CellLinesMeta$V1[CellLinesMeta$V1=='RPMI8826']='RPMI8226'
CellLinesMeta$V1[CellLinesMeta$V1=='NCIH322']='NCIH322M'
setdiff(as.character(CellLinesMeta$V1),rownames(PCAHeatMap))

# Confluency
# First tried algorithm for confluency values - Suspension cells run separately
setdiff(as.character(CellLinesMeta$V1),labels(Heatamp_clust_))
removingSuspensionCells=which(CellLinesMeta$V3<100)
Confluency=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,CellLinesMeta$V3[removingSuspensionCells],'Confluency',t(CellLinesMeta$V1[removingSuspensionCells]),LOG_10=FALSE,MinSizeCLuster=18)
pvalposlist=TreeEnrichment_Stats(Heatamp_clust_,CellLinesMeta$V3[removingSuspensionCells],t(CellLinesMeta$V1[removingSuspensionCells]),LOG_10=FALSE,MinSizeCLuster=18)
finalmatDrug2_corr=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,Confluency[[1]],Confluency[[2]],'Confluency',TRUE)

# Suspension Cells
ConfluencySusCells=TreeEnrichment_Qual_Stat(Heatamp_clust_,CellLinesMeta$V1[CellLinesMeta$V3>100],CellLinesMeta$V1)

# Batch
batch=unique(CellLinesMeta$V2)
pvalSubTree={}
posSubTree={}
for(i in 1:length(batch)){
  BatchTmp=TreeEnrichment_Qual_Stat(Heatamp_clust_,CellLinesMeta$V1[CellLinesMeta$V2==batch[i]],CellLinesMeta$V1)
  pvalSubTree=rbind(pvalSubTree,BatchTmp$PVal)
  posSubTree=cbind(posSubTree,BatchTmp$posOnTree)
}
BatchEnrich=list("PVal" = pvalSubTree, "posOnTree"=posSubTree)
Batch=sprintf("Batch %d",batch)

# Tissue
tissue=unique(CellLinesMeta$V4)
pvalSubTree={}
posSubTree={}
for(i in 1:length(tissue)){
  TissueTmp=TreeEnrichment_Qual_Stat(Heatamp_clust_,CellLinesMeta$V1[CellLinesMeta$V4==tissue[i]],CellLinesMeta$V1)
  pvalSubTree=rbind(pvalSubTree,TissueTmp$PVal)
  posSubTree=cbind(posSubTree,TissueTmp$posOnTree)
}
TissueEnrich=list("PVal" = pvalSubTree, "posOnTree"=posSubTree)

## Metadata from CCLE
CellLinesAnnotation=read.csv("../../data/CellLinesMeta_CCLE.csv", header = TRUE,stringsAsFactors = FALSE)
CellLinesAnnotation$NameAstra=fixlc(CellLinesAnnotation$NameAstra)
length(CellLinesAnnotation$NameAstra[!(CellLinesAnnotation$NameAstra=='')])
CellLinesAnnotation$NameAstra[is.na(CellLinesAnnotation$NameAstra)]=''

cellLineCLLE=unique(CellLinesAnnotation$NameAstra)
cellLineCLLE=cellLineCLLE[cellLineCLLE != ""] 
cellLineCLLE=fixlc(cellLineCLLE)
setdiff(cellLineCLLE,rownames(PCAHeatMap))

# Pathology
pathology=unique(CellLinesAnnotation$Pathology)
pathology<-pathology[!is.na(pathology)]
pathology<-pathology[!pathology=='<undefined>']
pvalSubTree={}
posSubTree={}
for(i in 1:length(pathology)){
  cellsSelected=CellLinesAnnotation$NameAstra[CellLinesAnnotation$Pathology==pathology[i]]
  cellsSelected=na.omit(cellsSelected)
  cellsSelected=fixlc(cellsSelected)
  cellsSelected=cellsSelected[cellsSelected != ""]
  cellsSelected=fixlc(cellsSelected)
  Pathotmp=TreeEnrichment_Qual_Stat(Heatamp_clust_,cellsSelected,cellLineCLLE)
  pvalSubTree=rbind(pvalSubTree,Pathotmp$PVal)
  posSubTree=cbind(posSubTree,Pathotmp$posOnTree)
}
PathologyEnrich=list("PVal" = pvalSubTree, "posOnTree"=posSubTree)

# Histology
histology=unique(CellLinesAnnotation$Histology)
histology<-histology[!is.na(histology)]
histology<-histology[!histology=='<undefined>']
pvalSubTree={}
posSubTree={}
for(i in 1:length(histology)){
  cellsSelected=CellLinesAnnotation$NameAstra[CellLinesAnnotation$Histology==histology[i]]
  cellsSelected=na.omit(cellsSelected)
  cellsSelected=fixlc(cellsSelected)
  cellsSelected=cellsSelected[cellsSelected != ""]
  cellsSelected=fixlc(cellsSelected)
  histotmp=TreeEnrichment_Qual_Stat(Heatamp_clust_,cellsSelected,cellLineCLLE)
  pvalSubTree=rbind(pvalSubTree,histotmp$PVal)
  posSubTree=cbind(posSubTree,histotmp$posOnTree)
}
HistologyEnrich=list("PVal" = pvalSubTree, "posOnTree"=posSubTree)

# Histology Subtype
histologysubtype=unique(CellLinesAnnotation$Hist_Subtype1)
histologysubtype<-histologysubtype[!is.na(histologysubtype)]
histologysubtype<-histologysubtype[!histologysubtype=='<undefined>']
pvalSubTree={}
posSubTree={}
for(i in 1:length(histologysubtype)){
  cellsSelected=CellLinesAnnotation$NameAstra[CellLinesAnnotation$Hist_Subtype1==histologysubtype[i]]
  cellsSelected=na.omit(cellsSelected)
  cellsSelected=fixlc(cellsSelected)
  cellsSelected=cellsSelected[cellsSelected != ""]
  cellsSelected=fixlc(cellsSelected)
  histosubtmp=TreeEnrichment_Qual_Stat(Heatamp_clust_,cellsSelected,cellLineCLLE)
  pvalSubTree=rbind(pvalSubTree,histosubtmp$PVal)
  posSubTree=cbind(posSubTree,histosubtmp$posOnTree)
}
HistologySubtypeEnrich=list("PVal" = pvalSubTree, "posOnTree"=posSubTree)

# Cancer Type
type1=unique(CellLinesAnnotation$type)
type1<-type1[!is.na(type1)]
pvalSubTree={}
posSubTree={}
for(i in 1:length(type1)){
  cellsSelected=CellLinesAnnotation$NameAstra[CellLinesAnnotation$type==type1[i]]
  cellsSelected=na.omit(cellsSelected)
  cellsSelected=fixlc(cellsSelected)
  cellsSelected=cellsSelected[cellsSelected != ""]
  cellsSelected=fixlc(cellsSelected)
  typetmp=TreeEnrichment_Qual_Stat(Heatamp_clust_,cellsSelected,cellLineCLLE)
  pvalSubTree=rbind(pvalSubTree,typetmp$PVal)
  posSubTree=cbind(posSubTree,typetmp$posOnTree)
}
TypeEnrich=list("PVal" = pvalSubTree, "posOnTree"=posSubTree)

# Ethnicity
ethnicity=unique(CellLinesAnnotation$inferred_ethnicity)
ethnicity<-ethnicity[!is.na(ethnicity)]
pvalSubTree={}
posSubTree={}
for(i in 1:length(ethnicity)){
  cellsSelected=CellLinesAnnotation$NameAstra[CellLinesAnnotation$inferred_ethnicity==ethnicity[i]]
  cellsSelected=na.omit(cellsSelected)
  cellsSelected=fixlc(cellsSelected)
  cellsSelected=cellsSelected[cellsSelected != ""]
  cellsSelected=fixlc(cellsSelected)
  ethnmp=TreeEnrichment_Qual_Stat(Heatamp_clust_,cellsSelected,cellLineCLLE)
  pvalSubTree=rbind(pvalSubTree,ethnmp$PVal)
  posSubTree=cbind(posSubTree,ethnmp$posOnTree)
}
EthnicityEnrich=list("PVal" = pvalSubTree, "posOnTree"=posSubTree)

# Mutation rate 
pos=which(CellLinesAnnotation$NameAstra != "")
cellsSelected=fixlc(CellLinesAnnotation$NameAstra[pos])
MutationRateEnrich=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,CellLinesAnnotation$mutRate[pos],'MutationRate',cellsSelected,LOG_10=FALSE,MinSizeCLuster=18)


###### Genomics Enrichment ###### 
## To load data from Li et al. 2019: https://www.nature.com/articles/s41591-019-0404-8 from William R. Sellers's lab

# Mutation data
load('../../data/mut_features.Rdata')
rownames(mut_features)
CellLines_Mut=sub("_.*", "", rownames(mut_features))
rownames(mut_features)=CellLines_Mut

CellLines_Mut_InAstra=match(rownames(PCAHeatMap),CellLines_Mut)
CellLines_Mut_InAstra=fixln(na.omit(CellLines_Mut_InAstra))

mut_features_Astra=mut_features[CellLines_Mut_InAstra,]
pvalSubTree={}
posSubTree={}
for(i in 1:dim(mut_features)[2]){
  posMut=which(mut_features[CellLines_Mut_InAstra,i]==1)
  Muttmp=TreeEnrichment_Qual_Stat(Heatamp_clust_,names(posMut),CellLines_Mut[CellLines_Mut_InAstra])
  pvalSubTree=rbind(pvalSubTree,Muttmp$PVal)
  posSubTree=cbind(posSubTree,Muttmp$posOnTree)
}
MutationEnrich=list("PVal" = pvalSubTree, "posOnTree"=posSubTree)
MutationNames=colnames(mut_features)

# Methylation data
load('../../data/meth_features.Rdata')
rownames(meth_features)
CellLines_Met=sub("_.*", "", rownames(meth_features))
rownames(meth_features)=CellLines_Met
CellLines_Met_InAstra=match(rownames(PCAHeatMap),CellLines_Met)
CellLines_Met_InAstra=fixln(na.omit(CellLines_Met_InAstra))
MethylationEnrich=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,meth_features[CellLines_Met_InAstra,],colnames(meth_features),CellLines_Met[CellLines_Met_InAstra],LOG_10=FALSE,MinSizeCLuster=18)

# CNV data
load('../../data/CN_features.Rdata')
rownames(CN_features)
CellLines_CN=sub("_.*", "", rownames(CN_features))
rownames(CN_features)=CellLines_CN

CellLines_CN_InAstra=match(rownames(PCAHeatMap),CellLines_CN)
CellLines_CN_InAstra=fixln(na.omit(CellLines_CN_InAstra))

CNEnrich=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,CN_features[CellLines_CN_InAstra,],colnames(CN_features),CellLines_CN[CellLines_CN_InAstra],LOG_10=FALSE,MinSizeCLuster=18)

###### Transcriptomics Enrichment ######
# Data taken from CCLE

transcriptomicsAstra=readMat('../../data/transcriptomicsAstra.mat')
transcriptomicsAstra=transcriptomicsAstra$transcriptomicsAstra
GeneSelectedTrans=readMat('../../data/GeneSelectedTrans.mat')
GeneSelectedTransEnsembl <- sapply(GeneSelectedTransEnsembl[[1]], function(n) n[[1]])
CellLineSelectedTrans=readMat('../../data/CellLineSelectedTrans.mat')
CellLineSelectedTrans <- sapply(CellLineSelectedTrans[[1]], function(n) n[[1]])

TranscriptsEnrich=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,transcriptomicsAstra,GeneSelectedTrans,CellLineSelectedTrans,LOG_10=FALSE,MinSizeCLuster=18)

###### Transcription Factors Enrichment ######
# Data taken from Saez Rodriguez lab

TF_activity=readMat('../../data/TF_activity_Astra.mat')
TF_activity_Astra=TF_activity$TF.activity.Astra
TF_list=readMat('../../data/TF_list.mat')
TF_list <- sapply(TF_list[[1]], function(n) n[[1]])
TF_cellline=readMat('../../data/TF_cellline.mat')
TF_cellline <- sapply(TF_cellline[[1]], function(n) n[[1]])

TFEnrich=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,t(TF_activity_Astra),TF_list,TF_cellline,LOG_10=FALSE,MinSizeCLuster=18)

###### EMT Enrichment ######
# Data taken from Rajapakse et al. 2018

EMT_TransitionAstra=readMat('../../data/EMT_TransitionAstra.mat')
EMT_TransitionAstra=t(EMT_TransitionAstra$EMT.TransitionAstra)
EMT_cellline=readMat('../../data/EMT_cellline.mat')
EMT_cellline <- sapply(EMT_cellline[[1]], function(n) n[[1]])

EMTEnrich=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,EMT_TransitionAstra,'EMT',EMT_cellline,LOG_10=FALSE,MinSizeCLuster=18)

###### Proteomics RRPA ######

Protein_Levels_Astra=readMat('../../data/Protein_Levels_Astra.mat')
Protein_Levels_Astra=Protein_Levels_Astra$Protein.Levels.Astra
Protein_list=readMat('../../data/Protein_list.mat')
Protein_list <- sapply(Protein_list[[1]], function(n) n[[1]])
Protein_cells_Astra=readMat('../../data/Protein_cells_Astra.mat')
Protein_cells_Astra <- sapply(Protein_cells_Astra[[1]], function(n) n[[1]])

ProteinEnrich=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,t(Protein_Levels_Astra),Protein_list,Protein_cells_Astra,LOG_10=FALSE,MinSizeCLuster=18)

###### Progeny - Signaling Activity ######
# Data taken from the Saez Rodriguez

Prog_Activity_Astra=readMat('../../data/Prog_Activity_Astra.mat')
Prog_Activity_Astra=Prog_Activity_Astra$Prog.Activity.Astra
Prog_Pathways=readMat('../../data/Prog_Pathways.mat')
Prog_Pathways <- sapply(Prog_Pathways[[1]], function(n) n[[1]])
Proj_cellline=readMat('../../data/Proj_cellline.mat')
Proj_cellline <- sapply(Proj_cellline[[1]], function(n) n[[1]])

ProgEnrich=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,t(Prog_Activity_Astra),Prog_Pathways,Proj_cellline,LOG_10=FALSE,MinSizeCLuster=18)

###### Growth Rate  ######
# Data taken from CCLE achilles (Achilles_cell_line_growth_rate.csv)

GrowthRate_Astra=readMat('../../data/GrowthRate_Levels_Astra.mat')
GrowthRate_Astra <- sapply(GrowthRate_Astra[[1]], function(n) n[[1]])
GrowthRate_cellline=readMat('../../data/GrowthRate_cells_Astra.mat')
GrowthRate_cellline <- sapply(GrowthRate_cellline[[1]], function(n) n[[1]])

GREnrich=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,t(GrowthRate_Astra),'GrowthRate',GrowthRate_cellline,LOG_10=FALSE,MinSizeCLuster=18)


###### Proteomics - Mass Spec ######
# Data Taken from Nusinow et al. 2020
 
Protein_list_MS=readMat('../../data/Protein_list_MS.mat')
Protein_list_MS <- sapply(Protein_list_MS[[1]], function(n) n[[1]])
# Fix format and missing protein names
proteins_to_remove=c()
for(i in 1:length(Protein_list_MS)){
  if(length(Protein_list_MS[[i]])==0){
    proteins_to_remove=c(proteins_to_remove,i)
  }
}
Protein_list_MS2=Protein_list_MS[-proteins_to_remove]
Protein_list_MS2=unlist(Protein_list_MS2)

Protein_Levels_MS_Astra=readMat('../../data/Protein_Levels_MS_Astra.mat')
Protein_Levels_MS_Astra=Protein_Levels_MS_Astra$Protein.Levels.MS.Astra
Protein_Levels_MS_Astra2=t(Protein_Levels_MS_Astra[,-proteins_to_remove])

Protein_cells_MS_Astra=readMat('../../data/Protein_cells_MS_Astra.mat')
Protein_cells_MS_Astra <- sapply(Protein_cells_MS_Astra[[1]], function(n) n[[1]])

ProtMSEnrich=TreeEnrichment_Wrapper_ALL(Heatamp_clust_,t(Protein_Levels_MS_Astra2),Protein_list_MS2,Protein_cells_MS_Astra,LOG_10=FALSE,MinSizeCLuster=18)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Correction for multiple hypothesis testing - merging all omics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Merge Pval - do Benjamini-Hochberg
pvalMerged=c(ConfluencySusCells[[1]],Confluency[[1]],BatchEnrich[[1]],TissueEnrich[[1]],PathologyEnrich[[1]],HistologyEnrich[[1]],HistologySubtypeEnrich[[1]],TypeEnrich[[1]],EthnicityEnrich[[1]],MutationRateEnrich[[1]],MutationEnrich[[1]],MethylationEnrich[[1]],CNEnrich[[1]],TFEnrich[[1]],EMTEnrich[[1]],ProteinEnrich[[1]],ProgEnrich[[1]],TranscriptsEnrich[[1]],ProtMSEnrich[[1]],GREnrich[[1]])
posNA=is.na(pvalMerged)
pvalMerged_noNA=pvalMerged[!posNA]
pvalMerged_noNA_cor=p.adjust(abs(pvalMerged_noNA), method = "BH")
pvalMerged_noNA_cor=(as.integer(pvalMerged_noNA>=0)*2-1)*pvalMerged_noNA_cor
pvalMerged_corr_all=pvalMerged
pvalMerged_corr_all[!posNA]=pvalMerged_noNA_cor
length(pvalMerged_corr_all)

colourPanel=brewer.pal(n = 10, name = "PiYG")
colourPanel[2]='#D8147A'
colourPanel[3]='#FD7A86'
colourPanel[4]='#FFABAB'
colourPanel[6]='white'
colourPanel[5]='white'
colourPanel[7]='#ABFFAB'
colourPanel[8]='#79F979'
colourPanel[9]='#19BB19'

corpval_ConfluencySusCells=matrix(pvalMerged_corr_all[1:length(ConfluencySusCells[[1]])],dim(ConfluencySusCells[[1]]))
TreeEnrichment_Sus=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_ConfluencySusCells,ConfluencySusCells[[2]],'Suspension',TRUE)
PosLast=length(ConfluencySusCells[[1]])

corpval_Confluency=matrix(pvalMerged_corr_all[PosLast+1:length(Confluency[[1]])],dim(Confluency[[1]]))
TreeEnrichment_Conf=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_Confluency,Confluency[[2]],'Confluency',TRUE)
PosLast=PosLast+length(Confluency[[1]])

corpval_BatchEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(BatchEnrich[[1]])],dim(BatchEnrich[[1]]))
TreeEnrichment_BatchEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_BatchEnrich,BatchEnrich[[2]],Batch,TRUE)
PosLast=PosLast+length(BatchEnrich[[1]])

corpval_TissueEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(TissueEnrich[[1]])],dim(TissueEnrich[[1]]))
TreeEnrichment_TissueEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_TissueEnrich,TissueEnrich[[2]],tissue,TRUE)
PosLast=PosLast+length(TissueEnrich[[1]])

corpval_PathologyEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(PathologyEnrich[[1]])],dim(PathologyEnrich[[1]]))
TreeEnrichment_PathologyEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_PathologyEnrich,PathologyEnrich[[2]],pathology,TRUE)
PosLast=PosLast+length(PathologyEnrich[[1]])

corpval_HistologyEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(HistologyEnrich[[1]])],dim(HistologyEnrich[[1]]))
TreeEnrichment_Histology=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_HistologyEnrich,HistologyEnrich[[2]],histology,TRUE)
PosLast=PosLast+length(HistologyEnrich[[1]])

corpval_HistologySubtypeEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(HistologySubtypeEnrich[[1]])],dim(HistologySubtypeEnrich[[1]]))
TreeEnrichment_HistologySubtypeEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_HistologySubtypeEnrich,HistologySubtypeEnrich[[2]],histologysubtype,TRUE)
PosLast=PosLast+length(HistologySubtypeEnrich[[1]])

corpval_TypeEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(TypeEnrich[[1]])],dim(TypeEnrich[[1]]))
TreeEnrichment_TypeEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_TypeEnrich,TypeEnrich[[2]],type1,TRUE)
PosLast=PosLast+length(TypeEnrich[[1]])

corpval_EthnicityEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(EthnicityEnrich[[1]])],dim(EthnicityEnrich[[1]]))
TreeEnrichment_EthnicityEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_EthnicityEnrich,EthnicityEnrich[[2]],ethnicity,TRUE)
PosLast=PosLast+length(EthnicityEnrich[[1]])

corpval_MutationRateEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(MutationRateEnrich[[1]])],dim(MutationRateEnrich[[1]]))
TreeEnrichment_MutationRateEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_MutationRateEnrich,MutationRateEnrich[[2]],'Mutation Rate',TRUE)
PosLast=PosLast+length(MutationRateEnrich[[1]])

corpval_MutationEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(MutationEnrich[[1]])],dim(MutationEnrich[[1]]))
TreeEnrichment_MutationEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_MutationEnrich,MutationEnrich[[2]],colnames(mut_features),FALSE)

# To get only significant
signMut=which(abs(TreeEnrichment_MutationEnrich)>1, arr.ind = TRUE)
signMut2=unique(signMut[,1])

# Because only one gene
matMut=rbind(TreeEnrichment_MutationEnrich[signMut2,],TreeEnrichment_MutationEnrich[signMut2,])
rownames(matMut)=rownames(signMut)[1:2]
heatmap.2(matMut,dendrogram='none', Rowv=TRUE, Colv=FALSE,trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.8,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray", margins=c(5,8))
PosLast=PosLast+length(MutationEnrich[[1]])

corpval_MethylationEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(MethylationEnrich[[1]])],dim(MethylationEnrich[[1]]))
TreeEnrichment_MethylationEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_MethylationEnrich,MethylationEnrich[[2]],colnames(meth_features),FALSE)

# To get only significant
signMet=which(abs(TreeEnrichment_MethylationEnrich)>1, arr.ind = TRUE)
signMet2=unique(signMet[,1])
heatmap.2(TreeEnrichment_MethylationEnrich[signMet2,],dendrogram='none', Rowv=TRUE, Colv=FALSE,hclustfun = function(x) hclust(x, method="ward.D2"),trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.2,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(5,8))
write.csv(rownames(TreeEnrichment_MethylationEnrich[signMet2,]), file = "Methyl.csv")

PosLast=PosLast+length(MethylationEnrich[[1]])

corpval_CNEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(CNEnrich[[1]])],dim(CNEnrich[[1]]))
TreeEnrichment_CNEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_CNEnrich,CNEnrich[[2]],colnames(CN_features),FALSE)
# To get only significant
signCN=which(abs(TreeEnrichment_CNEnrich)>1, arr.ind = TRUE)
signCN2=unique(signCN[,1])
heatmap.2(TreeEnrichment_CNEnrich[c(signCN2,signCN2),],dendrogram='none', Rowv=TRUE, Colv=FALSE,trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.8,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(5,8))
PosLast=PosLast+length(CNEnrich[[1]])

corpval_TFEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(TFEnrich[[1]])],dim(TFEnrich[[1]]))
TreeEnrichment_TFEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_TFEnrich,TFEnrich[[2]],TF_list,FALSE)
# To get only significant
signTF=which(abs(TreeEnrichment_TFEnrich)>1, arr.ind = TRUE)
signTF2=unique(signTF[,1])
heatmap.2(TreeEnrichment_TFEnrich[signTF2,],dendrogram='none', Rowv=TRUE, Colv=FALSE,hclustfun = function(x) hclust(x, method="average"),trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.3,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(5,8))
write.csv(rownames(TreeEnrichment_TFEnrich[signTF2,]), file = "TF.csv")
PosLast=PosLast+length(TFEnrich[[1]])

corpval_EMTEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(EMTEnrich[[1]])],dim(EMTEnrich[[1]]))
TreeEnrichment_EMTEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_EMTEnrich,EMTEnrich[[2]],'EMT',TRUE)
PosLast=PosLast+length(EMTEnrich[[1]])

corpval_ProteinEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(ProteinEnrich[[1]])],dim(ProteinEnrich[[1]]))
TreeEnrichment_ProteinEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_ProteinEnrich,ProteinEnrich[[2]],Protein_list,FALSE)
# To get only significant
signProt=which(abs(TreeEnrichment_ProteinEnrich)>1, arr.ind = TRUE)
signProt2=unique(signProt[,1])
heatmap.2(TreeEnrichment_ProteinEnrich[signProt2,],dendrogram='none', Rowv=TRUE, Colv=FALSE,trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.8,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(5,8))
PosLast=PosLast+length(ProteinEnrich[[1]])

corpval_ProgEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(ProgEnrich[[1]])],dim(ProgEnrich[[1]]))
TreeEnrichment_ProgEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_ProgEnrich,ProgEnrich[[2]],Prog_Pathways,TRUE)
PosLast=PosLast+length(ProgEnrich[[1]])

corpval_TranscriptsEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(TranscriptsEnrich[[1]])],dim(TranscriptsEnrich[[1]]))
TreeEnrichment_TranscriptsEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_TranscriptsEnrich,TranscriptsEnrich[[2]],GeneSelectedTrans,FALSE)
# To get only significant
signTranscripts=which(abs(TreeEnrichment_TranscriptsEnrich)>2, arr.ind = TRUE)
signTranscripts2=unique(signTranscripts[,1])
heatmap.2(TreeEnrichment_TranscriptsEnrich[signTranscripts2,],dendrogram='none', Rowv=TRUE, Colv=FALSE,trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.1,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(5,8))
#write.csv(rownames(TreeEnrichment_TranscriptsEnrich[signTranscripts2,]), file = "Transcripts.csv")
PosLast=PosLast+length(TranscriptsEnrich[[1]])

## Add New proteomics + Growth rate
corpval_ProteinMSEnrich=matrix(pvalMerged_corr_all[PosLast+1:length(ProtMSEnrich[[1]])],dim(ProtMSEnrich[[1]]))
TreeEnrichment_ProteinMSEnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_ProteinMSEnrich,ProtMSEnrich[[2]],Protein_list_MS2,FALSE)
# To get only significant
signProteinMS=which(abs(TreeEnrichment_ProteinMSEnrich)>1, arr.ind = TRUE)
signProteinMS2=unique(signProteinMS[,1])
heatmap.2(TreeEnrichment_ProteinMSEnrich[signProteinMS2,],dendrogram='none', Rowv=TRUE, Colv=FALSE,trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.1,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(5,8))
PosLast=PosLast+length(ProtMSEnrich[[1]])

# Growth Rate GREnrich
corpval_GREnrich=matrix(pvalMerged_corr_all[PosLast+1:length(GREnrich[[1]])],dim(GREnrich[[1]]))
TreeEnrichment_GREnrich=TreeEnrichment_Visualisation_ALL(Heatamp_clust_,corpval_GREnrich,GREnrich[[2]],'Growth Rate',TRUE)
PosLast=PosLast+length(GREnrich[[1]])

# Are you at the end?
PosLast==length(pvalMerged_corr_all)

save.image(file='New_TreeEnrichment_PathwayScore_WithProtMS_GR.RData')

AllTrait=c('Suspension','Confluence',Batch,tissue,pathology,histology,histologysubtype,type1,ethnicity,'Mutation Rate',colnames(mut_features),colnames(meth_features),colnames(CN_features),TF_list,'EMT',Protein_list,Prog_Pathways,GeneSelectedTrans,ProtMSEnrich)

# To merge table
MergedDataEnrich=rbind(TreeEnrichment_Conf,TreeEnrichment_BatchEnrich,TreeEnrichment_Sus,TreeEnrichment_TissueEnrich,rep(0,180),matMut,rep(0,180),TreeEnrichment_CNEnrich[signCN2,])
MergedDataEnrichClean=MergedDataEnrich[c(1,9,10,17,23,24,27,28),]

TreeEnrichment_MethylationEnrichClean=TreeEnrichment_MethylationEnrich[signMet2,]

# To find trait that impacts the 3 main clusters
MethylToKeep=findTraitType(Heatamp_clust,TreeEnrichment_MethylationEnrichClean)
heatmap.2(TreeEnrichment_MethylationEnrichClean[MethylToKeep,],dendrogram='none', Rowv=TRUE, Colv=FALSE,hclustfun = function(x) hclust(x, method="ward.D2"),trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.2,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(5,8))
write.csv(rownames(TreeEnrichment_MethylationEnrichClean[MethylToKeep,]), file = "Methylation3Types.csv")

TreeEnrichment_TranscriptsClean=TreeEnrichment_TranscriptsEnrich[signTranscripts2,]
TranscriptsToKeep=findTraitType(Heatamp_clust,TreeEnrichment_TranscriptsClean)
a=heatmap.2(TreeEnrichment_TranscriptsClean[TranscriptsToKeep,],dendrogram='none', Rowv=TRUE, Colv=FALSE,distfun = function(x) dist(x, method="manhattan"),hclustfun = function(x) hclust(x, method="mcquitty"),trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.2,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(5,8))
write.csv(rownames(TreeEnrichment_TranscriptsClean[TranscriptsToKeep,]), file = "Transcripts3Types.csv")

TreeEnrichment_ProteinEnrichClean=TreeEnrichment_ProteinEnrich[signProt2,]
ProteinToKeep=findTraitType(Heatamp_clust,TreeEnrichment_ProteinEnrichClean)
b=heatmap.2(TreeEnrichment_ProteinEnrichClean[ProteinToKeep,],dendrogram='none', Rowv=TRUE, Colv=FALSE,distfun = function(x) dist(x, method="manhattan"),hclustfun = function(x) hclust(x, method="mcquitty"),trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(5,8))
write.csv(rownames(TreeEnrichment_ProteinEnrichClean[ProteinToKeep,]), file = "Protein3Types.csv")

TreeEnrichment_TFEnrichClean=TreeEnrichment_TFEnrich[signTF2,]
TFToKeep=findTraitType(Heatamp_clust,TreeEnrichment_TFEnrichClean)
c=heatmap.2(TreeEnrichment_TFEnrichClean[TFToKeep,],dendrogram='none', Rowv=TRUE, Colv=FALSE,distfun = function(x) dist(x, method="manhattan"),hclustfun = function(x) hclust(x, method="mcquitty"),trace='none',col=colourPanel,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(5,8))
write.csv(rownames(TreeEnrichment_TFEnrichClean[TFToKeep,]), file = "TF3Types.csv")

# Merged
MergedDataEnrichClean
TreeEnrichment_MethylationEnrichClean[MethylToKeep,]
TreeEnrichment_TranscriptsClean[TranscriptsToKeep,]
TreeEnrichment_ProteinEnrichClean[ProteinToKeep,]
TreeEnrichment_TFEnrichClean[TFToKeep,]
TreeEnrichment_ProgEnrich[c(7,9),]
TreeEnrichment_EMTEnrich

MergedDataAll=rbind(MergedDataEnrichClean,rep(0,180),TreeEnrichment_MethylationEnrichClean[MethylToKeep,],rep(0,180),TreeEnrichment_TranscriptsClean[TranscriptsToKeep[a[["rowInd"]]],],rep(0,180),TreeEnrichment_ProteinEnrichClean[ProteinToKeep[b[["rowInd"]]],],rep(0,180),TreeEnrichment_TFEnrichClean[TFToKeep[c[["rowInd"]]],],rep(0,180),TreeEnrichment_ProgEnrich[c(7,9),],rep(0,180),TreeEnrichment_EMTEnrich)

lwid = c(1.5,4)
lhei = c(1.5,4,1)
lmat = rbind(c(0,3),c(2,1),c(0,4))

COLOR1=c("#BB7E98","pink3",'pink2','pink', 'white',"white", 'steelblue1',"steelblue2",'steelblue3','steelblue4')

heatmap.2(MergedDataAll,dendrogram='none', Rowv=FALSE, Colv=FALSE,hclustfun = function(x) hclust(x, method="ward.D2"),trace='none',col=COLOR1,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.1,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",margins=c(3,7))
heatmap.2(MergedDataAll,dendrogram='none', Rowv=FALSE, Colv=FALSE,hclustfun = function(x) hclust(x, method="ward.D2"),trace='none',col=COLOR1,symkey=F, breaks=seq(-5,5,1),key=T,keysize = 1, density.info="none",cexRow=0.3,cexCol=0.4,key.xlab='-Log10(Corr.Pval)', na.color="gray",lmat = lmat, lwid = lwid, lhei = lhei)


### For Publication - Export
# All traits - significant and not significant
MergedData_Publication=rbind(TreeEnrichment_Conf,TreeEnrichment_BatchEnrich,TreeEnrichment_Sus,rep(0,180),TreeEnrichment_TissueEnrich,TreeEnrichment_PathologyEnrich,TreeEnrichment_Histology,TreeEnrichment_HistologySubtypeEnrich,TreeEnrichment_TypeEnrich,TreeEnrichment_EthnicityEnrich,TreeEnrichment_MutationRateEnrich,rep(0,180),TreeEnrichment_MutationEnrich,rep(0,180),TreeEnrichment_CNEnrich,rep(0,180),TreeEnrichment_MethylationEnrich,rep(0,180),TreeEnrichment_TranscriptsEnrich,rep(0,180),TreeEnrichment_ProteinEnrich,rep(0,180),TreeEnrichment_TFEnrich,rep(0,180),TreeEnrichment_ProgEnrich,rep(0,180),TreeEnrichment_EMTEnrich)
MergedData_Publication_2=10^(-abs(MergedData_Publication))
MergedData_Publication_3=(as.integer(MergedData_Publication>=0)*2-1)*MergedData_Publication_2


write.csv(MergedData_Publication_3,'TreeEnrichment_OmicsIntegration_All.csv')

# Only the significant
MergedDataAll=rbind(MergedDataEnrichClean,rep(0,180),TreeEnrichment_MutationEnrich[signMut2,],rep(0,180),TreeEnrichment_CNEnrich[signCN2,],rep(0,180),TreeEnrichment_MethylationEnrich[signMet2,],rep(0,180),TreeEnrichment_TranscriptsEnrich[signTranscripts2,],rep(0,180),TreeEnrichment_ProteinEnrich[signProt2,],rep(0,180),TreeEnrichment_TFEnrich[signTF2,],rep(0,180),TreeEnrichment_ProgEnrich[c(7,9),],rep(0,180),TreeEnrichment_EMTEnrich)
MergedDataAll_2=10^(-abs(MergedDataAll))
MergedDataAll_3=(as.integer(MergedDataAll>=0)*2-1)*MergedDataAll_2

write.csv(MergedDataAll_3,'TreeEnrichment_OmicsIntegration_Significant.csv')





