
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lipid Labeling - 13C enrichment of selected lipids   
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd('MetabCCLs/code/Lipids/LipidsLabeling')

library(readxl)
library(gplots)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(reshape2)
library(ComplexHeatmap)
source('data_summary.R')

integration_results_new_008 <- read_excel("../../../data/lipidomicsLabelingData.xls")

# Keep lipid Name, formula and position
Lipids_Formula=integration_results_new_008[1,]
posLipids=which(Lipids_Formula[1,]!='')
integration_results_new_008=integration_results_new_008[c(-1,-2),]
colnames(integration_results_new_008[,posLipids])=Lipids_Formula[1,posLipids]
integration_results_new_008[,3:21]=lapply(integration_results_new_008[,3:21],function(x) as.numeric(as.character(x)))

# Metabolic types - cell lines
integration_results_new_008$Typing[regexec('OVCAR3',integration_results_new_008$...1)==1]='Type 2'
integration_results_new_008$Typing[regexec('SW620',integration_results_new_008$...1)==1]='Type 2'
integration_results_new_008$Typing[regexec('SKMEL5',integration_results_new_008$...1)==1]='Type 1'
integration_results_new_008$Typing[regexec('HS578T',integration_results_new_008$...1)==1]='Type 1'
integration_results_new_008$Typing[regexec('NCIH460',integration_results_new_008$...1)==1]='Type 2'
integration_results_new_008$Typing[regexec('OVCAR5',integration_results_new_008$...1)==1]='Type 1'
integration_results_new_008$Typing[regexec('SF539',integration_results_new_008$...1)==1]='Type 1'
integration_results_new_008$Typing[regexec('HCT15',integration_results_new_008$...1)==1]='Type 2'
integration_results_new_008$Typing[regexec('T47D',integration_results_new_008$...1)==1]='Type 2'

integration_results_new_008$CellLines=sub('_._..*','',integration_results_new_008$...1)
integration_results_new_008$CellLinesReplicate=sub(' M.*','',integration_results_new_008$...1)
integration_results_new_008$CellLinesReplicate=as.factor(integration_results_new_008$CellLinesReplicate)


###  Compute Mass distribution vector and Fractional Contribution (FC) ###
FCPval=c()
FCval=c()
FClipid=c()
for (i in posLipids){
  integration_results_new_008$Enrichment=rep(0,length(integration_results_new_008$...1))
  integration_results_new_008$SumArea=rep(0,length(integration_results_new_008$...1))

  integration_results_new_009=integration_results_new_008
  print(i)
  
  integration_results_new_009$Enrichment
  CellUniq=unique(integration_results_new_009$CellLinesReplicate)
  for(j in 1:length(CellUniq)){
    pos=which(integration_results_new_009$CellLinesReplicate==CellUniq[j])
    Sum_area=sum(integration_results_new_009[pos,i],na.rm = T)
    integration_results_new_009$SumArea[pos]=rep(Sum_area,length(pos))
  }
  
  integration_results_new_009$Enrichment=integration_results_new_009[,i]/integration_results_new_009$SumArea
  integration_results_new_009=integration_results_new_009 %>% group_by(CellLinesReplicate) %>% mutate(SumAreaCheck=sum(Enrichment,na.rm = T))
  
  integration_13CG=integration_results_new_009[regexec('13CG',integration_results_new_009$CellLinesReplicate)>1,]
  
  integration_13CG$Typing=as.factor(integration_13CG$Typing)
  integration_13CG$...2=as.factor(integration_13CG$...2)
  integration_13CG$Enrichment=as.numeric(unlist(integration_13CG$Enrichment))
  
  df3=data_summary(integration_13CG,'Enrichment',c('Typing','...2'))
  
  # Convert dose to a factor variable
  df3$...2=round(as.numeric(sub('_M','',df3$...2)))
  df3$Isotop=as.factor(df3$...2)
  head(df3)
  df4=df3[!is.na(df3$Enrichment),]
  
  b <-  ggplot(df4, aes(x=Isotop, y=Enrichment, fill=Typing)) + 
    geom_bar(stat="identity", position=position_dodge()) +
    labs(x = "Isotopologues")+
    ggtitle(colnames(integration_results_new_009)[i])+
    geom_errorbar(aes(ymin=Enrichment-sd, ymax=Enrichment+sd), width=.2,position=position_dodge(.9))  + theme_minimal()+  scale_fill_manual(values=c("blue3", "red3"))
  
  b
  
  integration_13CG$Enrichment[is.na(integration_13CG$Enrichment)]=0
  
  if(sum(df4$Enrichment,na.rm =T)!=0){
    Carbons=sub('C','',Lipids_Formula[i])
    Carbons=as.numeric(sub('H.*','',Carbons))
    
    integration_13CG$...2=round(as.numeric(sub('_M','',integration_13CG$...2)))
    #FC=integration_results_002_13CGlucose  %>% group_by(CellLinesReplicate) %>% summarize(frac=sum(Enrichment*...2,na.rm=T)/Carbons,Typing=dplyr::first(Typing))
    
    
    CellUniq=unique(integration_13CG$CellLinesReplicate)
    FC <- data.frame(frac=double(),
                     CellLinesReplicate=character(), 
                     Typing=character())
    for(j in 1:length(CellUniq)){
      pos=which(integration_13CG$CellLinesReplicate==CellUniq[j])
      frac=sum(integration_13CG$Enrichment[pos]*integration_13CG$...2[pos],na.rm=T)/Carbons
      val=data.frame(frac,CellLinesReplicate=CellUniq[j],Typing=integration_13CG$Typing[pos][1])
      FC=rbind(FC,val)
    }
    
    p<-ggplot(FC,aes(x = Typing, y=frac,fill=Typing)) + 
      geom_boxplot(width=0.5)+ 
      theme_classic()+
      # xlab("class") +
      xlab("Metabolic Types") +
      ylab("Fractional contribution")+
      ggtitle(paste0(colnames(integration_results_new_008)[i],'   N = ',length(FC$frac)))+ 
      scale_fill_manual(values=c("blue3", "red3"))+
      stat_compare_means(method='t.test')
    p
    
    compFC=t.test(FC$frac[FC$Typing=='Type 1'],FC$frac[FC$Typing=='Type 2'])
    FCPval=cbind(FCPval,compFC$p.value)
    FCval=cbind(FCval,FC$frac)
    FClipid=cbind(FClipid,colnames(integration_results_new_008)[i])
    print(p)
    
  }
}



###  To plot fractional contribution - Figure 5E  ###

rownames(FCval)=FC$CellLinesReplicate
colnames(FCval)=FClipid

FCval_2=melt(data =FCval, id.vars = "Lipids", measure.vars = FClipid)
FCval_2$Typing[regexec('OVCAR3',FCval_2$Var1)==1]='Type 2'
FCval_2$Typing[regexec('SW620',FCval_2$Var1)==1]='Type 2'
FCval_2$Typing[regexec('SKMEL5',FCval_2$Var1)==1]='Type 1'
FCval_2$Typing[regexec('HS578T',FCval_2$Var1)==1]='Type 1'
FCval_2$Typing[regexec('NCIH460',FCval_2$Var1)==1]='Type 2'
FCval_2$Typing[regexec('OVCAR5',FCval_2$Var1)==1]='Type 1'
FCval_2$Typing[regexec('SF539',FCval_2$Var1)==1]='Type 1'
FCval_2$Typing[regexec('HCT15',FCval_2$Var1)==1]='Type 2'
FCval_2$Typing[regexec('T47D',FCval_2$Var1)==1]='Type 2'

FCval_2$Var2=sub(';.*','',FCval_2$Var2)
FCval_2$Var2=sub('\\[.*','',FCval_2$Var2)
FCval_2_stat=data_summary(FCval_2,'value',c('Typing','Var2'))

FCval_2_stat%>% 
  #filter(Var2 != "TAG 50:1 -correct RT ") %>%
  filter(Var2 != 'TAG 50:1 -old RT') %>%
  filter(Var2 != 'LPE 16:0 ') %>%
  ggplot(aes(x=Typing,y=value*100,fill = Typing))+
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=100*(value-sd), ymax=100*(value+sd)), width=.2,position=position_dodge(.9)) +
  facet_wrap(~Var2,ncol = 15)+
  scale_fill_manual(values=c("blue3", "red3"))+
  stat_compare_means(data=subset(FCval_2,!Var2 %in% c( "TAG 50:1 -old RT",'LPE 16:0 ')),method='t.test',label.x.npc = "center",label = "p.signif",symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns")))+
  theme_minimal()+
  ylab('Fractional Labeling [%]')+
  xlab('')



# To plot Single Cell line
###  MDV Selected Cell line - Figure 5D  ###

i=10 # TAG 50:1 Selected
integration_results_new_008$Enrichment=rep(0,length(integration_results_new_008$...1))
integration_results_new_008$SumArea=rep(0,length(integration_results_new_008$...1))
integration_results_new_009=integration_results_new_008

integration_results_new_009$Enrichment
CellUniq=unique(integration_results_new_009$CellLinesReplicate)
for(j in 1:length(CellUniq)){
  pos=which(integration_results_new_009$CellLinesReplicate==CellUniq[j])
  Sum_area=sum(integration_results_new_009[pos,i],na.rm = T)
  integration_results_new_009$SumArea[pos]=rep(Sum_area,length(pos))
}

integration_results_new_009$Enrichment=integration_results_new_009[,i]/integration_results_new_009$SumArea
integration_results_new_009=integration_results_new_009 %>% group_by(CellLinesReplicate) %>% mutate(SumAreaCheck=sum(Enrichment,na.rm = T))
integration_13CG=integration_results_new_009[regexec('13CG',integration_results_new_009$CellLinesReplicate)>1,]
integration_13CG$Typing=as.factor(integration_13CG$Typing)
integration_13CG$...2=as.factor(integration_13CG$...2)
integration_13CG$Enrichment=as.numeric(unlist(integration_13CG$Enrichment))
df3=data_summary(integration_13CG,'Enrichment',c('CellLines','...2'))

# Convert dose to a factor variable
df3$...2=round(as.numeric(sub('_M','',df3$...2)))
df3$Isotop=as.factor(df3$...2)
head(df3)
df4=df3[!is.na(df3$Enrichment),]
Carbons=sub('C','',Lipids_Formula[i])
Carbons=as.numeric(sub('H.*','',Carbons))

df4$CellLines=sub('T47D_13CG','T47D',df4$CellLines)
df4 %>%
  filter(CellLines == "T47D") %>%
  ggplot(mapping=aes(x=Isotop, y=Enrichment, fill=CellLines)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Isotopologues")+
  ggtitle('TAG 50:1')+
  geom_errorbar(aes(ymin=Enrichment-sd, ymax=Enrichment+sd), width=.2,position=position_dodge(.9))  + theme_minimal()+  scale_fill_manual(values=c("red3"))+
  scale_x_discrete( breaks=seq(1,Carbons,2)) 

ggsave('TAG_T47D.pdf',plot = last_plot(),width = 7.54, height = 2.61)

df4$CellLines=sub('HS578T_13CG','HS578T',df4$CellLines)
df4 %>%
  filter(CellLines == "HS578T") %>%
  ggplot(mapping=aes(x=Isotop, y=Enrichment, fill=CellLines)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "Isotopologues")+
  ggtitle("TAG 50:1")+
  geom_errorbar(aes(ymin=Enrichment-sd, ymax=Enrichment+sd), width=.2,position=position_dodge(.9))  + theme_minimal()+  scale_fill_manual(values=c("blue3"))+
  scale_x_discrete( breaks=seq(1,Carbons,2)) 

ggsave('TAG_HS578T.pdf',plot = last_plot(),width = 7.54, height = 2.61)

# Compute mean fractional contribution of selected cell lines
# For HS578T
x=FCval[7:(7+5),7]
mean(x)
sem<-sd(x)/sqrt(length(x))

# For T47D
x=FCval[48:53,7]
mean(x)
sem<-sd(x)/sqrt(length(x))
