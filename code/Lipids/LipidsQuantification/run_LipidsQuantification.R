

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lipid Quantification using internal standard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd('MetabCCLs/code/Lipids/LipidsQuantification')

library(ggpubr)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(tidyverse)
library(unikn)
library(circlize)
library(colorspace)
library(GetoptLong)
library(gplots)

## Load Lipidomics Data
Cancer_metabotypes_pos_quant_AO <-
  read.csv("../../../data/lipidomicsData.csv")
Cancer_metabotypes_pos_quant_AO$Class = trimws(Cancer_metabotypes_pos_quant_AO$Class, which =
                                                 'both')
colnames(Cancer_metabotypes_pos_quant_AO) = sub('_NoLab.*', '', colnames(Cancer_metabotypes_pos_quant_AO))

## Cell Line types
# Not using HCT15 and T47D because spike in standard were not added at same time - absolute values not comparable
posType1 = which(!is.na(match(
  colnames(Cancer_metabotypes_pos_quant_AO),
  c('SF539', 'OVCAR5', 'HS578T', 'SKMEL5')
)))
colnames(Cancer_metabotypes_pos_quant_AO)[posType1]
posType2 = which(!is.na(match(
  colnames(Cancer_metabotypes_pos_quant_AO),
  c('OVCAR3', 'SW620', 'NCIH460')
)))
colnames(Cancer_metabotypes_pos_quant_AO)[posType2]

# Lipid class
Sphingolipid = c('HexCer', 'SM')
Phospholipid = c('CL', 'PC', 'LPC', 'Ether PC', 'PE', 'LPE', 'PS', 'PG')
Glycerolipids = c('TAG')
Merged = c(Sphingolipid, Phospholipid, Glycerolipids)

# Colors - Lipid Class
my_pal <- seecol(pal_unikn_pair, 20)
my_pal2 = my_pal[1:11]
my_pal2[c(1, 2)] = my_pal[c(19, 18)]
my_pal2[c(3)] = my_pal[c(16)]
my_pal2[c(4, 5)] = my_pal[c(14, 12)]
my_pal2[c(6)] = my_pal[c(11)]
my_pal2[c(7, 8)] = my_pal[c(9, 6)]
my_pal2[c(9, 10)] = my_pal[c(5, 3)]
my_pal2[c(11)] = my_pal[c(1)]
names(my_pal2) = Merged

## Compute Sum ALL species
sumType1 = colSums(Cancer_metabotypes_pos_quant_AO[, posType1], na.rm = TRUE)
sumType2 = colSums(Cancer_metabotypes_pos_quant_AO[, posType2], na.rm = TRUE)
AllSpecies = data.frame(SumALLSpecies = c(sumType1, sumType2),
                        Typing = c(rep('Type 1', length(sumType1)), rep('Type 2', length(sumType2))))

## Compute sum all class
sumALLType1 = aggregate(
  Cancer_metabotypes_pos_quant_AO[, posType1],
  by = list(Category = Cancer_metabotypes_pos_quant_AO$Class),
  FUN = sum
)
sumALLType2 = aggregate(
  Cancer_metabotypes_pos_quant_AO[, posType2],
  by = list(Category = Cancer_metabotypes_pos_quant_AO$Class),
  FUN = sum
)

# Number of each lipids
Cancer_metabotypes_pos_quant_AO[, 1:5] %>% count(Class)
dplyr::summarise(Cancer_metabotypes_pos_quant_AO, number = distinct(Class))
sumALLType1_t = tidyr::gather(sumALLType1, 2:25, key = 'CellLines', value =
                                'SUMGroup')
sumALLType1_t = cbind(sumALLType1_t, Typing = rep('Type 1', length(sumALLType1_t)))
sumALLType2_t = tidyr::gather(sumALLType2, 2:19, key = 'CellLines', value =
                                'SUMGroup')
sumALLType2_t = cbind(sumALLType2_t, Typing = rep('Type 2', length(sumALLType2_t)))

sumALLTyping = rbind(sumALLType1_t, sumALLType2_t)

#### To plot only few examples - Sup Figure 6 ###
sumALLTyping %>%
  filter(Category != 'HexCer') %>%
  filter(Category != 'LPE') %>%
  filter(Category != 'LPC') %>%
  filter(Category != 'PE') %>%
  filter(Category != 'TAG') %>%
  filter(Category != 'Ether PC') %>%
  filter(Category != 'PC') %>%
  ggplot(aes(
    x = Typing,
    y = log10(SUMGroup),
    fill = Typing
  )) +
  geom_boxplot() +
  geom_jitter(width = 0.3,
              alpha = 1,
              size = 1) +
  xlab("Typing") +
  ylab("log10(Sum - pmole/1e6 cell)") +
  facet_wrap( ~ Category, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("blue3", "red3")) +
  stat_compare_means(
    method = 't.test',
    label.x.npc = "center",
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  theme_minimal() +
  ggtitle("Total Lipids per Class") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(legend.position = "bottom")

#### To plot only few examples - Figure 5C ###
sumALLTyping %>%
  filter(Category != 'CL') %>%
  filter(Category != 'PG') %>%
  filter(Category != 'PS') %>%
  filter(Category != 'SM') %>%
  ggplot(aes(
    x = Typing,
    y = log10(SUMGroup),
    fill = Typing
  )) +
  geom_boxplot() +
  geom_jitter(width = 0.3,
              alpha = 1,
              size = 1) +
  xlab("Typing") +
  ylab("log10(Sum - pmole/1e6 cell)") +
  facet_wrap( ~ Category, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("blue3", "red3")) +
  stat_compare_means(
    method = 't.test',
    label.x.npc = "center",
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  theme_minimal() +
  ggtitle("Total Lipids per Class") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(legend.position = "bottom")


###  Compare Total Lipid Content - Figure 5B  ###
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))

sumALLTyping$Category = trimws(sumALLTyping$Category, which = 'both')
sumALLTyping = sumALLTyping %>%
  distinct(Typing, SUMGroup, CellLines, Category) %>%
  mutate(Category = fct_relevel(Category, Merged)) %>%
  arrange(Category)

pval_ttest_total = t.test(sumType1, sumType2)

p <- ggbarplot(
  sumALLTyping,
  x = "Typing",
  y = "SUMGroup",
  add = c("mean_se"),
  color = "Category",
  fill = "Category",
  palette = my_pal2,
  alpha = 0.8
)
p + theme(legend.position = c(1, 0.2)) + theme_minimal() +
  ylab("Sum(pmole/million of cells)") +
  xlab("") +
  ggtitle(paste0("Total Lipids pval= ", pval_ttest_total$p.value))



### Bubble plot - FoldChange - Figure 5A ###
# Mean concentration of each lipids
meanConcentration = rowMeans(Cancer_metabotypes_pos_quant_AO[, c(posType1, posType2)],
                             na.rm = FALSE,
                             dims = 1)
FC = log2(
  rowMeans(Cancer_metabotypes_pos_quant_AO[, posType1]) / rowMeans(Cancer_metabotypes_pos_quant_AO[, posType2])
)
stat <-
  sapply(1:nrow(Cancer_metabotypes_pos_quant_AO), function(i)
    t.test(as.numeric(as.character(
      unlist(Cancer_metabotypes_pos_quant_AO[i, posType1])
    )), as.numeric(as.character(
      unlist(Cancer_metabotypes_pos_quant_AO[i, posType2])
    )))[c("p.value")])
stat_corr = p.adjust(stat, method = "fdr")
df_FC_Pval = data.frame(
  LipidClass = Cancer_metabotypes_pos_quant_AO$Class,
  FoldChange = FC,
  CorrPval = stat_corr,
  meanConcentration = meanConcentration,
  LipidName = Cancer_metabotypes_pos_quant_AO[, 1]
)

df_FC_Pval = df_FC_Pval %>%
  distinct(LipidClass, FoldChange, CorrPval, meanConcentration, LipidName) %>%
  mutate(LipidClass = fct_relevel(LipidClass, Merged)) %>%
  arrange(LipidClass)

# Alphabetic order
df_FC_Pval = cbind(df_FC_Pval, Position = length(stat_corr):1)

g = ggplot(
  df_FC_Pval,
  aes(
    x = FoldChange,
    y = Position,
    size = meanConcentration,
    color = LipidClass
  )
) +
  geom_point(alpha = 0.4, stat = "identity") +
  geom_point(data = df_FC_Pval[df_FC_Pval$CorrPval < 0.05 &
                                 abs(df_FC_Pval$FoldChange) > 1, ],
             aes(x = FoldChange, y = Position),
             alpha = 1) +
  geom_vline(
    aes(xintercept = 0),
    color = "black",
    linetype = "dashed",
    size = 0.5
  ) + theme_minimal() + theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  geom_vline(
    aes(xintercept = 1),
    color = "blue",
    linetype = "dashed",
    size = 0.5
  ) +
  geom_vline(
    aes(xintercept = -1),
    color = "red",
    linetype = "dashed",
    size = 0.5
  ) +
  # scale_size_area()+
  # scale_size_binned()+
  scale_size(
    range = c(3, 25),
    name = "mean Concentration",
    breaks = c(
      min(meanConcentration),
      mean(meanConcentration) - sd(meanConcentration),
      mean(meanConcentration),
      mean(meanConcentration) + sd(meanConcentration),
      max(meanConcentration)
    )
  ) +
  # scale_size_manual(values =c(1,4,5,6,7,8,9,10), name="mean Concentration",breaks = c(5,10))+
  ylab("Lipid Species") +
  xlab("log2(Fold Change)") + guides(colour = guide_legend(override.aes = list(size =
                                                                                 5))) +
  # scale_color_viridis(discrete = TRUE, option = "D", alpha = 1, begin =0, end = 1, direction = -1)
  scale_colour_manual(values = my_pal2)
g


### Number of Carbon vs. Double Bonds - Sup. Figure 6 ###
LipidLCass = unique(Cancer_metabotypes_pos_quant_AO$Class)
# For lipid class
for (k in 1:length(LipidLCass)) {
  pdf(file = paste0(LipidLCass[k], "_Double_Carbon.pdf"))
  
  # Lipid Pos
  posClass = LipidLCass[k] == Cancer_metabotypes_pos_quant_AO$Class
  Cancer_metabotypes_Lipid = aggregate(Cancer_metabotypes_pos_quant_AO[posClass, 6:57],
                                       by = list(Cancer_metabotypes_pos_quant_AO[posClass, 3]),
                                       sum) ### CHANGE here for
  
  type1 = c('SF539', 'OVCAR5', 'HS578T', 'SKMEL5')
  posType1_new = grep(paste(type1, collapse = "|"),
                      colnames(Cancer_metabotypes_Lipid))
  
  type2 = c('OVCAR3', 'SW620', 'NCIH460')
  posType2_new = grep(paste(type2, collapse = "|"),
                      colnames(Cancer_metabotypes_Lipid))
  
  FC_Lipid = log2(rowMeans(Cancer_metabotypes_Lipid[, posType1_new]) / rowMeans(Cancer_metabotypes_Lipid[, posType2_new]))
  
  df_FC_Lipid = data.frame(LipidName = Cancer_metabotypes_Lipid$Group.1, FoldChange =
                             FC_Lipid)
  
  # Get number of Carbon and double bond
  df_FC_Lipid$LipidName = str_replace(df_FC_Lipid$LipidName, 'd', '')
  df_FC_Lipid$LipidName = str_replace(df_FC_Lipid$LipidName, 'e', '')
  CarbonDouble = strsplit(as.character(trimws(df_FC_Lipid$LipidName, 'both')), ':')
  CarbonDouble = matrix(unlist(CarbonDouble),
                        byrow = T,
                        nrow = length(CarbonDouble))
  
  DoubleBonds = min(CarbonDouble[, 2]):max(CarbonDouble[, 2])
  Carbons = min(CarbonDouble[, 1]):max(CarbonDouble[, 1])
  
  # Data Frame full of NA
  df <-
    data.frame(matrix(ncol = length(DoubleBonds), nrow = length(Carbons)))
  row.names(df) = Carbons
  names(df) = DoubleBonds
  
  
  for (i in 1:dim(df_FC_Lipid)[1]) {
    posCD = CarbonDouble[i, ]
    posC = which(names(df) == posCD[2])
    posD = which(row.names(df) == posCD[1])
    df[posD, posC] = df_FC_Lipid$FoldChange[i]
  }
  
  redblue_scale = redblue(8)
  redblue_scale[2] = 'red1'
  redblue_scale[8] = 'blue3'
  redblue_scale[7] = 'blue1'
  redblue_scale[1] = 'red3'
  
  col_fun_2 = colorRamp2(c(-6, -4, -2, 0, 0, 2, 4, 6), c(redblue_scale))
  heatmap_lipid = Heatmap(
    as.matrix(df),
    name = "log2(Fold Change)",
    col = colorRamp2(
      breaks = c(-6, -4, -2, 0, 0, 2, 4, 6),
      colors = c(redblue_scale)
    ),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    column_names_side = "top",
    row_names_side = "left",
    column_title = paste0(LipidLCass[k], " - Double bond number"),
    row_title = "Carbon Number",
    heatmap_legend_param = list(
      at = c(-6, -4, -2, 0, 2, 4, 6),
      labels = c('-6  Type 2 ' , -4, -2, 0, 2, 4, '6  Type 1')
    ),
    cell_fun = function(j, i, x, y, w, h, col) {
      # add text to each grid
      grid.text(if (!is.na(df[i, j])) {
        round(df[i, j])
      }, x, y)
    }
  )
  draw(heatmap_lipid)
  dev.off()
}


### Lipid Unsaturation Index ###
LipidLCass = unique(Cancer_metabotypes_pos_quant_AO$Class)

# Empty dataframe for result
sat_index_df_merged = data.frame()

numAcylChain = c(4, 2, 1, 1, 2, 2, 2, 2, 2, 2, 3)
for (k in 1:length(LipidLCass)) {
  # For all lipid class
  posClass = LipidLCass[k] == Cancer_metabotypes_pos_quant_AO$Class
  
  # Index is the weighted average - Get cell lines
  Cancer_metabotypes_Lipid_index = Cancer_metabotypes_pos_quant_AO[posClass, 6:57]
  
  formula_lipid = str_replace(Cancer_metabotypes_pos_quant_AO[posClass, 3], 'd', '')
  formula_lipid = str_replace(formula_lipid, 'e', '')
  CarbonDouble = strsplit(as.character(trimws(formula_lipid, 'both')), ':')
  CarbonDouble = matrix(unlist(CarbonDouble),
                        byrow = T,
                        nrow = length(CarbonDouble))
  # wEIGHTED AVERAGE
  saturation_index = apply(Cancer_metabotypes_Lipid_index,
                           2,
                           weighted.mean,
                           x = as.numeric(CarbonDouble[, 2]))
  
  # Per Acyl Chain
  saturation_index = saturation_index / numAcylChain[k]
  
  type1 = c('SF539', 'OVCAR5', 'HS578T', 'SKMEL5')
  posType1_new = grep(paste(type1, collapse = "|"), names(saturation_index))
  
  type2 = c('OVCAR3', 'SW620', 'NCIH460')
  posType2_new = grep(paste(type2, collapse = "|"), names(saturation_index))
  
  sat_index_df = data.frame(
    saturation_index = c(saturation_index[posType1_new], saturation_index[posType2_new]),
    Typing = c(rep('Type 1', length(posType1_new)), rep('Type 2', length(posType2_new))),
    LipidClass = rep(LipidLCass[k], length(posType2_new) + length(posType1_new)),
    CellLines = c(names(saturation_index)[posType1_new], names(saturation_index)[posType2_new])
  )
  sat_index_df_merged = rbind(sat_index_df_merged, sat_index_df)
}

### Lipid Unsaturation Index per Class - Sup. Figure 7A ###
sat_index_df_merged %>%
  ggplot(aes(x = Typing, y = saturation_index, fill = Typing)) +
  geom_boxplot() +
  geom_jitter(width = 0.3,
              alpha = 1,
              size = 1) +
  xlab("") +
  ylab("Average double bond per acyl chain") +
  facet_wrap( ~ LipidClass, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("blue3", "red3")) +
  stat_compare_means(
    method = 't.test',
    label.x.npc = "center",
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  #stat_compare_means(method='t.test',label.x.npc = "center")+
  theme_minimal() +
  ggtitle("Lipid unsaturation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

### Lipid Unsaturation Index global - Sup. Figure 7B ###
sat_index_df_merged_mean = sat_index_df_merged %>%
  group_by(CellLines, Typing) %>%
  summarise(saturation_index = mean(saturation_index),
            n = n())
# Box plot
sat_index_df_merged_mean %>%
  ggplot(aes(x = Typing, y = saturation_index, fill = Typing)) +
  geom_boxplot() +
  geom_jitter(width = 0.3,
              alpha = 1,
              size = 1) +
  xlab("") +
  ylab("Average double bond per acyl chain") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("blue3", "red3")) +
  stat_compare_means(
    method = 't.test',
    label.x.npc = "center",
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  #stat_compare_means(method='t.test',label.x.npc = "center")+
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Lipid unsaturation all Class") +
  theme(plot.title = element_text(hjust = 0.5))




### Look into Acylchain length ###
LipidLCass = unique(Cancer_metabotypes_pos_quant_AO$Class)

# Where to store results
chain_index_df_merged = data.frame()

numAcylChain = c(4, 2, 1, 1, 2, 2, 2, 2, 2, 2, 3)
for (k in 1:length(LipidLCass)) {
  # For all lipid class
  posClass = LipidLCass[k] == Cancer_metabotypes_pos_quant_AO$Class
  
  # Index is the weighted average - Get cell lines
  Cancer_metabotypes_Lipid_index = Cancer_metabotypes_pos_quant_AO[posClass, 6:57]
  
  formula_lipid = str_replace(Cancer_metabotypes_pos_quant_AO[posClass, 3], 'd', '')
  formula_lipid = str_replace(formula_lipid, 'e', '')
  CarbonDouble = strsplit(as.character(trimws(formula_lipid, 'both')), ':')
  CarbonDouble = matrix(unlist(CarbonDouble),
                        byrow = T,
                        nrow = length(CarbonDouble))
  # wEIGHTED AVERAGE
  carbon_index = apply(Cancer_metabotypes_Lipid_index,
                       2,
                       weighted.mean,
                       x = as.numeric(CarbonDouble[, 1]))
  
  # Per Acyl Chain
  carbon_index = carbon_index / numAcylChain[k]
  
  type1 = c('SF539', 'OVCAR5', 'HS578T', 'SKMEL5')
  posType1_new = grep(paste(type1, collapse = "|"), names(carbon_index))
  
  type2 = c('OVCAR3', 'SW620', 'NCIH460')
  posType2_new = grep(paste(type2, collapse = "|"), names(carbon_index))
  
  chain_index_df = data.frame(
    carbon_index = c(carbon_index[posType1_new], carbon_index[posType2_new]),
    Typing = c(rep('Type 1', length(posType1_new)), rep('Type 2', length(posType2_new))),
    LipidClass = rep(LipidLCass[k], length(posType2_new) + length(posType1_new)),
    CellLines = c(names(carbon_index)[posType1_new], names(carbon_index)[posType2_new])
  )
  chain_index_df_merged = rbind(chain_index_df_merged, chain_index_df)
}

chain_index_df_merged %>%
  ggplot(aes(x = Typing, y = carbon_index , fill = Typing)) +
  geom_boxplot() +
  geom_jitter(width = 0.3,
              alpha = 1,
              size = 1) +
  xlab("") +
  ylab("Average length per acyl chain") +
  facet_wrap( ~ LipidClass, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("blue3", "red3")) +
  stat_compare_means(
    method = 't.test',
    label.x.npc = "center",
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  #stat_compare_means(method='t.test',label.x.npc = "center")+
  theme_minimal() +
  ggtitle("Acyl Chain Length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Merged Class
chain_index_df_merged_mean = chain_index_df_merged %>%
  group_by(CellLines, Typing) %>%
  summarise(carbon_index = mean(carbon_index), n = n())

# Box plot
chain_index_df_merged_mean %>%
  ggplot(aes(x = Typing, y = carbon_index, fill = Typing)) +
  geom_boxplot() +
  geom_jitter(width = 0.3,
              alpha = 1,
              size = 1) +
  xlab("") +
  ylab("Average length per acyl chain") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("blue3", "red3")) +
  stat_compare_means(
    method = 't.test',
    label.x.npc = "center",
    label = "p.signif",
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  #stat_compare_means(method='t.test',label.x.npc = "center")+
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Acyl Chain Length") +
  theme(plot.title = element_text(hjust = 0.5))
