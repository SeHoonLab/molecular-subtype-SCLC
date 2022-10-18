.libPaths("/home/jjg/tools/RprofileLibpath/")
require(ggplot2)
library(reshape2)
library(ggsci)
library("xlsx")
library(tidyverse)
library(ggsci)
library(cowplot)
library(EnhancedVolcano)
library(ggsci)
library(fgsea)
library(scales)
library(singscore)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
require(survminer)
require(survival)

axis_theme = theme(axis.text.x = element_text(size=10, face="bold", color = "black"),
              axis.title.x = element_text(size=10, face="bold", color = "black"),
        axis.text.y = element_text(size=10, face="bold", color = "black"),
        axis.title.y = element_text(size=10, face="bold", color = "black"), panel.border = element_rect(size = 1))

work_dir = "/BiO2/users/jjg/SCLC/I_SCLC/Paper/221011_Github_code/Fig4/"
SCLC_meta = read.xlsx("/BiO2/users/jjg/SCLC/I_SCLC/Paper/221011_Github_code/Essential_check_Edited.xlsx", header = TRUE, sheetName = "Sheet1")

SCLC_meta$SGI_ID = lapply(SCLC_meta$SGI_ID, function(x) {strsplit(gsub("^ ", "", x), " ")[[1]][1]}) %>% unlist() %>% as.character()
SCLC_subtype_WTS_meta = SCLC_meta[!is.na(SCLC_meta$SCLC_subtype) & c(SCLC_meta$WTS_QC_Result %in% c("Pass", "1")) & grepl("includ", SCLC_meta$study.inclusion),]
SCLC_subtype_WTS_meta$NE_subtype = lapply(SCLC_subtype_WTS_meta$SCLC_subtype, function(x) {
                                                           if(x=="A" | x=="N" | x=="trans"){
                                                                return("NE")
                                                           }else if(x=="P" | x=="TN"){
                                                                return("Non-NE")
                                                           }
}) %>% as.character()
SCLC_subtype_WTS_meta$SCLC_subtype = factor(SCLC_subtype_WTS_meta$SCLC_subtype, levels = c("A", "trans", "N", "P", "TN"))

pathways.hallmark <- gmtPathways("/BiO2/users/jjg/EGS/Figure/Fig2/GeneSet/db/h.all.v7.4.symbols.gmt")
pathway.kegg = gmtPathways("/BiO2/users/jjg/EGS/Figure/Fig2/GeneSet/db/c2.cp.kegg.v7.4.symbols.gmt")

SCLC_cpm =read.table("/BiO2/users/jjg/SCLC/I_SCLC/Paper/221011_Github_code/expr_matrix/log2_CPM_n226.txt", sep = "\t")

NE_25_genelist = c("BEX1","ASCL1","INSM1","CHGA","TAGLN3","KIF5C","CRMP1","SCG3","SYT4","RTN1","MYT1","SYP","KIF1A","TMSB15A","SYN1","SYT11","RUNDC3A","TFF3","CHGB","FAM57B","SH3GL2","BSN","SEZ6","TMSB15B","CELF3")
nonNE_25_genelist = c("RAB27B","TGFBR2","SLC16A5","S100A10","ITGB4","YAP1","LGALS3","EPHA2","S100A16","PLAU","ABCC3","ARHGDIB","CYR61","PTGES","CCND1","IFITM2","IFITM3","AHNAK","CAV2","TACSTD2","TGFBI","EMP1","CAV1","ANXA1","MYOF")

dataPath = "/BiO2/users/jjg/SCLC/I_SCLC/220214_JTO_NK_CD8T_review/NK_scoring/R implementation/data/"
##------- Read in the NK signature (from Supplementary Table S1)
nk_signature <- read.csv(paste0(dataPath, "Cursons_Guimaraes_NKsignature_CIR_2019.csv"),
 stringsAsFactors = F)

nk <- as.character(nk_signature$HGNC.Symbol[nk_signature$Cursons.Guimaraes.sigGene == "TRUE"])
T_cell_Inflamed = c("PSMB10", "HLA-DQA1", "HLA-DRB1", "CMKLR1", "HLA-E", "NKG7", "CD8A", "CCL5", "CXCL9", "CD27", "CXCR6", "IDO1", "STAT1", "TIGIT", "LAG3", "CD274", "PDCD1LG2", "CD276")
SupFig5B = read.table("/BiO2/users/jjg/SCLC/I_SCLC/220404_WT_silent_non_damaging_and_damage/3.WTS_with_IHC_182/CancerCell_SupFig5B.txt", header = FALSE, sep = "\t")
Exp.count = read.table("./expr_matrix/Exp.count.txt", header = TRUE, sep = '\t')

######################################################
################################# fsea plot
######################################################

set.seed(100)
N_Atrans_DEG = readRDS("Fig_S5B_edgeR_N_Atrans_DEG.Rds")
N_Atrans_DEG = N_Atrans_DEG$table
#write.table(N_Atrans_DEG, "/BiO2/users/jjg/SCLC/I_SCLC/Paper/220608/Fig4/220609_plot_N_versus_A&trans_DEG.txt", col.names = TRUE, row.names = TRUE, sep = '\t', quote = FALSE)
#N_Atrans_DEG$Log2FC = log2(exp(N_Atrans_DEG$logFC)) ####logical, if TRUE then log2 values are returned.

Filt_DEG = N_Atrans_DEG[N_Atrans_DEG$PValue<=1,]

Filt_DEG$fcsign <- sign(Filt_DEG$logFC)
Filt_DEG$logP=-log10(Filt_DEG$PValue)
Filt_DEG$metric= Filt_DEG$logP/Filt_DEG$fcsign


Rank_DEG = data.frame(rownames(Filt_DEG), Filt_DEG$metric)
colnames(Rank_DEG) = c("Gene", "Rank")


fgseaRes <- fgsea(pathways = pathways.hallmark,
                  stats    = deframe(Rank_DEG),
                  minSize  = 15,
                  maxSize  = 500)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaRes$pathway = gsub("HALLMARK_", "", fgseaRes$pathway)
fgseaRes$pathway = gsub("_", " ", fgseaRes$pathway)

fgseaRes$log10pval = -log10(fgseaRes$pval)
fgsea_N_Atrans = fgseaRes

p = ggplot(fgseaRes[grepl("NOTCH", fgseaRes$pathway) | (!is.na(fgseaRes$pval) & abs(fgseaRes$NES)>=1 & fgseaRes$pval < 0.1),]) +
  geom_bar(stat = 'identity', aes(reorder(pathway, NES), NES, fill = log10pval)) +
  coord_flip()  +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0("Hallmark pathways NES from GSEA(N vs A&trans)")) +
  theme_bw() + scale_fill_gradient(low = "#ff0000", high = "#0000FF") + axis_theme +  scale_fill_gradientn(limits = c(0,3), colours=c("blue", "red"))

ggsave(paste0(work_dir, "Fig_4B_N_Atrans_DEG_pvalue_fgsea.png"), plot = p, width = 9, height=5)



#P_vs A&trans

set.seed(100)
P_Atrans_DEG = readRDS("Fig_S5A_edgeR_P_Atrans_DEG.Rds")
P_Atrans_DEG = P_Atrans_DEG$table

Filt_DEG = P_Atrans_DEG[P_Atrans_DEG$PValue<=1,]

Filt_DEG$fcsign <- sign(Filt_DEG$logFC)
Filt_DEG$logP=-log10(Filt_DEG$PValue)
Filt_DEG$metric= Filt_DEG$logP/Filt_DEG$fcsign


Rank_DEG = data.frame(rownames(Filt_DEG), Filt_DEG$metric)
colnames(Rank_DEG) = c("Gene", "Rank")

fgseaRes <- fgsea(pathways = pathways.hallmark,
                  stats    = deframe(Rank_DEG),
                  minSize  = 15,
                  maxSize  = 500)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaRes$pathway = gsub("HALLMARK_", "", fgseaRes$pathway)
fgseaRes$pathway = gsub("_", " ", fgseaRes$pathway)
fgseaRes$log10pval = -log10(fgseaRes$pval)
fgsea_P_Atrans = fgseaRes


p = ggplot(fgseaRes[grepl("NOTCH", fgseaRes$pathway) | c(!is.na(fgseaRes$pval) & abs(fgseaRes$NES)>=1 & fgseaRes$pval < 0.1) & !(fgseaRes$pathway%in% c("KRAS SIGNALING DN", "ESTROGEN RESPONSE EARLY")),]) +
  geom_bar(stat = 'identity', aes(reorder(pathway, NES), NES, fill = log10pval)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0("Hallmark pathways NES from GSEA(P vs A&trans)")) +
  theme_bw() + scale_fill_gradient(low = "#ff0000", high = "#0000FF") + axis_theme + scale_fill_gradientn(limits = c(0,3), colours=c("blue", "red"))

ggsave(paste0(work_dir, "Fig_4A_P_A&trans_DEG_pvalue_fgsea.png"), plot = p, widt=8.4, height=5)




#########


SCLC_cpm2 = SCLC_cpm[,SCLC_subtype_WTS_meta$WTS_ID]




SCLC_cpm_Rank <- rankGenes(SCLC_cpm2)


SingScore <- simpleScore(
  rankData = SCLC_cpm_Rank,
  upSet = NE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$NE_25_genelist = rescale(SingScore$TotalScore, to = c(-4, 4))

SingScore <- simpleScore(
  rankData = SCLC_cpm_Rank,
  upSet = nonNE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$nonNE_25_genelist = rescale(SingScore$TotalScore, to = c(-4, 4))


################################################
##### Fig 4C



SCLC_subtype_WTS_meta_NE_P = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$NE_subtype=="NE" | SCLC_subtype_WTS_meta$SCLC_subtype == "P",]

pt.matrix = SCLC_cpm2[,SCLC_subtype_WTS_meta_NE_P$WTS_ID]
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- rownames(SCLC_cpm2); colnames(pt.matrix) = colnames(SCLC_cpm2)[SCLC_subtype_WTS_meta$NE_subtype=="NE" | SCLC_subtype_WTS_meta$SCLC_subtype == "P"]



tmp = SCLC_subtype_WTS_meta_NE_P
tmp = tmp[order(match(tmp$SCLC_subtype, c("A", "trans", "N", "P"))),]
pt.matrix = pt.matrix[,tmp$WTS_ID]
annotation = data.frame(tmp$SCLC_subtype, tmp$NE_25_genelist, tmp$nonNE_25_genelist)
rownames(annotation) = tmp$WTS_ID
colnames(annotation) = c("SCLC_subtype", "NE Signature", "non-NE Signature")
annotation = annotation[colnames(pt.matrix),]
annotation = data.frame(annotation)
colnames(annotation) = c("SCLC_subtype", "NE Signature", "non-NE Signature")
rownames(annotation) = colnames(pt.matrix)
annotation = annotation[,c("NE Signature", "non-NE Signature", "SCLC_subtype")]

col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white","red"))


annotation$SCLC_subtype = factor(annotation$SCLC_subtype, level=c("A", "trans", "N", "P"))
colAnn <- HeatmapAnnotation(df = annotation,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col = list('NE Signature' = col_fun,
             'non-NE Signature' = col_fun,
             SCLC_subtype = c("A" = "#BC3C29", "trans" = "#0072B5", "N" = "#E18727", "P"="#20854E")))


gene_annotation = SupFig5B[SupFig5B[,2]%in%rownames(pt.matrix),c(2,3,4)]
tmp = gene_annotation[,1]
rownames(gene_annotation) = tmp
gene_annotation = gene_annotation[,-1]
colnames(gene_annotation) = c("Gene_classification", "MYC_target")

#### Ordering only
#gene_annotation = gene_annotation[c("DLL1", "DLL3", "DLL4", "FBXW7", "LNX1", "LFNG", "NOTCH1", "NOTCH2", "NOTCH3", "SOX9", "HEY2", "HEYL", "HES1", "MAMLD1", "JAG1"),]
gene_annotation = gene_annotation[c("NOTCH1", "NOTCH2", "NOTCH3", "SOX9", "HEY2", "HEYL", "HES1", "MAMLD1", "JAG1", "DLL1", "DLL3", "DLL4", "FBXW7", "LNX1", "LFNG"),]

pvalue_list = lapply(rownames(gene_annotation), function(x) {wilcox.test(as.numeric(SCLC_cpm2[x,SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "trans")]), as.numeric(SCLC_cpm2[x,SCLC_subtype_WTS_meta$SCLC_subtype=="N"]))$p.value}) %>% unlist %>% as.numeric() %>% round(.,3)
high_low = lapply(rownames(gene_annotation), function(x) {
                         if(median(as.numeric(SCLC_cpm2[x,SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "trans")])) >= median(as.numeric(SCLC_cpm2[x,SCLC_subtype_WTS_meta$SCLC_subtype=="N"]))){
                                 return("MYC-repressed")
                         }else{return("MYC-induced")}
  }) %>% unlist() %>% as.character()

gene_annotation = gene_annotation[pvalue_list<=0.5,]
high_low = high_low[pvalue_list<=0.5]
pvalue_list = pvalue_list[pvalue_list<=0.5]

gene_annotation$MYC_target = gsub("_", " ", gene_annotation$MYC_target)

rowAnn = rowAnnotation(df=gene_annotation,
                       annotation_height = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col=list(Gene_classification=c("Receptors" = "#812417", "Targets" = "#C8D0C5", "Ligands" = "#968179", "Transducers" = "#3E2F26", "Inhibitors" = "#00008b"),
           MYC_target = c("Non annotation" = "#ededff", "MYC target in RPM GEMM" = "#ff0000", "MYC induced" = "#9400D3", "MYC repressed" = "#00ff00")))

columne_split_NE = lapply(colnames(pt.matrix), function(x) {
                                  match_subtype = as.character(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$WTS_ID==x,"SCLC_subtype"]);
                                  return(match_subtype)
  }) %>% unlist() %>% as.character()
columne_split_NE = factor(columne_split_NE, levels = c("A", "trans", "N", "P"))
pt.matrix2 = data.frame(pt.matrix[,annotation$SCLC_subtype=="A"], pt.matrix[,annotation$SCLC_subtype=="trans"][,rev(order(pt.matrix["ASCL1",annotation$SCLC_subtype=="trans"]))], pt.matrix[,annotation$SCLC_subtype=="N"], pt.matrix[,annotation$SCLC_subtype=="P"])

set.seed(100)
plot_heatmap = Heatmap(pt.matrix2[rownames(gene_annotation),], name = paste0("Z-scale"), top_annotation = colAnn, cluster_columns = FALSE, right_annotation = rowAnn, cluster_rows = FALSE, show_column_names = FALSE, row_order = rownames(gene_annotation), row_split = high_low, column_split = columne_split_NE, row_gap = unit(5, "mm"), column_gap = unit(5, "mm"))

pdf(paste0(work_dir, "Fig4C_left_Heatmap_and_annotation_NE_and_P_subtype.pdf"), width = 10, height = 8)
print(plot_heatmap)
dev.off()


################################################################################

axis_theme = theme(axis.text.x = element_text(size=10, face="bold", color = "black"),
              axis.title.x = element_text(size=10, face="bold", color = "black"),
        axis.text.y = element_text(size=10, face="bold", color = "black"),
        axis.title.y = element_text(size=10, face="bold", color = "black"), panel.border = element_rect(size = 1))
my_comparisons <- list(c("A", "trans"), c("trans", "N"))




SCLC_cpm_Rank <- rankGenes(SCLC_cpm2)

SingScore <- simpleScore(
  rankData = SCLC_cpm_Rank,
  upSet = rownames(gene_annotation)[gene_annotation$MYC_target=="MYC induced"],
  downSet = rownames(gene_annotation)[gene_annotation$MYC_target=="MYC repressed"],
  centerScore = T,
  knownDirection = T
  )
SCLC_subtype_WTS_meta$Net_notch_score = SingScore$TotalScore

pvalue = wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "trans"),"Net_notch_score"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("N", "P"),"Net_notch_score"])$p.value %>% round(.,4)
#p = ggplot(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$NE_subtype=="NE",],
p = ggplot(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "trans", "N", "P"),],
           aes(x=SCLC_subtype, y=Net_notch_score, fill = SCLC_subtype)) + geom_boxplot(outlier.shape = NA) + ggtitle("upSet: MYC induced, downSet: MYC repressed")
p = p + ylab("Singscore") + xlab("") + ylim(-0.4, 0.7) + guides(fill="none")
p = p + theme_bw() + scale_fill_nejm() + axis_theme #+ stat_compare_means(comparisons = my_comparisons, method = "wilcox")
ggsave(paste0(work_dir, "Fig4C_right_Heatmap_Results_upSet_MYC_induced_downSet_MYC_repressed_singscore.png"), plot = p, width = 4, height = 4)

wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "trans"),"Net_notch_score"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("N"),"Net_notch_score"])$p.value
wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "trans"),"Net_notch_score"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("P"),"Net_notch_score"])$p.value






################################################
##### Fig 4D


P_Atrans_filt_gene = fgsea_P_Atrans[sign(fgsea_P_Atrans$NES) == sign(fgsea_N_Atrans$NES) * -1 & c(grepl("INTER", fgsea_P_Atrans$pathway) | grepl("ALLOGRAFT REJECTION", fgsea_P_Atrans$pathway)) , "leadingEdge"] %>% unlist() %>% as.character()


P_Atrans_filt_gene = P_Atrans_filt_gene[P_Atrans_filt_gene%in%rownames(P_Atrans_DEG)[P_Atrans_DEG$PValue<0.1]]




N_Atrans_filt_gene = fgsea_N_Atrans[sign(fgsea_N_Atrans$NES) == sign(fgsea_P_Atrans$NES) * -1 & c(grepl("INTER", fgsea_N_Atrans$pathway) | grepl("ALLOGRAFT REJECTION", fgsea_N_Atrans$pathway)) , "leadingEdge"] %>% unlist() %>% as.character()
N_Atrans_filt_gene = N_Atrans_filt_gene[N_Atrans_filt_gene%in%rownames(N_Atrans_DEG)[N_Atrans_DEG$PValue<0.1]]



Hallmark_ImmuneRelated_gene = P_Atrans_filt_gene[P_Atrans_filt_gene%in%N_Atrans_filt_gene] %>% unique()
Hallmark_ImmuneRelated_gene = Hallmark_ImmuneRelated_gene[Hallmark_ImmuneRelated_gene%in%rownames(SCLC_cpm)]
Hallmark_ImmuneRelated_gene = c(Hallmark_ImmuneRelated_gene, "HLA-DQA1", "HLA-C", "IFIH1", "IFIT3", "IL4R")




SCLC_subtype_WTS_meta_NE = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$NE_subtype=="NE" | SCLC_subtype_WTS_meta$SCLC_subtype == "P",]

pt.matrix = SCLC_cpm2[,SCLC_subtype_WTS_meta_NE$WTS_ID]
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- rownames(SCLC_cpm2); colnames(pt.matrix) = colnames(SCLC_cpm2)[SCLC_subtype_WTS_meta$NE_subtype=="NE" | SCLC_subtype_WTS_meta$SCLC_subtype == "P"]



tmp = SCLC_subtype_WTS_meta_NE
tmp = tmp[order(match(tmp$SCLC_subtype, c("A", "trans", "N", "P"))),]
pt.matrix = pt.matrix[,tmp$WTS_ID]
annotation = data.frame(tmp$SCLC_subtype, tmp$NE_25_genelist, tmp$nonNE_25_genelist)
rownames(annotation) = tmp$WTS_ID
colnames(annotation) = c("SCLC_subtype", "NE Signature", "non-NE Signature")
annotation = annotation[colnames(pt.matrix),]
annotation = data.frame(annotation)
colnames(annotation) = c("SCLC_subtype", "NE Signature", "non-NE Signature")
rownames(annotation) = colnames(pt.matrix)
annotation = annotation[,c("NE Signature", "non-NE Signature", "SCLC_subtype")]

col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white","red"))


annotation$SCLC_subtype = factor(annotation$SCLC_subtype, level=c("A", "trans", "N", "P"))
colAnn <- HeatmapAnnotation(df = annotation,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col = list('NE Signature' = col_fun,
             'non-NE Signature' = col_fun,
             SCLC_subtype = c("A" = "#BC3C29", "trans" = "#0072B5", "N" = "#E18727", "P"="#20854E")))




columne_split_NE = lapply(colnames(pt.matrix), function(x) {
                                  match_subtype = as.character(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$WTS_ID==x,"SCLC_subtype"]);
                                  return(match_subtype)
  }) %>% unlist() %>% as.character()
columne_split_NE = factor(columne_split_NE, levels = c("A", "trans", "N", "P"))
pt.matrix2 = data.frame(pt.matrix[,annotation$SCLC_subtype=="A"], pt.matrix[,annotation$SCLC_subtype=="trans"][,rev(order(pt.matrix["ASCL1",annotation$SCLC_subtype=="trans"]))], pt.matrix[,annotation$SCLC_subtype=="N"], pt.matrix[,annotation$SCLC_subtype=="P"])

gene_order = c(Hallmark_ImmuneRelated_gene[Hallmark_ImmuneRelated_gene%in%pathways.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE],
               Hallmark_ImmuneRelated_gene[Hallmark_ImmuneRelated_gene%in%pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE & !Hallmark_ImmuneRelated_gene%in%pathways.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE],
               Hallmark_ImmuneRelated_gene[Hallmark_ImmuneRelated_gene%in%pathways.hallmark$HALLMARK_ALLOGRAFT_REJECTION & !Hallmark_ImmuneRelated_gene%in%pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE])

print(sum(pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE =="IL4R"))

gene_annotation = c(rep("INTERFERON_ALPHA", length(Hallmark_ImmuneRelated_gene[Hallmark_ImmuneRelated_gene%in%pathways.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE])),
                    rep("INTERFERON_GAMMA", length(Hallmark_ImmuneRelated_gene[Hallmark_ImmuneRelated_gene%in%pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE & !Hallmark_ImmuneRelated_gene%in%pathways.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE])),
                    rep("ALLOGRAFT_REJECTION", length(Hallmark_ImmuneRelated_gene[Hallmark_ImmuneRelated_gene%in%pathways.hallmark$HALLMARK_ALLOGRAFT_REJECTION & !Hallmark_ImmuneRelated_gene%in%pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE])))
gene_annotation = data.frame(gene_annotation)
rownames(gene_annotation) = gene_order

rowAnn = rowAnnotation(df=gene_annotation,
                       annotation_height = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'), #annotation_text=anno_text(paste0(pvalue_list,"_",high_low, "     ")),
  col=list(Gene_classification=c("INTERFERON_ALPHA" = "#812417", "INTERFERON_GAMMA" = "#C8D0C5", "ALLOGRAFT_REJECTION" = "#00008b")))


pdf(paste0(work_dir, "Fig4D_left_Immune_related_heatmap_NE_and_P.pdf"), width = 11, height = 7.5)
Heatmap(pt.matrix2[rownames(gene_annotation),], name = paste0("Z-scale"), top_annotation = colAnn, cluster_columns = FALSE,  cluster_rows = FALSE, show_column_names = FALSE, column_split = columne_split_NE, row_gap = unit(5, "mm"), column_gap = unit(5, "mm"), row_order = rownames(gene_annotation), left_annotation = rowAnn)
dev.off()




SCLC_cpm_Rank <- rankGenes(SCLC_cpm2)

SingScore <- simpleScore(
  rankData = SCLC_cpm_Rank,
  upSet = Hallmark_ImmuneRelated_gene[Hallmark_ImmuneRelated_gene%in%rownames(SCLC_cpm2)],
  centerScore = T,
  knownDirection = T
  )
SCLC_subtype_WTS_meta$Hallmark_immune = SingScore$TotalScore



pvalue = wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "trans"),"Hallmark_immune"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("N", "P"),"Hallmark_immune"])$p.value %>% round(.,4)
p = ggplot(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "trans", "N", "P"),],
           aes(x=SCLC_subtype, y=Hallmark_immune, fill = SCLC_subtype)) + geom_boxplot(outlier.shape = NA)
p = p + ylab("Signature score") + xlab("") + guides(fill="none") + ylim(-0.3, 0.5)
p = p + theme_bw() + scale_fill_nejm() + axis_theme
ggsave(paste0(work_dir, "Fig4D_right_heatmap_Immune_related_signature_singscore.png"), plot = p, width = 4, height = 4)

wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "trans"),"Hallmark_immune"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("N"),"Hallmark_immune"])$p.value
wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "trans"),"Hallmark_immune"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("P"),"Hallmark_immune"])$p.value




###############
#####



SingScore <- simpleScore(
  rankData = SCLC_cpm_Rank,
  upSet = NE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$NE_25_genelist = SingScore$TotalScore

SingScore <- simpleScore(
  rankData = SCLC_cpm_Rank,
  upSet = nonNE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$nonNE_25_genelist = SingScore$TotalScore



gene_visualize_order = c(T_cell_Inflamed[!(T_cell_Inflamed%in%c("GZMA", "PRF1", nk))], "GZMA", "PRF1", nk[nk%in%T_cell_Inflamed], nk[!(nk%in%c("GZMA", "PRF1", T_cell_Inflamed))])
gene_visualize_order = gene_visualize_order[gene_visualize_order%in%rownames(pt.matrix)]

InflamedScore <- simpleScore(
  rankData = SCLC_cpm_Rank,
  upSet = gene_visualize_order,
  centerScore = T,
  knownDirection = T
  )
SCLC_subtype_WTS_meta$InflamedScore = InflamedScore$TotalScore


#(Y축) non-NE score vs. (X축) net-notch score
##########################################################################################
#################################### ASCL1 trans subtype
############################################################################


#(Y축) non-NE score vs. (X축) net-notch score

p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("A", "trans"), ], x = "Net_notch_score", y = "nonNE_25_genelist",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.3, label.y=0.2, label.sep = "\n"),
   xlim = c(-0.3, 0.4), ylim = c(-0.3, 0.2)
   )
ggsave(paste0(work_dir, "Fig_S5C_ASCL1_trans_subtype_x_Net_NOTCH_y_Non-NE_sig.png"), plot = p, width = 7, height = 7)


# (Y축) Hallmark immune score vs. (X축) non-NE score

p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("A", "trans"), ], x = "nonNE_25_genelist", y = "Hallmark_immune",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.1, label.y=0.3, label.sep = "\n"),
   xlim = c(-0.3, 0.2), ylim = c(-0.2, 0.3)
   )
ggsave(paste0(work_dir, "Fig_S5C_ASCL1_trans_subtype_x_Non-NE_y_Hallmark_immune.png"), plot = p, width = 7, height = 7)




# (Y축) Hallmark immune score vs. (X축) NetNOTCH score

p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("A", "trans"), ], x = "Net_notch_score", y = "Hallmark_immune",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.3, label.y=0.3, label.sep = "\n"),
   xlim = c(-0.3, 0.4), ylim = c(-0.2, 0.3)
   )
ggsave(paste0(work_dir, "Fig_S5C_ASCL1_trans_subtype_x_Net_NOTCH_y_Hallmark_immune.png"), plot = p, width = 7, height = 7)


###########################################################################################
#################################### NEUROD1 subtype
############################################################################




#(Y축) non-NE score vs. (X축) net-notch score

p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("N"), ], x = "Net_notch_score", y = "nonNE_25_genelist",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.3, label.y=0.2, label.sep = "\n"),
   xlim = c(-0.3, 0.4), ylim = c(-0.3, 0.2)
   )
ggsave(paste0(work_dir, "Fig_S5D_NEUROD1_subtype_x_Net_NOTCH_y_Non-NE_sig.png"), plot = p, width = 7, height = 7)


# (Y축) Hallmark immune score vs. (X축) non-NE score

p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("N"), ], x = "nonNE_25_genelist", y = "Hallmark_immune",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.1, label.y=0.3, label.sep = "\n"),
   xlim = c(-0.3, 0.2), ylim = c(-0.2, 0.3)
   )
ggsave(paste0(work_dir, "Fig_S5D_NEUROD1_subtype_x_Non-NE_y_Hallmark_immune.png"), plot = p, width = 7, height = 7)




# (Y축) Hallmark immune score vs. (X축) NetNOTCH score

p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("N"), ], x = "Net_notch_score", y = "Hallmark_immune",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.3, label.y=0.3, label.sep = "\n"),
   xlim = c(-0.3, 0.4), ylim = c(-0.2, 0.3)
   )
ggsave(paste0(work_dir, "Fig_S5D_NEUROD1_subtype_x_Net_NOTCH_y_Hallmark_immune.png"), plot = p, width = 7, height = 7)






###############################



gene_list = c(nk, T_cell_Inflamed) %>% unique()
#gene_list = T_cell_Inflamed
UnionScore <- simpleScore(
  rankData = SCLC_cpm_Rank,
  upSet = gene_list,
  centerScore = T,
  knownDirection = T
  )
SCLC_subtype_WTS_meta$UnionScore = UnionScore$TotalScore



##############
pt.matrix = SCLC_cpm2
pt.matrix <- apply(pt.matrix,2,function(x){(x-mean(x))/sd(x)})
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))

rownames(pt.matrix) <- rownames(SCLC_cpm2); colnames(pt.matrix) = colnames(SCLC_cpm2)


gene_visualize_order = c(T_cell_Inflamed[!(T_cell_Inflamed%in%c("GZMA", "PRF1", nk))], "GZMA", "PRF1", nk[nk%in%T_cell_Inflamed], nk[!(nk%in%c("GZMA", "PRF1", T_cell_Inflamed))])
gene_visualize_order = gene_visualize_order[gene_visualize_order%in%rownames(pt.matrix)]
gene_annotation = c(nk[!(nk%in%T_cell_Inflamed)], nk[nk%in%T_cell_Inflamed], T_cell_Inflamed[!(T_cell_Inflamed%in%nk)])

gene_annotation = data.frame(gene_annotation, lapply(gene_annotation, function(x) {
  if(sum(c("GZMA", "PRF1")%in%x)>=1){
          return("Cytolytic T cell")
  }else if(sum(nk[nk%in%T_cell_Inflamed]%in%x)>=1){
                       return("nk&TcellInflamed")
  }else if(sum(nk%in%x)>=1){
          return("nk")
  }else if(sum(T_cell_Inflamed%in%x)>=1){
          return("T cell Inflamed")
  }}) %>% unlist() %>% as.character()
)

gene_annotation = gene_annotation[gene_annotation[,1]%in%rownames(pt.matrix),]
tmp = gene_annotation[,1]
gene_annotation = data.frame(gene_annotation[,2])
rownames(gene_annotation) = tmp
colnames(gene_annotation) = c("Gene_Set")

annotation = data.frame(SCLC_subtype_WTS_meta$SCLC_subtype)
rownames(annotation) = SCLC_subtype_WTS_meta$WTS_ID
colnames(annotation) = c("SCLC_subtype")
annotation = annotation[colnames(pt.matrix),]
annotation = data.frame(annotation)
colnames(annotation) = c("SCLC_subtype")
rownames(annotation) = colnames(pt.matrix)


MFP = read.table("/BiO2/users/jjg/tools/MFP/test/results2//SCLC_test.txt")
MFP = MFP[rownames(annotation),]
annotation$MFP = factor(MFP, levels = c("IE", "IE/F", "F", "D"))


MFP = read.table("/BiO2/users/jjg/tools/MFP/test/results2//SCLC_score.txt")
MFP = MFP[rownames(annotation),]
annotation$Effector_cells = MFP$Effector_cells
annotation$Antitumor_cytokines = MFP$Antitumor_cytokines
annotation$NK_cells=MFP$NK_cells
annotation$T_cells=MFP$T_cells
annotation$B_cells=MFP$B_cells
annotation$Macrophages=MFP$Macrophages
annotation$Checkpoint_inhibition=MFP$Checkpoint_inhibition
annotation$Proliferation_rate=MFP$Proliferation_rate



col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white","red"))
viridis = colorRamp2(c(-2, 2), c("#fde725", "#440154"))

annotation = annotation[,c("MFP", "Effector_cells", "Antitumor_cytokines", "NK_cells", "T_cells", "B_cells", "Macrophages", "Checkpoint_inhibition", "Proliferation_rate", "SCLC_subtype")]

annotation= annotation[,c("MFP", "Proliferation_rate", "Checkpoint_inhibition", "B_cells", "Effector_cells", "NK_cells", "T_cells", "Antitumor_cytokines", "Macrophages", "SCLC_subtype")]
annotation$SCLC_subtype = factor(annotation$SCLC_subtype, level=c("A", "trans", "N", "P", "TN"))
colAnn <- HeatmapAnnotation(df = annotation,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col = list(SCLC_subtype = c("A" = "#BC3C29", "trans" = "#0072B5", "N" = "#E18727", "P" = "#20854E", "TN" = "#7876B1"),
             MFP = c("IE" = "#c24857", "IE/F" = "#cb7131", "F" = "#5cad68", "D" = "#204637"),
             Effector_cells=col_fun,
             Antitumor_cytokines=col_fun,
             NK_cells=col_fun,
             T_cells=col_fun,
             B_cells=col_fun,
             Macrophages=col_fun,
             Checkpoint_inhibition=col_fun,
             Proliferation_rate=col_fun
             ))

rowAnn = rowAnnotation(df=gene_annotation,
                       annotation_height = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col=list(Gene_Set=c("T cell Inflamed" = "#812417", "Cytolytic T cell" = "#C8D0C5", "nk&TcellInflamed" = "#968179", "nk" = "#3E2F26")))



################ start plotting
k_cluster = 2
i_distance_method = "maximum"

set.seed(100)
sample_scale_heatmap = Heatmap(pt.matrix[rownames(gene_annotation),], name = paste0("Z-scale"), top_annotation = colAnn, cluster_columns = TRUE, column_km = k_cluster, clustering_distance_columns = i_distance_method, left_annotation = rowAnn, cluster_rows = FALSE, show_column_names = FALSE, row_order = gene_visualize_order,
                               column_gap = unit(5, "mm"))
sample_scale_heatmap_ht = draw(sample_scale_heatmap); sample_scale_heatmap_tr_order = column_order(sample_scale_heatmap_ht)
#ggsave(paste0(work_dir, "sample_scale_", k_cluster, "_complexheatmap_", "pearson.png"), plot = sample_scale_heatmap_ht, width = 7, height = 7)
sample_scale_heatmap_ht = draw(sample_scale_heatmap); sample_scale_heatmap_tr_order = column_order(sample_scale_heatmap_ht)

sample_scale_heatmap_tr_order_name = names(sample_scale_heatmap_tr_order)
max_singscore_cluster = lapply(sample_scale_heatmap_tr_order_name, function(x) {
               sample_scale_heatmap_colnames = colnames(pt.matrix)[sample_scale_heatmap_tr_order[[x]]]
               return(mean(SCLC_subtype_WTS_meta$UnionScore[SCLC_subtype_WTS_meta$WTS_ID%in%sample_scale_heatmap_colnames]))
  }) %>% as.numeric() %>% which.max() %>% sample_scale_heatmap_tr_order_name[.]

surv_df_SCLC_subtype_WTS_meta = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$IO_line==1,]
tmp = rep("low_cluster", nrow(surv_df_SCLC_subtype_WTS_meta))
tmp[surv_df_SCLC_subtype_WTS_meta$WTS_ID%in%colnames(pt.matrix)[sample_scale_heatmap_tr_order[[max_singscore_cluster]]]] = "high_cluster"
surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster = tmp

pfs_fit <- survfit(Surv(IO_pfs, IO_pfs_event==1)~Immune_high_cluster, data = surv_df_SCLC_subtype_WTS_meta)#############################조심
median_pfs_fit = summary(pfs_fit)$table %>% data.frame()
median_pfs_fit = c(median_pfs_fit["Immune_high_cluster=high_cluster","median"], median_pfs_fit["Immune_high_cluster=low_cluster","median"])
pvalue = surv_pvalue(pfs_fit)$pval %>% round(., 5)


annotation = annotation[colnames(pt.matrix)[c(sample_scale_heatmap_tr_order[[1]], sample_scale_heatmap_tr_order[[2]])],]
colAnn <- HeatmapAnnotation(df = annotation,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col = list(SCLC_subtype = c("A" = "#BC3C29", "trans" = "#0072B5", "N" = "#E18727", "P" = "#20854E", "TN" = "#7876B1"),
             MFP = c("IE" = "#c24857", "IE/F" = "#cb7131", "F" = "#5cad68", "D" = "#204637"),
             Effector_cells=col_fun,
             Antitumor_cytokines=col_fun,
             NK_cells=col_fun,
             T_cells=col_fun,
             B_cells=col_fun,
             Macrophages=col_fun,
             Checkpoint_inhibition=col_fun,
             Proliferation_rate=col_fun
             ))

sample_scale_heatmap = Heatmap(pt.matrix[rownames(gene_annotation), c(sample_scale_heatmap_tr_order[[1]], sample_scale_heatmap_tr_order[[2]])],
                               column_split = c(rep("A", length(sample_scale_heatmap_tr_order[[1]])), rep("B", length(sample_scale_heatmap_tr_order[[2]]))),
                               name = paste0("Z-scale"), top_annotation = colAnn, cluster_columns = FALSE, clustering_distance_columns = i_distance_method, left_annotation = rowAnn, cluster_rows = FALSE, show_column_names = FALSE, row_order = gene_visualize_order,
                               column_gap = unit(5, "mm"))
SCLC_subtype_WTS_meta$Immune_cluster = rep("low_cluster", nrow(SCLC_subtype_WTS_meta))
SCLC_subtype_WTS_meta$Immune_cluster[SCLC_subtype_WTS_meta$WTS_ID %in% colnames(pt.matrix[,sample_scale_heatmap_tr_order[[2]]])] = "high_cluster"

pdf(paste0(work_dir, "Gene_Z_scale_", i_distance_method, "_K,", k_cluster, "_", pvalue, "_complexheatmap.pdf"), width = 14, height = 10)
print(sample_scale_heatmap)
dev.off()







surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster[surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster=="high_cluster"] = "Inflamed"
surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster[surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster=="low_cluster"] = "Non-Inflamed"
pfs_cox_fit <- coxph(Surv(IO_pfs, IO_pfs_event==1)~Immune_high_cluster, data = surv_df_SCLC_subtype_WTS_meta)
pfs_cox_fit_CI = exp(confint(pfs_cox_fit))
pfs_cox_fit = summary(pfs_cox_fit)

pfs_fit <- survfit(Surv(IO_pfs, IO_pfs_event==1)~Immune_high_cluster, data = surv_df_SCLC_subtype_WTS_meta)#############################조심

##################################
thickness = theme(legend.title = element_text(size = 13, color = "black", face = "bold"),
      legend.text = element_text(size = 13, color = "black", face = "bold"),
      axis.text.x = element_text(size = 13, color = "black", face = "bold"),
      axis.text.y = element_text(size = 13, color = "black", face = "bold"),
      axis.title.x = element_text(size = 13, color = "black", face = "bold"),
      axis.title.y = element_text(size = 13, color = "black", face = "bold"))
##################################


pdf(paste0(work_dir, "Gene_Z_scale_", i_distance_method,"_K,", k_cluster, "_", pvalue, "_cluster_sample_scale_matrix_Survplot.pdf"), width = 8, height = 6)
        p = ggsurvplot(
                fit = pfs_fit,
                size = 1.5,
                xlab = "Months since Treatment",
                ylab = "PFS",
                surv.median.line = "hv",
                palette = "nejm",
                xlim=c(0,13),
                ylim=c(0,1),
                break.time.by=3,
                #title=paste0("Immunotherapy PFS Analysis:\n",
        #                    "UnionGeneSet: T_cell_Inflamed&JTO_NK\n",
#                            "Cluster:", k_cluster, " ", i_distance_method, "\n",
#                            "High_median:", round(median_pfs_fit[1],1), "\n",
#                            "Low_median:", round(median_pfs_fit[2],1), "\n",
#                            "HR:", round(pfs_cox_fit$coefficients[,2], 1), ", 95%CI:", round(pfs_cox_fit_CI[,1],2),"-",round(pfs_cox_fit_CI[,2],2)),
                font.title=c("bold"), font.subtitle=c("italic"),
                legend.title="Strata",
                pval = T, pval.coord = c(1, 0.25),
                surv.scale = "percent",
                risk.table = TRUE, risk.table.height = 0.3, conf.int = F)
p$plot <- p$plot + thickness
print(p)
dev.off()




#########################################################
#########################################################

Responder_plot_df = data.frame("SCLC_subtype"=surv_df_SCLC_subtype_WTS_meta$SCLC_subtype, "IO_pfs_event"=surv_df_SCLC_subtype_WTS_meta$IO_pfs_event, "IO_pfs"=surv_df_SCLC_subtype_WTS_meta$IO_pfs, "Inflamed"=surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster, "WTS_ID"=surv_df_SCLC_subtype_WTS_meta$WTS_ID)
tmp = rep("DCB", nrow(Responder_plot_df))
#df$IO_pfs < 6 & df$IO_pfs_event == 1
tmp[(Responder_plot_df$IO_pfs<7 & surv_df_SCLC_subtype_WTS_meta$IO_pfs_event==1)] = "NDB"
Responder_plot_df$responder = tmp

Responder_plot_df$TME_subtype = lapply(Responder_plot_df$WTS_ID, function(x) {annotation[rownames(annotation)==x,"MFP"]}) %>% unlist() %>% as.character()
###################
data = lapply(c("DCB", "NDB"), function(x) {tmp = data.frame(table(Responder_plot_df$SCLC_subtype[Responder_plot_df$responder==x]));
                rownames(tmp) = tmp[,1]
                tmp=tmp[c("A", "trans", "N", "P", "TN"),2]}) %>% data.frame()
colnames(data) = c("DCB", "NDB")
rownames(data) = c("A", "trans", "N", "P", "TN")
data$SCLC_subtype = c("A", "trans", "N", "P", "TN")
data = data[c("A", "trans", "N"),]
fisher_pvalue = fisher.test(data[,c("DCB", "NDB")])
data = reshape2::melt(data)
data$SCLC_subtype = factor(data$SCLC_subtype, levels = c("A", "trans", "N"))

p=ggplot(data, aes(fill=SCLC_subtype, y=value, x=variable)) +
    geom_bar(position="fill", stat="identity") + theme_bw() + axis_theme + scale_fill_nejm() + scale_y_continuous(labels = function(x) paste0(x*100, "%"))
p = p + theme(axis.text.x = element_text(angle = 20, vjust = 0.5)) + xlab("") + ylab("") + ggtitle(paste0("Fisher p.value: ",round(fisher_pvalue$p.value, 3)))
ggsave(paste0(work_dir, "./responder/SCLC_subtype_proportion_by_responder.png"), plot=p, width = 4, height=4)


###################
data = lapply(c("DCB", "NDB"), function(x) {tmp = data.frame(table(Responder_plot_df$Inflamed[Responder_plot_df$responder==x]));
                rownames(tmp) = tmp[,1]
                tmp=tmp[c("Inflamed", "Non-Inflamed"),2]}) %>% data.frame()
colnames(data) = c("DCB", "NDB")
rownames(data) = c("Inflamed", "Non-Inflamed")
data$Inflamed = rownames(data)

fisher_pvalue = fisher.test(data[,c("DCB", "NDB")])

data = reshape2::melt(data)
data$Inflamed = factor(data$Inflamed, levels = c("Inflamed", "Non-Inflamed"))


p=ggplot(data, aes(fill=Inflamed, y=value, x=variable)) +
    geom_bar(position="fill", stat="identity") + theme_bw() + axis_theme + scale_fill_manual(values=c("#BC3C29", "#0072B5")) + scale_y_continuous(labels = function(x) paste0(x*100, "%"))
p = p + theme(axis.text.x = element_text(angle = 20, vjust = 0.5)) + xlab("") + ylab("") + ggtitle(paste0("Fisher p.value: ",round(fisher_pvalue$p.value, 3)))
ggsave(paste0(work_dir, "./responder/InflamedInfo_proportion_by_responder.png"), plot=p, width = 4, height=4)


###################
data = lapply(c("DCB", "NDB"), function(x) {tmp = data.frame(table(Responder_plot_df$TME_subtype[Responder_plot_df$responder==x]));
                rownames(tmp) = tmp[,1]
                rownames_ = c("IE", "IE/F", "F", "D")
                rownames_ = rownames_[rownames_%in%rownames(tmp)]
                tmp=tmp[rownames_,]})

data = merge(data[[1]], data[[2]], by = "Var1", all=TRUE)

colnames(data) = c("TME","DCB", "NDB")
data$DCB[is.na(data$DCB)] = 0
data$`NDB`[is.na(data$`NDB`)] = 0

fisher_pvalue = fisher.test(data[,c("DCB", "NDB")])

data = reshape2::melt(data)
data$TME = factor(data$TME, levels = c("IE", "IE/F", "F", "D"))


p=ggplot(data, aes(fill=TME, y=value, x=variable)) +
    geom_bar(position="fill", stat="identity") + theme_bw() + axis_theme + scale_fill_manual(values=c("#c24857", "#cb7131", "#5cad68", "#204637")) + scale_y_continuous(labels = function(x) paste0(x*100, "%"))
p = p + theme(axis.text.x = element_text(angle = 20, vjust = 0.5)) + xlab("") + ylab("") + ggtitle(paste0("Fisher p.value: ",round(fisher_pvalue$p.value, 3)))
ggsave(paste0(work_dir, "./responder/TME_subtype_proportion_by_responder.png"), plot=p, width = 3.5, height=4)


################## Fig 5E
MFP = read.table("/BiO2/users/jjg/tools/MFP/test/results2//SCLC_score.txt")
MFP = MFP[rownames(annotation),]
MFP = MFP[Responder_plot_df$WTS_ID,]

boxplot_df = data.frame(MFP, "responder"=Responder_plot_df$responder)
boxplot_df = reshape2::melt(boxplot_df)
boxplot_df$responder = factor(boxplot_df$responder, levels = c("DCB", "NDB"))
boxplot_df = boxplot_df[boxplot_df$variable %in% c("MHCI", "MHCII", "Effector_cells", "NK_cells", "T_cells", "CAF", "Matrix", "Matrix_remodeling", "Angiogenesis", "EMT_signature"),]
boxplot_df$variable = factor(boxplot_df$variable, levels = c("MHCI", "MHCII", "Effector_cells", "NK_cells", "T_cells", "CAF", "Matrix", "Matrix_remodeling", "Angiogenesis", "EMT_signature"))

p = ggplot(boxplot_df, aes(x = responder, y = value, fill=responder)) + geom_boxplot(outlier.shape = NA, width=0.7) + geom_point(position = position_jitter(seed = 42), alpha = 0.6) +  stat_compare_means(comparisons = list(c("DCB", "NDB"))) + facet_wrap(. ~ variable, ncol=5)
p = p + theme_bw() + axis_theme + scale_fill_manual(values=c("#BC3C29", "#204637")) + ylim(-2, 2.5)
ggsave(paste0(work_dir, "./responder/TME_subtype_signature_facet_boxplot.png"), plot=p, width = 12, height=6)





###############################
##### Supple Fig 6
################# 182 samples inflamed non-inflamed
library(edgeR)
library(org.Hs.eg.db)
library(dplyr); library(stringr); library(tibble); library(magrittr)


idfound <- Exp.count$gene_name %in% mappedRkeys(org.Hs.egSYMBOL)
length(idfound) #20345
Exp.count <- Exp.count[idfound,]
Exp.order <- order(rowSums(Exp.count %>% dplyr::select(-gene_name)), decreasing = TRUE)
Exp.count <- Exp.count[Exp.order,]
Exp.duplicated <- duplicated(Exp.count$gene_name)
Exp.count.removed <- Exp.count[Exp.duplicated,]
Exp.count <- Exp.count[!Exp.duplicated,]
library(magrittr); library(tibble)
Exp.count <- tibble::remove_rownames(Exp.count)
Exp.count %<>% tibble::column_to_rownames("gene_name")
head(Exp.count)




sample_group_info = SCLC_subtype_WTS_meta$Immune_cluster

Exp.count2 = Exp.count[,SCLC_subtype_WTS_meta$WTS_ID]
y <- DGEList(counts=Exp.count2,
             gene=rownames(Exp.count2),
             group=sample_group_info)
dim(y)
keep <- filterByExpr(y, group=sample_group_info)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")
y$samples
#eff.lib.size <- y$samples$lib.size*y$samples$norm.factors
Cpm.count <- edgeR::cpm(y, log=FALSE)
#boxplot(Cpm.count, col="gray", las=3)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

All_inflamed_versus_nonInflamed <- exactTest(y, pair=c("low_cluster", "high_cluster")) # compare groups 1 and 2
topTags(All_inflamed_versus_nonInflamed, n=10)


Inflamed_nonInflamed_DEG = All_inflamed_versus_nonInflamed$table
candidate_gene_list = c(pathway.kegg$KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION, pathway.kegg$KEGG_CHEMOKINE_SIGNALING_PATHWAY, pathway.kegg$KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY, pathway.kegg$KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY, pathway.kegg$KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION, pathways.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE, pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE, pathways.hallmark$HALLMARK_INFLAMMATORY_RESPONSE)

Inflamed_nonInflamed_DEG[rownames(Inflamed_nonInflamed_DEG) %in% candidate_gene_list & Inflamed_nonInflamed_DEG$logFC > 1 & Inflamed_nonInflamed_DEG$PValue < 0.05, ]

Immune_gene = c("CD74", "PTPRC", "ITGAL", "HLA-DMA", "BTK", "HLA-DQA1", "IL7R", "CCL3L3", "CD3D", "CD4", "IL4R", "CD79B", "CD274", "CD8A", "IL2RG","CD48", "CCL2", "IDO1")


Inflamed_nonInflamed_DEG$gene = rownames(Inflamed_nonInflamed_DEG)
keyvals = apply(Inflamed_nonInflamed_DEG, 1, function(x) {
                        if(as.numeric(x[3])< 0.05 & abs(as.numeric(x[1])) >= 1){
                                if(as.numeric(x[1]) >= 0.1){
                                        if(sum(x[4]%in%Immune_gene)>0){
                                                return('#800000')
                                        }else{return("#949494")}
                                }else if(as.numeric(x[1]) <= -0.1){
                                        return("#949494")
                                }
                        }else{return("#949494")}
                                              }) %>% unlist() %>% as.character()

names(keyvals)[keyvals == 'black'] <- 'A&trans NE'
#names(keyvals)[keyvals == '#c38e63'] <- 'nonhigh'# "#f6bdc0"
names(keyvals)[keyvals == '#800000'] <- 'Inflamed gene'
#names(keyvals)[keyvals == '#59788e'] <- 'nonlow' #"#b0dbf1"
names(keyvals)[keyvals == '#949494'] <- 'non-sig'


plot_Inflamed_nonInflamed_DEG = EnhancedVolcano(Inflamed_nonInflamed_DEG,
                lab = rownames(Inflamed_nonInflamed_DEG),
                 x = 'logFC', labSize=3,
                 y = 'PValue', ylim = c(0, 40),#xlim = c(-5.6, 5.6), ylim = c(0, 20),
                 pCutoff = 0.05,
                 pointSize = 1, selectLab=Immune_gene, colCustom=keyvals,
                 drawConnectors = TRUE,directionConnectors = "both",  max.overlaps=Inf, colAlpha = 0.9)

ggsave(paste0(work_dir, "182_samples_Inflamed_nonInflamed_VolcanoPlot.png"), plot = plot_Inflamed_nonInflamed_DEG, width = 5, height = 6)


Inflamed_nonInflamed_DEG = All_inflamed_versus_nonInflamed$table

Filt_DEG = Inflamed_nonInflamed_DEG[Inflamed_nonInflamed_DEG$PValue<=1,]

Filt_DEG$fcsign <- sign(Filt_DEG$logFC)
Filt_DEG$logP=-log10(Filt_DEG$PValue)
Filt_DEG$metric= Filt_DEG$logP/Filt_DEG$fcsign


Rank_DEG = data.frame(rownames(Filt_DEG), Filt_DEG$metric)
colnames(Rank_DEG) = c("Gene", "Rank")

fgseaRes <- fgsea(pathways = pathways.hallmark,
                  stats    = deframe(Rank_DEG),
                  minSize  = 15,
                  maxSize  = 500)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


fgseaRes$pathway = gsub("HALLMARK_", "", fgseaRes$pathway)
fgseaRes$pathway = gsub("_", " ", fgseaRes$pathway)
fgseaRes$log10pval = -log10(fgseaRes$pval)


p = ggplot(fgseaRes[grepl("NOTCH", fgseaRes$pathway) | c(!is.na(fgseaRes$pval) & abs(fgseaRes$NES)>=1 & fgseaRes$pval < 0.1),]) +
  geom_bar(stat = 'identity', aes(reorder(pathway, NES), NES, fill = log10pval)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0("Hallmark pathways NES from GSEA\n(Inflamed vs non-Inflamed) in 182samples")) +
  theme_bw() + scale_fill_gradient(low = "#ff0000", high = "#0000FF") + axis_theme + scale_fill_gradientn(limits = c(0,3), colours=c("blue", "red"))
ggsave(paste0(work_dir, "182_samples_Inflamed_nonInflamed_fgsea.png"), plot = p, widt=8.4, height=5)




########


MFP = read.table("/BiO2/users/jjg/tools/MFP/test/results2//SCLC_test.txt")
MFP = MFP[SCLC_subtype_WTS_meta$WTS_ID,]
SCLC_subtype_WTS_meta$TME_subtype = factor(MFP, levels = c("IE", "IE/F", "F", "D"))



surv_IO_SCLC_subtype_WTS_meta = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$IO_line == 1,]
surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype2 = as.character(surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype)
surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype2[surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype2 %in% c("A", "trans")] = "A&trans"
surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype2 = factor(surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype2, levels = c("A&trans", "N"))


pfs_fit <- survfit(Surv(IO_pfs, IO_pfs_event==1)~SCLC_subtype2, data = surv_IO_SCLC_subtype_WTS_meta)#############################조심
pvalue = surv_pvalue(pfs_fit)$pval %>% round(., 5)



pdf(paste0(work_dir, "pfs_by_A&trans_versus_N_Survplot_in_IO.pdf"), width = 8, height = 6)
        p = ggsurvplot(
                fit = pfs_fit,
                size = 1.5,
                xlab = "Months since Treatment",
                ylab = "PFS",
                surv.median.line = "hv",
                palette = "nejm",
                xlim=c(0,13),
                ylim=c(0,1),
                break.time.by=3,
                font.title=c("bold"), font.subtitle=c("italic"),
                legend.title="Strata",
                pval = T, pval.coord = c(1, 0.25),
                surv.scale = "percent",
                risk.table = TRUE, risk.table.height = 0.3, conf.int = F)
p$plot <- p$plot + thickness
print(p)
dev.off()





surv_IO_SCLC_subtype_WTS_meta$TME_subtype2 = as.character(surv_IO_SCLC_subtype_WTS_meta$TME_subtype)
surv_IO_SCLC_subtype_WTS_meta$TME_subtype2[surv_IO_SCLC_subtype_WTS_meta$TME_subtype2%in% c("IE/F","F")] = "IE/F & F"
surv_IO_SCLC_subtype_WTS_meta$TME_subtype2 = factor(surv_IO_SCLC_subtype_WTS_meta$TME_subtype2, levels = c("IE", "IE/F & F", "D"))

pfs_fit <- survfit(Surv(IO_pfs, IO_pfs_event==1)~TME_subtype2, data = surv_IO_SCLC_subtype_WTS_meta)#############################조심
pvalue = surv_pvalue(pfs_fit)$pval %>% round(., 5)



pdf(paste0(work_dir, "pfs_by_TME_subtype3_Survplot_in_IO.pdf"), width = 8, height = 6)
        p = ggsurvplot(
                fit = pfs_fit,
                size = 1.5,
                xlab = "Months since Treatment",
                ylab = "PFS",
                surv.median.line = "hv",
                palette = "nejm",
                xlim=c(0,13),
                ylim=c(0,1),
                break.time.by=3,
                font.title=c("bold"), font.subtitle=c("italic"),
                legend.title="Strata",
                pval = T, pval.coord = c(1, 0.25),
                surv.scale = "percent",
                risk.table = TRUE, risk.table.height = 0.3, conf.int = F)
p$plot <- p$plot + thickness
print(p)
dev.off()



p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("A", "trans"), ], x = "Net_notch_score", y = "InflamedScore",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.3, label.y=0.2, label.sep = "\n"),
   ylim = c(-0.4, 0.2)
   )
ggsave(paste0(work_dir, "Fig_S6E_A&trans_x_Net_NOTCH_y_InflamedScore.png"), plot = p, width = 7, height = 7)



#########################################################
#########################################################
paired_survdiff = data.frame()

i_percent = 0.67
library(cowplot)

i_quantile = quantile(SCLC_subtype_WTS_meta$Net_notch_score, i_percent)
SCLC_subtype_WTS_meta$NOTCH_high = lapply(SCLC_subtype_WTS_meta$Net_notch_score, function(x) {if(x>=i_quantile){return("High Notch")}else{return("Low Notch")}}) %>% unlist() %>% as.character()
SCLC_subtype_WTS_meta$NOTCH_high = factor(SCLC_subtype_WTS_meta$NOTCH_high, levels=c("Low Notch", "High Notch"))
#SCLC_subtype_WTS_meta$Immune_cluster = factor(SCLC_subtype_WTS_meta$Immune_cluster, levels=c("low_cluster", "high_cluster"))
SCLC_subtype_WTS_meta$Immune_cluster = as.character(SCLC_subtype_WTS_meta$Immune_cluster)
SCLC_subtype_WTS_meta$Immune_cluster[SCLC_subtype_WTS_meta$Immune_cluster=="high_cluster"] = "Inflamed"
SCLC_subtype_WTS_meta$Immune_cluster[SCLC_subtype_WTS_meta$Immune_cluster=="low_cluster"] = "Non-Inflamed"
SCLC_subtype_WTS_meta$Immune_cluster = factor(SCLC_subtype_WTS_meta$Immune_cluster, levels = c("Non-Inflamed", "Inflamed"))


##############################################################

thickness = theme(legend.title = element_text(size = 13, color = "black", face = "bold"),
      legend.text = element_text(size = 13, color = "black", face = "bold"),
      axis.text.x = element_text(size = 13, color = "black", face = "bold"),
      axis.text.y = element_text(size = 13, color = "black", face = "bold"),
      axis.title.x = element_text(size = 13, color = "black", face = "bold"),
      axis.title.y = element_text(size = 13, color = "black", face = "bold"))


two_two_table = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("A", "trans"),]

p1 = ggplot(data =  data.frame(table(two_two_table$NOTCH_high, two_two_table$Immune_cluster)), mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq), size=13, fontface="bold"), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") + ggtitle("A&trans: 144")

two_two_table = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("N"),]

p2 = ggplot(data =  data.frame(table(two_two_table$NOTCH_high, two_two_table$Immune_cluster)), mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq), size=13, fontface="bold"), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") + ggtitle("N: 19")


two_two_table = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("P"),]
p3 = ggplot(data =  data.frame(table(two_two_table$NOTCH_high, two_two_table$Immune_cluster)), mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq), size=13, fontface="bold"), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") + ggtitle("P: 15")

two_two_table = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("TN"),]
p4 = ggplot(data =  data.frame(table(two_two_table$NOTCH_high, two_two_table$Immune_cluster)), mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq), size=13, fontface="bold"), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") + ggtitle("TN type: 4")

p5 = ggplot(data =  data.frame(table(SCLC_subtype_WTS_meta$NOTCH_high, SCLC_subtype_WTS_meta$Immune_cluster)), mapping = aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq), size=13, fontface="bold"), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none") + ggtitle("All samples: 182")

pdf(paste0(work_dir, "top_", (1-i_percent),"_NOTCH_and_Immune_cluster_distribution.pdf"), width = 9, height = 12)
print(plot_grid(p1 + thickness, p2 + thickness ,p3 + thickness, p4 + thickness, p5 + thickness, ncol = 2))
dev.off()


