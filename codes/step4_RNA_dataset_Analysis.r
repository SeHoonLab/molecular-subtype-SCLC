LibPath = paste0(getwd(), "/LibPath/")
.libPaths(LibPath)
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

work_dir = "../Figures/step4/"
system(paste0("mkdir ", work_dir))

SCLC_meta = read.xlsx("../resource/Essential_check_Edited.xlsx", header = TRUE, sheetName = "Sheet1")
Signature_df = read.xlsx("../ref_signature/Signature gene list.xlsx", header = TRUE, sheetName="Sheet1")

SCLC_meta$SGI_ID = lapply(SCLC_meta$SGI_ID, function(x) {strsplit(gsub("^ ", "", x), " ")[[1]][1]}) %>% unlist() %>% as.character()
SCLC_subtype_WTS_meta = SCLC_meta[!is.na(SCLC_meta$SCLC_subtype) & c(SCLC_meta$WTS_QC_Result %in% c("Pass", "1")) & grepl("includ", SCLC_meta$study.inclusion),]
SCLC_subtype_WTS_meta$NE_subtype = lapply(SCLC_subtype_WTS_meta$SCLC_subtype, function(x) {
                                                           if(x=="A" | x=="N" | x=="AN"){
                                                                return("NE")
                                                           }else if(x=="P" | x=="TN"){
                                                                return("Non-NE")
                                                           }
}) %>% as.character()
SCLC_subtype_WTS_meta$SCLC_subtype = factor(SCLC_subtype_WTS_meta$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))
#system("wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.4/h.all.v7.4.symbols.gmt -P ../resource/")
pathways.hallmark <- gmtPathways("../resource/h.all.v7.4.symbols.gmt")
#system("wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.4/c2.cp.kegg.v7.4.symbols.gmt -P ../resource/")
pathway.kegg = gmtPathways("../resource/c2.cp.kegg.v7.4.symbols.gmt")

SCLC_TPM =read.table("../resource/log2_TPM_n226.txt", sep = "\t")

NE_25_genelist = strsplit(Signature_df[Signature_df$Signature.name=="Neuroendocrine",2], ", ")[[1]]
nonNE_25_genelist = strsplit(Signature_df[Signature_df$Signature.name=="Non-Neuroendocrine",2], ", ")[[1]]

nk <- strsplit(Signature_df[Signature_df$Signature.name=="NK signature",2], ", ")[[1]]
T_cell_Inflamed = strsplit(Signature_df[Signature_df$Signature.name=="T cell Inflamed",2], ", ")[[1]]

Exp.count = read.table("../resource/Exp.count.txt", header = TRUE, sep = '\t')

######################################################
################################# fsea plot
######################################################

set.seed(100)
N_AAN_DEG = readRDS("../Figures/step1/Fig_S5B_edgeR_N_AAN_DEG.Rds")
N_AAN_DEG = N_AAN_DEG$table

Filt_DEG = N_AAN_DEG[N_AAN_DEG$PValue<=1,]

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
fgsea_N_AAN = fgseaRes

p = ggplot(fgseaRes[grepl("NOTCH", fgseaRes$pathway) | (!is.na(fgseaRes$pval) & abs(fgseaRes$NES)>=1 & fgseaRes$pval < 0.1),]) +
  geom_bar(stat = 'identity', aes(reorder(pathway, NES), NES, fill = log10pval)) +
  coord_flip()  +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0("Hallmark pathways NES from GSEA(N vs A&AN)")) +
  theme_bw() + scale_fill_gradient(low = "#ff0000", high = "#0000FF") + axis_theme +  scale_fill_gradientn(limits = c(0,3), colours=c("blue", "red"))

#ggsave(paste0(work_dir, "Fig4B_N_AAN_DEG_pvalue_fgsea.png"), plot = p, width = 9, height=5)
pdf(paste0(work_dir, "Fig4B_N_AAN_DEG_pvalue_fgsea.pdf"),  width = 9, height=5)
print(p)
dev.off()


#P_vs A&AN

set.seed(100)
P_AAN_DEG = readRDS("../Figures/step1/Fig_S5A_edgeR_P_AAN_DEG.Rds")
P_AAN_DEG = P_AAN_DEG$table

Filt_DEG = P_AAN_DEG[P_AAN_DEG$PValue<=1,]

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
fgsea_P_AAN = fgseaRes


p = ggplot(fgseaRes[grepl("NOTCH", fgseaRes$pathway) | c(!is.na(fgseaRes$pval) & abs(fgseaRes$NES)>=1 & fgseaRes$pval < 0.1) & !(fgseaRes$pathway%in% c("KRAS SIGNALING DN", "ESTROGEN RESPONSE EARLY")),]) +
  geom_bar(stat = 'identity', aes(reorder(pathway, NES), NES, fill = log10pval)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0("Hallmark pathways NES from GSEA(P vs A&AN)")) +
  theme_bw() + scale_fill_gradient(low = "#ff0000", high = "#0000FF") + axis_theme + scale_fill_gradientn(limits = c(0,3), colours=c("blue", "red"))

#ggsave(paste0(work_dir, "Fig4A_P_A&AN_DEG_pvalue_fgsea.png"), plot = p, widt=8.4, height=5)

pdf(paste0(work_dir, "Fig4A_P_A&AN_DEG_pvalue_fgsea.pdf"), widt=8.4, height=5)
print(p)
dev.off()


#########


SCLC_TPM2 = SCLC_TPM[,SCLC_subtype_WTS_meta$WTS_ID]




SCLC_TPM_Rank <- rankGenes(SCLC_TPM2)


SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = NE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$NE_25_genelist = rescale(SingScore$TotalScore, to = c(-4, 4))

SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = nonNE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$nonNE_25_genelist = rescale(SingScore$TotalScore, to = c(-4, 4))


################################################
##### Fig 4C



SCLC_subtype_WTS_meta_NE_P = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$NE_subtype=="NE" | SCLC_subtype_WTS_meta$SCLC_subtype == "P",]

pt.matrix = SCLC_TPM2[,SCLC_subtype_WTS_meta_NE_P$WTS_ID]
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- rownames(SCLC_TPM2); colnames(pt.matrix) = colnames(SCLC_TPM2)[SCLC_subtype_WTS_meta$NE_subtype=="NE" | SCLC_subtype_WTS_meta$SCLC_subtype == "P"]



tmp = SCLC_subtype_WTS_meta_NE_P
tmp = tmp[order(match(tmp$SCLC_subtype, c("A", "AN", "N", "P"))),]
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


annotation$SCLC_subtype = factor(annotation$SCLC_subtype, level=c("A", "AN", "N", "P"))
colAnn <- HeatmapAnnotation(df = annotation,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col = list('NE Signature' = col_fun,
             'non-NE Signature' = col_fun,
             SCLC_subtype = c("A" = "#BC3C29", "AN" = "#0072B5", "N" = "#E18727", "P"="#20854E")))

gene_annotation = apply(Signature_df[Signature_df$Signature.name=="Myc related Notch",c(2,3)], 2, function(x) {strsplit(x, ", ")[[1]]}) %>% data.frame()
gene_annotation = gene_annotation[gene_annotation[,1]%in%rownames(pt.matrix),]
rownames(gene_annotation) = gene_annotation[,1]
colnames(gene_annotation) = c("GeneName","MYC_target")

#### Ordering only
annotation_order = c("NOTCH1", "NOTCH2", "NOTCH3", "SOX9", "HEY2", "HEYL", "HES1", "MAMLD1", "JAG1", "DLL1", "DLL3", "DLL4", "FBXW7", "LNX1", "LFNG")
gene_annotation = gene_annotation[annotation_order,]

pvalue_list = lapply(rownames(gene_annotation), function(x) {wilcox.test(as.numeric(SCLC_TPM2[x,SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN")]), as.numeric(SCLC_TPM2[x,SCLC_subtype_WTS_meta$SCLC_subtype=="N"]))$p.value}) %>% unlist %>% as.numeric() %>% round(.,3)
high_low = lapply(rownames(gene_annotation), function(x) {
                         if(mean(as.numeric(SCLC_TPM2[x,SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN")])) >= mean(as.numeric(SCLC_TPM2[x,SCLC_subtype_WTS_meta$SCLC_subtype=="N"]))){
                                 return("MYC-repressed")
                         }else{return("MYC-induced")}
  }) %>% unlist() %>% as.character()

gene_annotation = gene_annotation[pvalue_list<=0.5,]
high_low = high_low[pvalue_list<=0.5]
pvalue_list = pvalue_list[pvalue_list<=0.5]

gene_annotation = data.frame(gene_annotation[,2])
rownames(gene_annotation) = annotation_order
colnames(gene_annotation) = c("MYC_target")



gene_annotation$MYC_target = gsub("_", " ", gene_annotation$MYC_target)

rowAnn = rowAnnotation(df=gene_annotation,
                       annotation_height = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col=list(MYC_target = c("Non annotation" = "#ededff", "MYC target in RPM GEMM" = "#ff0000", "MYC induced" = "#9400D3", "MYC repressed" = "#00ff00")))

columne_split_NE = lapply(colnames(pt.matrix), function(x) {
                                  match_subtype = as.character(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$WTS_ID==x,"SCLC_subtype"]);
                                  return(match_subtype)
  }) %>% unlist() %>% as.character()
columne_split_NE = factor(columne_split_NE, levels = c("A", "AN", "N", "P"))
pt.matrix2 = data.frame(pt.matrix[,annotation$SCLC_subtype=="A"], pt.matrix[,annotation$SCLC_subtype=="AN"][,rev(order(pt.matrix["ASCL1",annotation$SCLC_subtype=="AN"]))], pt.matrix[,annotation$SCLC_subtype=="N"], pt.matrix[,annotation$SCLC_subtype=="P"])

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
my_comparisons <- list(c("A", "AN"), c("AN", "N"))




SCLC_TPM_Rank <- rankGenes(SCLC_TPM2)

SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = rownames(gene_annotation)[gene_annotation$MYC_target=="MYC induced"],
  downSet = rownames(gene_annotation)[gene_annotation$MYC_target=="MYC repressed"],
  centerScore = T,
  knownDirection = T
  )
SCLC_subtype_WTS_meta$Net_notch_score = SingScore$TotalScore

pvalue = wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN"),"Net_notch_score"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("N", "P"),"Net_notch_score"])$p.value %>% round(.,4)
p = ggplot(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN", "N", "P"),],
           aes(x=SCLC_subtype, y=Net_notch_score, fill = SCLC_subtype)) + geom_boxplot(outlier.shape = NA) + ggtitle("upSet: MYC induced, downSet: MYC repressed")
p = p + ylab("Singscore") + xlab("") + ylim(-0.4, 0.7) + guides(fill="none")
p = p + theme_bw() + scale_fill_nejm() + axis_theme 
#ggsave(paste0(work_dir, "Fig4C_right_Heatmap_Results_upSet_MYC_induced_downSet_MYC_repressed_singscore.png"), plot = p, width = 4, height = 4)

wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN"),"Net_notch_score"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("N"),"Net_notch_score"])$p.value
wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN"),"Net_notch_score"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("P"),"Net_notch_score"])$p.value

pdf(paste0(work_dir, "Fig4C_right_Heatmap_Results_upSet_MYC_induced_downSet_MYC_repressed_singscore.pdf"), width = 4, height = 4)
print(p)# + ggtitle(paste0("A&AN vs N : ", wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN"),"Net_notch_score"],
#                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("N"),"Net_notch_score"])$p.value,"\n",
#			 "A&AN vs P : ", wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN"),"Net_notch_score"],
#                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("P"),"Net_notch_score"])$p.value)))
dev.off()




################################################
##### Fig 4D


P_AAN_filt_gene = fgsea_P_AAN[sign(fgsea_P_AAN$NES) == sign(fgsea_N_AAN$NES) * -1 & c(grepl("INTER", fgsea_P_AAN$pathway) | grepl("ALLOGRAFT REJECTION", fgsea_P_AAN$pathway)) , "leadingEdge"] %>% unlist() %>% as.character()


P_AAN_filt_gene = P_AAN_filt_gene[P_AAN_filt_gene%in%rownames(P_AAN_DEG)[P_AAN_DEG$PValue<0.1]]




N_AAN_filt_gene = fgsea_N_AAN[sign(fgsea_N_AAN$NES) == sign(fgsea_P_AAN$NES) * -1 & c(grepl("INTER", fgsea_N_AAN$pathway) | grepl("ALLOGRAFT REJECTION", fgsea_N_AAN$pathway)) , "leadingEdge"] %>% unlist() %>% as.character()
N_AAN_filt_gene = N_AAN_filt_gene[N_AAN_filt_gene%in%rownames(N_AAN_DEG)[N_AAN_DEG$PValue<0.1]]



Hallmark_ImmuneRelated_gene = P_AAN_filt_gene[P_AAN_filt_gene%in%N_AAN_filt_gene] %>% unique()
Hallmark_ImmuneRelated_gene = Hallmark_ImmuneRelated_gene[Hallmark_ImmuneRelated_gene%in%rownames(SCLC_TPM)]
Hallmark_ImmuneRelated_gene = c(Hallmark_ImmuneRelated_gene, "HLA-DQA1", "HLA-C", "IFIH1", "IFIT3", "IL4R")




SCLC_subtype_WTS_meta_NE = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$NE_subtype=="NE" | SCLC_subtype_WTS_meta$SCLC_subtype == "P",]

pt.matrix = SCLC_TPM2[,SCLC_subtype_WTS_meta_NE$WTS_ID]
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- rownames(SCLC_TPM2); colnames(pt.matrix) = colnames(SCLC_TPM2)[SCLC_subtype_WTS_meta$NE_subtype=="NE" | SCLC_subtype_WTS_meta$SCLC_subtype == "P"]



tmp = SCLC_subtype_WTS_meta_NE
tmp = tmp[order(match(tmp$SCLC_subtype, c("A", "AN", "N", "P"))),]
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


annotation$SCLC_subtype = factor(annotation$SCLC_subtype, level=c("A", "AN", "N", "P"))
colAnn <- HeatmapAnnotation(df = annotation,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col = list('NE Signature' = col_fun,
             'non-NE Signature' = col_fun,
             SCLC_subtype = c("A" = "#BC3C29", "AN" = "#0072B5", "N" = "#E18727", "P"="#20854E")))




columne_split_NE = lapply(colnames(pt.matrix), function(x) {
                                  match_subtype = as.character(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$WTS_ID==x,"SCLC_subtype"]);
                                  return(match_subtype)
  }) %>% unlist() %>% as.character()
columne_split_NE = factor(columne_split_NE, levels = c("A", "AN", "N", "P"))
pt.matrix2 = data.frame(pt.matrix[,annotation$SCLC_subtype=="A"], pt.matrix[,annotation$SCLC_subtype=="AN"][,rev(order(pt.matrix["ASCL1",annotation$SCLC_subtype=="AN"]))], pt.matrix[,annotation$SCLC_subtype=="N"], pt.matrix[,annotation$SCLC_subtype=="P"])

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
  gap = unit(1, 'mm'), 
  col=list(Gene_classification=c("INTERFERON_ALPHA" = "#812417", "INTERFERON_GAMMA" = "#C8D0C5", "ALLOGRAFT_REJECTION" = "#00008b")))


pdf(paste0(work_dir, "Fig4D_left_Immune_related_heatmap_NE_and_P.pdf"), width = 11, height = 7.5)
Heatmap(pt.matrix2[rownames(gene_annotation),], name = paste0("Z-scale"), top_annotation = colAnn, cluster_columns = FALSE,  cluster_rows = FALSE, show_column_names = FALSE, column_split = columne_split_NE, row_gap = unit(5, "mm"), column_gap = unit(5, "mm"), row_order = rownames(gene_annotation), left_annotation = rowAnn)
dev.off()




SCLC_TPM_Rank <- rankGenes(SCLC_TPM2)

SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = Hallmark_ImmuneRelated_gene[Hallmark_ImmuneRelated_gene%in%rownames(SCLC_TPM2)],
  centerScore = T,
  knownDirection = T
  )
SCLC_subtype_WTS_meta$Hallmark_immune = SingScore$TotalScore



pvalue = wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN"),"Hallmark_immune"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("N", "P"),"Hallmark_immune"])$p.value %>% round(.,4)
p = ggplot(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN", "N", "P"),],
           aes(x=SCLC_subtype, y=Hallmark_immune, fill = SCLC_subtype)) + geom_boxplot(outlier.shape = NA)
p = p + ylab("Signature score") + xlab("") + guides(fill="none") + ylim(-0.3, 0.5)
p = p + theme_bw() + scale_fill_nejm() + axis_theme
#ggsave(paste0(work_dir, "Fig4D_right_heatmap_Immune_related_signature_singscore.png"), plot = p, width = 4, height = 4)

wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN"),"Hallmark_immune"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("N"),"Hallmark_immune"])$p.value
wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN"),"Hallmark_immune"],
                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("P"),"Hallmark_immune"])$p.value

pdf(paste0(work_dir, "Fig4D_right_heatmap_Immune_related_signature_singscore.pdf"), width = 4, height = 4)
print(p)# + ggtitle(paste0("A&AN vs N : ", wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN"),"Hallmark_immune"],
#                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("N"),"Hallmark_immune"])$p.value,"\n",
#			 "A&AN vs P : ", wilcox.test(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN"),"Hallmark_immune"],
#                     SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("P"),"Hallmark_immune"])$p.value)))
dev.off()


###############
#####



SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = NE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$NE_25_genelist = SingScore$TotalScore

SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = nonNE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$nonNE_25_genelist = SingScore$TotalScore


# Add ordering rule
gene_visualize_order = c(T_cell_Inflamed[!(T_cell_Inflamed%in%c("GZMA", "PRF1", nk))], "GZMA", "PRF1", nk[nk%in%T_cell_Inflamed], nk[!(nk%in%c("GZMA", "PRF1", T_cell_Inflamed))])
gene_visualize_order = gene_visualize_order[gene_visualize_order%in%rownames(pt.matrix)]

InflamedScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = gene_visualize_order,
  centerScore = T,
  knownDirection = T
  )
SCLC_subtype_WTS_meta$InflamedScore = InflamedScore$TotalScore

###############################

SCLC_TPM_Rank <- rankGenes(SCLC_TPM2)

gene_list = c(nk, T_cell_Inflamed) %>% unique()
#gene_list = T_cell_Inflamed
UnionScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = gene_list,
  centerScore = T,
  knownDirection = T
  )
SCLC_subtype_WTS_meta$UnionScore = UnionScore$TotalScore



##############
pt.matrix = SCLC_TPM2
pt.matrix <- apply(pt.matrix,2,function(x){(x-mean(x))/sd(x)})
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))

rownames(pt.matrix) <- rownames(SCLC_TPM2); colnames(pt.matrix) = colnames(SCLC_TPM2)


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

MFP = read.table("../resource/MFP_results/SCLC_TME_annotation.txt")

MFP = MFP[rownames(annotation),]
annotation$MFP = factor(MFP, levels = c("IE", "IE/F", "F", "D"))



col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white","red"))
viridis = colorRamp2(c(-2, 2), c("#fde725", "#440154"))

annotation = annotation[,c("MFP", "SCLC_subtype")]

annotation= annotation[,c("MFP", "SCLC_subtype")]
annotation$SCLC_subtype = factor(annotation$SCLC_subtype, level=c("A", "AN", "N", "P", "TN"))
colAnn <- HeatmapAnnotation(df = annotation,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col = list(SCLC_subtype = c("A" = "#BC3C29", "AN" = "#0072B5", "N" = "#E18727", "P" = "#20854E", "TN" = "#7876B1"),
             MFP = c("IE" = "#c24857", "IE/F" = "#cb7131", "F" = "#5cad68", "D" = "#204637")
             ))

rowAnn = rowAnnotation(df=gene_annotation,
                       annotation_height = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col=list(Gene_Set=c("T cell Inflamed" = "#812417", "Cytolytic T cell" = "#C8D0C5", "nk&TcellInflamed" = "#968179", "nk" = "#3E2F26")))



################ start plotting
k_cluster = 2
i_distance_method = "canberra"; i_cluster_medhod = "average"


set.seed(100)
sample_scale_heatmap = Heatmap(pt.matrix[rownames(gene_annotation),], name = paste0("Z-scale"), top_annotation = colAnn, cluster_columns = TRUE, column_split = k_cluster, clustering_distance_columns = i_distance_method, left_annotation = rowAnn, cluster_rows = FALSE, show_column_names = FALSE, row_order = gene_visualize_order,
                               column_gap = unit(5, "mm"), clustering_method_columns = i_cluster_medhod)
sample_scale_heatmap_ht = draw(sample_scale_heatmap); sample_scale_heatmap_tr_order = column_order(sample_scale_heatmap_ht)
#ggsave(paste0(work_dir, "sample_scale_", k_cluster, "_complexheatmap_", "pearson.png"), plot = sample_scale_heatmap_ht, width = 7, height = 7)
sample_scale_heatmap_ht = draw(sample_scale_heatmap); sample_scale_heatmap_tr_order = column_order(sample_scale_heatmap_ht)

pdf(paste0(work_dir,"TME_annotated_inflamed_clustering_complexheatmap1.pdf"), width = 14, height = 10)
print(sample_scale_heatmap_ht)
dev.off()

sample_scale_heatmap_tr_order_name = seq(1, length(sample_scale_heatmap_tr_order))
max_singscore_cluster = lapply(sample_scale_heatmap_tr_order_name, function(x) {
               sample_scale_heatmap_colnames = colnames(pt.matrix)[sample_scale_heatmap_tr_order[[x]]]
               return(mean(SCLC_subtype_WTS_meta$UnionScore[SCLC_subtype_WTS_meta$WTS_ID%in%sample_scale_heatmap_colnames]))
  }) %>% as.numeric() %>% which.max() %>% sample_scale_heatmap_tr_order_name[.]

surv_df_SCLC_subtype_WTS_meta = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$IO_line==1,]
tmp = rep("low_cluster", nrow(surv_df_SCLC_subtype_WTS_meta))
tmp[surv_df_SCLC_subtype_WTS_meta$WTS_ID%in%colnames(pt.matrix)[sample_scale_heatmap_tr_order[[max_singscore_cluster]]]] = "high_cluster"
surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster = tmp

pfs_fit <- survfit(Surv(IO_pfs, IO_pfs_event==1)~Immune_high_cluster, data = surv_df_SCLC_subtype_WTS_meta)#############################조심
cox_fit = coxph(Surv(ctx_1st_pfs, ctx_1st_pfs_event==1)~Immune_high_cluster, data = surv_df_SCLC_subtype_WTS_meta)
median_pfs_fit = summary(pfs_fit)$table %>% data.frame()
median_pfs_fit = c(median_pfs_fit["Immune_high_cluster=high_cluster","median"], median_pfs_fit["Immune_high_cluster=low_cluster","median"])
pvalue = surv_pvalue(pfs_fit)$pval %>% round(., 5)

if(pvalue >= 0.05){next}

annotation = annotation[colnames(pt.matrix)[c(sample_scale_heatmap_tr_order[[1]], sample_scale_heatmap_tr_order[[2]])],]
colAnn <- HeatmapAnnotation(df = annotation,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm'),
  col = list(SCLC_subtype = c("A" = "#BC3C29", "AN" = "#0072B5", "N" = "#E18727", "P" = "#20854E", "TN" = "#7876B1"),
             MFP = c("IE" = "#c24857", "IE/F" = "#cb7131", "F" = "#5cad68", "D" = "#204637")
             ))

sample_scale_heatmap = Heatmap(pt.matrix[rownames(gene_annotation), c(sample_scale_heatmap_tr_order[[1]], sample_scale_heatmap_tr_order[[2]])],
                               column_split = c(rep("A", length(sample_scale_heatmap_tr_order[[1]])), rep("B", length(sample_scale_heatmap_tr_order[[2]]))),
                               name = paste0("Z-scale"), top_annotation = colAnn, cluster_columns = FALSE, clustering_distance_columns = i_distance_method, left_annotation = rowAnn, cluster_rows = FALSE, show_column_names = FALSE, row_order = gene_visualize_order,
                               column_gap = unit(5, "mm"))
SCLC_subtype_WTS_meta$Immune_cluster = rep("low_cluster", nrow(SCLC_subtype_WTS_meta))
SCLC_subtype_WTS_meta$Immune_cluster[SCLC_subtype_WTS_meta$WTS_ID %in% colnames(pt.matrix[,sample_scale_heatmap_tr_order[[max_singscore_cluster]]])] = "high_cluster"

N_count = table(SCLC_subtype_WTS_meta$Immune_cluster, SCLC_subtype_WTS_meta$SCLC_subtype)
N_count = N_count["high_cluster", "N"]
if(N_count>=5){next}


P_count = table(SCLC_subtype_WTS_meta$Immune_cluster, SCLC_subtype_WTS_meta$SCLC_subtype)
P_count = P_count["high_cluster", "P"]

pdf(paste0(work_dir, pvalue,"_", N_count,"_", round(N_count/19,2),"_", round(P_count/15, 2), "_", i_distance_method,"_", i_cluster_medhod, "_",
	  round(table(SCLC_subtype_WTS_meta$Immune_cluster)["high_cluster"]/nrow(SCLC_subtype_WTS_meta), 2),
	  "_TME_annotated_inflamed_clustering_complexheatmap.pdf"), width = 14, height = 10)
print(sample_scale_heatmap)
dev.off()



surv_df_SCLC_subtype_WTS_meta$SCLC_subtype2 = as.character(surv_df_SCLC_subtype_WTS_meta$SCLC_subtype)
surv_df_SCLC_subtype_WTS_meta$SCLC_subtype2[surv_df_SCLC_subtype_WTS_meta$SCLC_subtype2 %in% c("A", "AN")] = "A&AN"
pfs_fit <- survfit(Surv(IO_pfs, IO_pfs_event==1)~SCLC_subtype2, data = surv_df_SCLC_subtype_WTS_meta)#############################조심
cox_fit <- coxph(Surv(IO_pfs, IO_pfs_event==1)~SCLC_subtype2, data = surv_df_SCLC_subtype_WTS_meta)




surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster[surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster=="high_cluster"] = "Inflamed"
surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster[surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster=="low_cluster"] = "Non-Inflamed"


pfs_fit <- survfit(Surv(IO_pfs, IO_pfs_event==1)~Immune_high_cluster, data = surv_df_SCLC_subtype_WTS_meta)#############################조심
cox_fit <- coxph(Surv(IO_pfs, IO_pfs_event==1)~Immune_high_cluster, data = surv_df_SCLC_subtype_WTS_meta)
cox_summary = cbind(data.frame(summary(cox_fit)$coefficients), data.frame(summary(cox_fit)$conf.int))

##################################
thickness = theme(legend.title = element_text(size = 13, color = "black", face = "bold"),
      legend.text = element_text(size = 13, color = "black", face = "bold"),
      axis.text.x = element_text(size = 13, color = "black", face = "bold"),
      axis.text.y = element_text(size = 13, color = "black", face = "bold"),
      axis.title.x = element_text(size = 13, color = "black", face = "bold"),
      axis.title.y = element_text(size = 13, color = "black", face = "bold"))
##################################


pdf(paste0(work_dir, "Fig5C_", pvalue,"_", N_count, "_", i_distance_method,"_", i_cluster_medhod,"_first_line_IO_sample_Survplot_by_inflamed.pdf"), width = 8, height = 6)
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

###################

Responder_plot_df = data.frame("SCLC_subtype"=surv_df_SCLC_subtype_WTS_meta$SCLC_subtype, "IO_pfs_event"=surv_df_SCLC_subtype_WTS_meta$IO_pfs_event, "IO_pfs"=surv_df_SCLC_subtype_WTS_meta$IO_pfs, "Inflamed"=surv_df_SCLC_subtype_WTS_meta$Immune_high_cluster, "WTS_ID"=surv_df_SCLC_subtype_WTS_meta$WTS_ID)
tmp = rep("DCB", nrow(Responder_plot_df))
#df$IO_pfs < 6 & df$IO_pfs_event == 1
tmp[(Responder_plot_df$IO_pfs<7 & surv_df_SCLC_subtype_WTS_meta$IO_pfs_event==1)] = "NDB"
Responder_plot_df$responder = tmp

Responder_plot_df$TME_subtype = lapply(Responder_plot_df$WTS_ID, function(x) {annotation[rownames(annotation)==x,"MFP"]}) %>% unlist() %>% as.character()

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
pdf(paste0(work_dir, "Fig5D_", pvalue,"_", N_count, "_", i_distance_method,"_", i_cluster_medhod,"_InflamedInfo_proportion_by_responder.pdf"), width = 4, height = 4)
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
                tmp=tmp[c("A", "AN", "N", "P", "TN"),2]}) %>% data.frame()
colnames(data) = c("DCB", "NDB")
rownames(data) = c("A", "AN", "N", "P", "TN")
data$SCLC_subtype = c("A", "AN", "N", "P", "TN")
data = data[c("A", "AN", "N"),]
fisher_pvalue = fisher.test(data[,c("DCB", "NDB")])
data = reshape2::melt(data)
data$SCLC_subtype = factor(data$SCLC_subtype, levels = c("A", "AN", "N"))

p=ggplot(data, aes(fill=SCLC_subtype, y=value, x=variable)) +
    geom_bar(position="fill", stat="identity") + theme_bw() + axis_theme + scale_fill_nejm() + scale_y_continuous(labels = function(x) paste0(x*100, "%"))
p = p + theme(axis.text.x = element_text(angle = 20, vjust = 0.5)) + xlab("") + ylab("") + ggtitle(paste0("Fisher p.value: ",round(fisher_pvalue$p.value, 3)))
ggsave(paste0(work_dir, "Fig5D_SCLC_subtype_proportion_by_responder.png"), plot=p, width = 4, height=4)


###################

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
ggsave(paste0(work_dir, "Fig_S6D_TME_subtype_proportion_by_responder.png"), plot=p, width = 3.5, height=4)



###############################
MFP = read.table("../resource/MFP_results/SCLC_TME_annotation.txt")

MFP = MFP[SCLC_subtype_WTS_meta$WTS_ID,]
SCLC_subtype_WTS_meta$TME_subtype = factor(MFP, levels = c("IE", "IE/F", "F", "D"))



surv_IO_SCLC_subtype_WTS_meta = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$IO_line == 1,]
surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype2 = as.character(surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype)
surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype2[surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype2 %in% c("A", "AN")] = "A&AN"
surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype2 = factor(surv_IO_SCLC_subtype_WTS_meta$SCLC_subtype2, levels = c("A&AN", "N"))


pfs_fit <- survfit(Surv(IO_pfs, IO_pfs_event==1)~SCLC_subtype2, data = surv_IO_SCLC_subtype_WTS_meta)#############################조심
pvalue = surv_pvalue(pfs_fit)$pval %>% round(., 5)



pdf(paste0(work_dir, "Fig_S6E_pfs_by_A&AN_versus_N_Survplot_in_IO.pdf"), width = 8, height = 6)
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







MFP = read.table("../resource/MFP_results/SCLC_TME_score.txt")
MFP = MFP[rownames(annotation),]

annotation$CAF = MFP$CAF
annotation$Matrix = MFP$Matrix
annotation$Matrix_remodeling = MFP$Matrix_remodeling


for(i_sig in c("CAF", "Matrix", "EMT_signature")){
        tmp_MFP = data.frame(MFP, "SCLC_subtype" = annotation$SCLC_subtype)
        tmp_MFP = tmp_MFP[,c(i_sig, "SCLC_subtype")]
        tmp_MFP = tmp_MFP[tmp_MFP$SCLC_subtype!="TN",]
        tmp_MFP$SCLC_subtype = factor(tmp_MFP$SCLC_subtype, levels = c("A", "AN", "N", "P"))
        colnames(tmp_MFP) = c("score", "SCLC_subtype")
        p = ggplot(tmp_MFP, aes(x=SCLC_subtype, y = score, fill = SCLC_subtype)) + geom_boxplot() + theme_bw()+ scale_fill_nejm() + xlab("SCLC subtype") +
                ylab(paste0(gsub("_signature", "", i_sig), " signature score")) + axis_theme + ylim(-2, 3) + theme(legend.position="none") + xlab("")
        ggsave(paste0(work_dir, "Fig_S6B_", i_sig,".png"), plot =p, width = 3.2, height = 3.8)
        print(wilcox.test(tmp_MFP$score[tmp_MFP$SCLC_subtype %in% c("A", "AN")], tmp_MFP$score[tmp_MFP$SCLC_subtype %in% c("N")])$p.value)
}




#################################### ASCL1 trans subtype
############################################################################
SCLC_subtype_WTS_meta$Inflamed_cluster = as.character(SCLC_subtype_WTS_meta$Immune_cluster)
SCLC_subtype_WTS_meta$Inflamed_cluster[SCLC_subtype_WTS_meta$Inflamed_cluster=="high_cluster"] = "Inflamed"
SCLC_subtype_WTS_meta$Inflamed_cluster[SCLC_subtype_WTS_meta$Inflamed_cluster=="low_cluster"] = "Non-Inflamed"
SCLC_subtype_WTS_meta$Inflamed_cluster = factor(SCLC_subtype_WTS_meta$Inflamed_cluster, levels = c("Inflamed", "Non-Inflamed"))

SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = NE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$NE_25_genelist = SingScore$TotalScore

SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = nonNE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$nonNE_25_genelist = SingScore$TotalScore



p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("A", "AN"), ], x = "Net_notch_score", y = "InflamedScore", fill = "Inflamed_cluster",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.2, label.y=0.1, label.sep = "\n"), palette = c("#E72644", "#86CFEC"),
   xlim = c(-0.3, 0.3), ylim = c(-0.4, 0.1))
ggsave(paste0(work_dir, "Fig_S6C_edit_ASCL1_AN_subtype_x_Net_NOTCH_y_InflamedScore_sig.png"), plot = p, width = 7, height = 7)



p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("A", "AN"), ], x = "NE_25_genelist", y = "InflamedScore", fill = "Inflamed_cluster",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.25, label.y=0.08, label.sep = "\n"), palette = c("#E72644", "#86CFEC"),
   xlim = c(0.0, 0.3), ylim = c(-0.4, 0.1))
ggsave(paste0(work_dir, "Fig_S6C_edit_ASCL1_AN_subtype_x_NE_y_InflamedScore.png"), plot = p, width = 7, height = 7)




p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("N"), ], x = "Net_notch_score", y = "InflamedScore", fill = "Inflamed_cluster",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.4, label.y=0.05, label.sep = "\n"), palette = c("#E72644", "#86CFEC"),
   ylim = c(-0.4, 0.1))
ggsave(paste0(work_dir, "Fig_S6D_edit_NEUROD1_subtype_x_Net_NOTCH_y_InflamedScore.png"), plot = p, width = 7, height = 7)



p = ggscatter(SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("N"), ], x = "NE_25_genelist", y = "InflamedScore", fill = "Inflamed_cluster",
   color = "black", shape = 21, size = 3, # Points color, shape and size
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE, # Add confidence interval
   cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
   cor.coeff.args = list(method = "pearson", label.x = 0.15, label.y=0.05, label.sep = "\n"), palette = c("#E72644", "#86CFEC"),
   ylim = c(-0.4, 0.1))
ggsave(paste0(work_dir, "Fig_S6D_edit_NEUROD1_subtype_x_NE_y_InflamedScore.png"), plot = p, width = 7, height = 7)


