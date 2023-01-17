LibPath = paste0(getwd(), "/LibPath/")
.libPaths(LibPath)
## load library and matrix file
library(NMF)
library(sva)
library(RColorBrewer)
library(GSVA)
library(tweeDEseq)
library(ComplexHeatmap)
library(pracma)
library(dplyr); library(stringr); library(tibble); library(magrittr)
library(ggplot2)
library(scales)
library(caret)
library(org.Hs.eg.db)
library(ggpubr)
library(fgsea)
library(circlize)
library(singscore)
library("xlsx")

ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Accuracy", percent_format()(m$overall[1]),
                   "Kappa", percent_format()(m$overall[2]))
  p <-
    ggplot(data = as.data.frame(m$table) ,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
    theme(legend.position = "none") +
    ggtitle(mytitle)
  return(p)
}
Signature_df = read.xlsx("../ref_signature/Signature gene list.xlsx", header = TRUE, sheetName="Sheet1")
NE_25_genelist = strsplit(Signature_df[Signature_df$Signature.name=="Neuroendocrine",2], ", ")[[1]]

Exp.count = read.table("../resource/Exp.count.txt", header = TRUE, sep = '\t')
SCLC_TPM =read.table("../resource/log2_TPM_n226.txt", sep = "\t")


## Preprocesing
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



## get CPM
Exp.count.melt <-reshape2::melt(Exp.count)
library(ggplot2)
library(edgeR)
library(EnhancedVolcano)
work_dir = "../Figures/step1/"
system(paste0("mkdir ", work_dir))
SCLC_meta = read.xlsx("../resource/Essential_check_Edited.xlsx", header = TRUE, sheetName = "Sheet1")

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

############################################
# write cpm matrix
############################################

gp <- rep(1,ncol(Exp.count))
Exp.count2 <- Exp.count
rownames(Exp.count2) <- toTable(org.Hs.egSYMBOL)$gene_id[match(rownames(Exp.count), toTable(org.Hs.egSYMBOL)$symbol)]
y <- DGEList(counts=Exp.count,
             gene=rownames(Exp.count),
             group=gp)
y$genes$ids <- rownames(Exp.count2)

keep = filterByExpr(y, gp,
                    min.count = 3, min.total.count = 10,
                    large.n = 2, min.prop = 0.2)
table(keep) #TRUE = 14554

y <- y[keep, ,keep.lib.size=FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ 0 + gp, data = y$samples)
y <- estimateDisp(y)
log2.CPM.count <- cpm(y, prior.count=2, log=TRUE)

# Not filtered CPM
y <- DGEList(counts=Exp.count,
             gene=rownames(Exp.count),
             group=gp)
y$genes$ids <- rownames(Exp.count2)
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ 0 + gp, data = y$samples)
y <- estimateDisp(y)
log2_CPM_n226_nonGeneFilter <- cpm(y, prior.count=2, log=TRUE)

SCLC_cpm = log2.CPM.count[,SCLC_subtype_WTS_meta$WTS_ID]
SCLC_TPM2 = SCLC_TPM[,SCLC_subtype_WTS_meta$WTS_ID]
############################################
# A versus N
############################################



sample_group_info = as.character(SCLC_subtype_WTS_meta$SCLC_subtype)

Exp.count2 = Exp.count[,SCLC_subtype_WTS_meta$WTS_ID]
y <- DGEList(counts=Exp.count2,
             gene=rownames(Exp.count2),
             group=sample_group_info)
dim(y)
keep <- filterByExpr(y, group=sample_group_info)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")
y$samples
eff.lib.size <- y$samples$lib.size*y$samples$norm.factors
Cpm.count <- edgeR::cpm(y, log=FALSE)
#boxplot(Cpm.count, col="gray", las=3)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
N_versus_A <- exactTest(y, pair=c("A","N")) # compare groups 1 and 2
topTags(N_versus_A, n=10)

plot_N_versus_A = EnhancedVolcano(N_versus_A$table,
                lab = rownames(N_versus_A$table),
                 x = 'logFC',
                 y = 'PValue',
                 pCutoff = 0.05,
                 pointSize = 1.5,
                 labSize = 5.0)
write.table(N_versus_A$table, paste0(work_dir, "plot_N_versus_A_DEG.txt"), col.names = TRUE, row.names = TRUE, sep = '\t', quote = FALSE)
saveRDS(N_versus_A, paste0(work_dir, "plot_N_versus_A_DEG.Rds"))


N_versus_A = N_versus_A$table
Filt_DEG = N_versus_A[N_versus_A$PValue<0.001, ]
Filt_DEG = Filt_DEG[rownames(Filt_DEG)%in%rownames(log2.CPM.count),]



SCLC_subtype_WTS_meta_NE = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$NE_subtype=="NE", ]

NE_SCLC_TPM2 = SCLC_TPM2[rownames(Filt_DEG),SCLC_subtype_WTS_meta_NE$WTS_ID]


calculate_group_centroid = function(Axis1, Axis2, group) {
        unique_group = unique(group)
        group_mean_Axis1 = lapply(unique_group, function(x) {
                                          return(mean(Axis1[group==x]))
        }) %>% as.numeric()
        group_mean_Axis2 = lapply(unique_group, function(x) {
                                          return(mean(Axis2[group==x]))
        }) %>% as.numeric()
        centroid_df = data.frame(group_mean_Axis1, group_mean_Axis2)
        rownames(centroid_df) = unique_group
        return(dist(centroid_df))
}



pca_df <- prcomp(t(NE_SCLC_TPM2),
                 center = T,
                 scale. = T)


Fig1D_PCA_df = data.frame(pca_df$x[,"PC1"], pca_df$x[, "PC2"], pca_df$x[,"PC3"], SCLC_subtype_WTS_meta_NE$SCLC_subtype)
colnames(Fig1D_PCA_df) = c("PC1", "PC2", "PC3", "SCLC_subtype")


p = ggscatter(Fig1D_PCA_df[Fig1D_PCA_df$PC1<40,], x="PC1", y = "PC2", shape = "SCLC_subtype", color = "SCLC_subtype", palette = c("#BC3C29", "#0072B5", "#E18727"), ellipse=TRUE,
         conf.int = FALSE, cor.coef = FALSE, ellipse.type = "convex")
ggsave(paste0(work_dir, "Fig1D_PCA_plot.pdf"), plot = p, width = 5, height = 5)

calculate_group_centroid(Fig1D_PCA_df[, "PC1"], Fig1D_PCA_df[,"PC2"], SCLC_subtype_WTS_meta_NE$SCLC_subtype)


##### Fig 1E right


DEG_dataframe = data.frame()

for(i_subtype in c("A", "AN", "N")){
        sample_group_info = as.character(SCLC_subtype_WTS_meta$SCLC_subtype)
        set.seed(100)
        sample_group_info[!(seq(1, length(sample_group_info)) %in% c(sample(seq(1,length(sample_group_info))[sample_group_info=="A"], 15),
                          sample(seq(1,length(sample_group_info))[sample_group_info=="AN"], 15),
                          sample(seq(1,length(sample_group_info))[sample_group_info=="N"], 15)))] = "except"

        sample_group_info[!(sample_group_info%in%c("A", "AN", "N"))] = "except"
        sample_group_info[sample_group_info!=i_subtype & sample_group_info!="except"] = "Others"
        Exp.count2 = Exp.count[,SCLC_subtype_WTS_meta$WTS_ID]
        y <- DGEList(counts=Exp.count2,
                     gene=rownames(Exp.count2),
                     group=sample_group_info)
        dim(y)
        keep <- filterByExpr(y, group=sample_group_info)
        y <- y[keep, , keep.lib.sizes=FALSE]
        y <- calcNormFactors(y, method="TMM")
        eff.lib.size <- y$samples$lib.size*y$samples$norm.factors
        Cpm.count <- edgeR::cpm(y, log=FALSE)
        y <- estimateCommonDisp(y)
        y <- estimateTagwiseDisp(y)
        DEG_object <- exactTest(y, pair=c("Others", i_subtype)) # compare groups 1 and 2
        DEG_df = DEG_object$table
        DEG_df$Gene = rownames(DEG_df)
        DEG_df$check = DEG_df$logFC#rep(1, nrow(DEG_df))
        colnames(DEG_df) = c(colnames(DEG_df)[seq(1,4)], i_subtype)
        assign(paste0("DEG_", i_subtype), DEG_df)
        DEG_df = DEG_df[DEG_df$logFC > 1.5,]
        DEG_df = DEG_df[order(DEG_df$PValue),]
	print(summary(DEG_df$PValue[seq(1,50)]))
        DEG_df = DEG_df[seq(1,50),c("Gene", i_subtype)]
        assign(paste0("DEG_", i_subtype), DEG_object$table)
        if(nrow(DEG_dataframe) > 0){
                DEG_dataframe = merge(DEG_dataframe, DEG_df, all = TRUE, by = "Gene")
        }else{
                DEG_dataframe = DEG_df
        }
}


######### Get all


identical(colnames(SCLC_TPM2), SCLC_subtype_WTS_meta$WTS_ID)
sample_group_info = as.character(SCLC_subtype_WTS_meta$SCLC_subtype)
print(table(sample_group_info))

######### return value

AverageExpr_DEG = lapply(DEG_dataframe$Gene[DEG_dataframe$Gene %in% rownames(SCLC_TPM2)], function(x) {
                x = as.character(x[1])
                expr_vector = SCLC_TPM2[x,] %>% unlist() %>% as.numeric()
               mean_A_type = mean(expr_vector[sample_group_info=="A"])
               mean_trans_type = mean(expr_vector[sample_group_info=="AN"])
               mean_N_type = mean(expr_vector[sample_group_info=="N"])
               return(c(mean_A_type, mean_trans_type, mean_N_type))
}) %>% data.frame() %>% t()
colnames(AverageExpr_DEG) = c("A", "AN", "N")

rownames(AverageExpr_DEG) = DEG_dataframe$Gene[DEG_dataframe$Gene %in% rownames(SCLC_TPM2)]
dend = as.dendrogram(hclust(dist(cor(AverageExpr_DEG))))
dd.reorder <- reorder(dend, c(3,1,2))

col_fun = colorRamp2(c(0.8, 0.9,1), c("#3a1b59", "#2E4A89", "#FEFA2B"))
print(cor(AverageExpr_DEG))
pdf(paste0(work_dir, "Fig1D_right_NeuroEndocrine_AverageExpr_DEG_correlation_heatmap.pdf"), width = 8, height = 7)
Heatmap(cor(AverageExpr_DEG), cluster_columns=rev(dd.reorder), cluster_rows=rev(dd.reorder), col = col_fun, name = "FC>1.5 & P value top50")
dev.off()




DEG_dataframe = data.frame()

for(i_subtype in c("A", "AN","N", "P")){
        sample_group_info = as.character(SCLC_subtype_WTS_meta$SCLC_subtype)
        set.seed(100)
        sample_group_info[!(seq(1, length(sample_group_info)) %in% c(sample(seq(1,length(sample_group_info))[sample_group_info=="A"], 15),
                          sample(seq(1,length(sample_group_info))[sample_group_info=="AN"], 15),
                          sample(seq(1,length(sample_group_info))[sample_group_info=="N"], 15),
                          sample(seq(1,length(sample_group_info))[sample_group_info=="P"], 15)))] = "except"

        sample_group_info[!(sample_group_info%in%c("A", "AN", "N", "P"))] = "except"
        sample_group_info[sample_group_info!=i_subtype & sample_group_info!="except"] = "Others"
        Exp.count2 = Exp.count[,SCLC_subtype_WTS_meta$WTS_ID]
        y <- DGEList(counts=Exp.count2,
                     gene=rownames(Exp.count2),
                     group=sample_group_info)
        dim(y)
        keep <- filterByExpr(y, group=sample_group_info)
        y <- y[keep, , keep.lib.sizes=FALSE]
        y <- calcNormFactors(y, method="TMM")
        eff.lib.size <- y$samples$lib.size*y$samples$norm.factors
        Cpm.count <- edgeR::cpm(y, log=FALSE)
        y <- estimateCommonDisp(y)
        y <- estimateTagwiseDisp(y)
        DEG_object <- exactTest(y, pair=c("Others", i_subtype)) # compare groups 1 and 2
        DEG_df = DEG_object$table
        DEG_df$Gene = rownames(DEG_df)
        DEG_df$check = DEG_df$logFC#rep(1, nrow(DEG_df))
        colnames(DEG_df) = c(colnames(DEG_df)[seq(1,4)], i_subtype)
        assign(paste0("DEG_", i_subtype), DEG_df)
        DEG_df = DEG_df[DEG_df$logFC > 1.25,]
        DEG_df = DEG_df[order(DEG_df$PValue),]
	print(summary(DEG_df$PValue[seq(1,50)]))
        DEG_df = DEG_df[seq(1,50),c("Gene", i_subtype)]
        assign(paste0("DEG_", i_subtype), DEG_object$table)
        if(nrow(DEG_dataframe) > 0){
                DEG_dataframe = merge(DEG_dataframe, DEG_df, all = TRUE, by = "Gene")
        }else{
                DEG_dataframe = DEG_df
        }
}




identical(colnames(SCLC_TPM2), SCLC_subtype_WTS_meta$WTS_ID)
sample_group_info = as.character(SCLC_subtype_WTS_meta$SCLC_subtype)
print(table(sample_group_info))

######### return value

AverageExpr_DEG = lapply(DEG_dataframe$Gene[DEG_dataframe$Gene %in% rownames(SCLC_TPM2)], function(x) {
                x = as.character(x[1])
                expr_vector = SCLC_TPM2[x,] %>% unlist() %>% as.numeric()
               mean_A_type = mean(expr_vector[sample_group_info=="A"])
               mean_trans_type = mean(expr_vector[sample_group_info=="AN"])
               mean_N_type = mean(expr_vector[sample_group_info=="N"])
               mean_P_type = mean(expr_vector[sample_group_info=="P"])
               mean_TN_type = mean(expr_vector[sample_group_info=="TN"])
               return(c(mean_A_type, mean_trans_type, mean_N_type, mean_P_type, mean_TN_type))
}) %>% data.frame() %>% t()
colnames(AverageExpr_DEG) = c("A", "AN", "N", "P", "TN")

rownames(AverageExpr_DEG) = DEG_dataframe$Gene[DEG_dataframe$Gene %in% rownames(SCLC_TPM2)]
dend = as.dendrogram(hclust(dist(cor(AverageExpr_DEG))))
dd.reorder <- reorder(dend, c(2, 1 , 3, 4))

col_fun = colorRamp2(c(0.6, 0.8,1), c("#3a1b59", "#2E4A89", "#FEFA2B"))
pdf(paste0(work_dir, "Fig1D_left_NeuroEndo_and_P_AverageExpr_DEG_correlation_heatmap.pdf"), width = 8, height = 8)
set.seed(0)
Heatmap(cor(AverageExpr_DEG), cluster_columns=TRUE, cluster_rows=TRUE, col = col_fun, name = "FC>1.25 & P value top50")
dev.off()

########################################
## N versus A&trans
############################

sample_group_info = lapply(as.character(SCLC_subtype_WTS_meta$SCLC_subtype), function(x) {if(x=="A" | x=="AN"){
               return("A&AN")
                 }else{return(x)}}) %>% as.character()

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

N_versus_AandAN <- exactTest(y, pair=c("A&AN","N")) # compare groups 1 and 2
topTags(N_versus_AandAN, n=10)

target_gene_of_N_vs_AAN = c("ASCL1", "CD74", "CXCL10", "CXCL9", "CCL18", "CXCL13", "IL7", "IL23A", "GZMH",
                                            "CLKNK", "CLECL1", "MCTP2", "COTL1", "ETS2", "SEC11C", "HEYL", "PTCRA", "DVL2",
                                              "NOTCH2", "COL9A1", "NEUROD1", "MMP8", "NOTCH1", "FGFR3","SERPINA10", "TIMP3")
N_AAN_DEG = N_versus_AandAN$table
N_AAN_DEG$gene = rownames(N_AAN_DEG)

keyvals = apply(N_AAN_DEG, 1, function(x) {
                        if(as.numeric(x[3])< 0.05 & abs(as.numeric(x[1])) >= 1){
                                if(as.numeric(x[1]) >= 0.1){
                                        if(sum(x[4]%in%target_gene_of_N_vs_AAN)>0){
                                                return("#1e2f97")
                                        }else{return("#949494")}
                                }else if(as.numeric(x[1]) <= -0.1){
                                        if(sum(x[4]%in%target_gene_of_N_vs_AAN)>0){
                                                return("#800000")
                                        }else{return("#949494")}
                                }
                        }else{return("#949494")}
                                              }) %>% unlist() %>% as.character()

names(keyvals)[keyvals == '#800000'] <- 'A&AN Immune'
#names(keyvals)[keyvals == '#c38e63'] <- 'nonhigh'# "#f6bdc0"
names(keyvals)[keyvals == "#1e2f97"] <- 'N Notch'
#names(keyvals)[keyvals == '#59788e'] <- 'nonlow' #"#b0dbf1"
names(keyvals)[keyvals == '#949494'] <- 'non-sig'
N_AAN_DEG$log2FC = log2(exp(N_AAN_DEG$logFC))

plot_N_AAN_DEG = EnhancedVolcano(N_AAN_DEG,
                lab = rownames(N_AAN_DEG),
                 x = 'log2FC', labSize=3,
                 y = 'PValue', xlim = c(-5.6, 5.6), ylim = c(0, 15),
                 pCutoff = 0.05, colCustom=keyvals,
                 pointSize = 1, selectLab=target_gene_of_N_vs_AAN,
                 drawConnectors = TRUE, directionConnectors = "both",  max.overlaps=Inf, colAlpha = 0.9)
ggsave(paste0(work_dir, "Fig_S5B_edgeR_N_AAN_DEG.png"), plot = plot_N_AAN_DEG, width = 5, height = 6)

plot_N_AAN_DEG = EnhancedVolcano(N_AAN_DEG,
                lab = rownames(N_AAN_DEG),
                 x = 'log2FC', labSize=3,
                 y = 'PValue', xlim = c(-5.6, 5.6), ylim = c(0, 15),
                 pCutoff = 0.05, colCustom=keyvals,
                 pointSize = 1, selectLab=target_gene_of_N_vs_AAN,
                 drawConnectors = FALSE, directionConnectors = "both",  max.overlaps=Inf, colAlpha = 0.9)
ggsave(paste0(work_dir, "Fig_S5B_edgeR_N_AAN_DEG_nodir.png"), plot = plot_N_AAN_DEG, width = 5, height = 6)


saveRDS(N_versus_AandAN, paste0(work_dir, "Fig_S5B_edgeR_N_AAN_DEG.Rds"))
########################################
pathways.hallmark <- gmtPathways("../resource/h.all.v7.4.symbols.gmt")
#system("wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.4/c2.cp.kegg.v7.4.symbols.gmt -P ../resource/")
pathway.kegg = gmtPathways("../resource/c2.cp.kegg.v7.4.symbols.gmt")

Hallmark_immune = pathways.hallmark[grepl("HALLMARK_I", names(pathways.hallmark))] %>% unlist() %>% as.character()
Hallmark_immune[Hallmark_immune %in% rownames(N_AAN_DEG)[N_AAN_DEG$logFC < (-0.5) & N_AAN_DEG$PValue < 0.001]]

NOTCH = c(pathway.kegg$KEGG_NOTCH_SIGNALING_PATHWAY, pathways.hallmark$HALLMARK_NOTCH_SIGNALING) %>% unique()
N_AAN_DEG[N_AAN_DEG$logFC>0 & rownames(N_AAN_DEG) %in% NOTCH,]

########################################
## P versus A&AN
############################

sample_group_info = lapply(as.character(SCLC_subtype_WTS_meta$SCLC_subtype), function(x) {if(x=="A" | x=="AN"){
               return("A&AN")
                 }else{return(x)}}) %>% as.character()

Exp.count2 = Exp.count[,SCLC_subtype_WTS_meta$WTS_ID]
y <- DGEList(counts=Exp.count2,
             gene=rownames(Exp.count2),
             group=sample_group_info)
dim(y)
keep <- filterByExpr(y, group=sample_group_info)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")
y$samples
Cpm.count <- edgeR::cpm(y, log=FALSE)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

P_versus_AandAN <- exactTest(y, pair=c("A&AN","P")) # compare groups 1 and 2
topTags(P_versus_AandAN, n=10)


target_gene_of_P_vs_AAN = c("ASCL1", "POU2F3", "NOTCH1", "NOTCH2",
                                            "HLA-DQA1", "IL7R", "IRF8", "CD4", "LYN", "IL17RB", "CD8A", "IL1R2", "PTGER4", "CCR3", "CXCL3", "CXCL2",
                                            "SAMHD1", "PIK3R5", "RNF166", "SPIC", "SPIB")


P_AAN_DEG = P_versus_AandAN$table
P_AAN_DEG$gene = rownames(P_AAN_DEG)
keyvals = apply(P_AAN_DEG, 1, function(x) {
                        if(as.numeric(x[3])< 0.05 & abs(as.numeric(x[1])) >= 1){
                                if(as.numeric(x[1]) >= 0.1){
                                        if(sum(x[4]%in%target_gene_of_P_vs_AAN)>0){
                                                return('#800000')
                                        }else{return("#949494")}
                                }else if(as.numeric(x[1]) <= -0.1){
                                        if(sum(x[4]%in%c(target_gene_of_P_vs_AAN, NE_25_genelist))>0){
                                                return('black')
                                        }else{return("#949494")}
                                }
                        }else{return("#949494")}
                                              }) %>% unlist() %>% as.character()

names(keyvals)[keyvals == 'black'] <- 'A&AN NE'
names(keyvals)[keyvals == '#800000'] <- 'P immune'
names(keyvals)[keyvals == '#949494'] <- 'non-sig'

VolcanoPlot_genelist = c(target_gene_of_P_vs_AAN, NE_25_genelist)
VolcanoPlot_genelist = VolcanoPlot_genelist[VolcanoPlot_genelist%in%rownames(P_AAN_DEG)]
VolcanoPlot_genelist = VolcanoPlot_genelist[P_AAN_DEG[VolcanoPlot_genelist,"PValue"] <0.008]
P_AAN_DEG$log2FC = log2(exp(P_AAN_DEG$logFC))

plot_P_AAN_DEG = EnhancedVolcano(P_AAN_DEG,
                lab = rownames(P_AAN_DEG),
                 x = 'log2FC', labSize=3,
                 y = 'PValue', ylim = c(0, 30), xlim = c(-10, 10),
                 pCutoff = 0.05,
                 pointSize = 1, selectLab=VolcanoPlot_genelist, colCustom=keyvals,
                 drawConnectors = TRUE,directionConnectors = "both",  max.overlaps=Inf, colAlpha = 0.9)
ggsave(paste0(work_dir, "Fig_S5A_edgeR_P_AAN_DEG.png"), plot = plot_P_AAN_DEG, width = 5, height = 6)
saveRDS(P_versus_AandAN, paste0(work_dir, "Fig_S5A_edgeR_P_AAN_DEG.Rds"))

###############################

plot_P_AAN_DEG = EnhancedVolcano(P_AAN_DEG,
                lab = rownames(P_AAN_DEG),
                 x = 'log2FC', labSize=3,
                 y = 'PValue', ylim = c(0, 30), xlim = c(-10, 10),
                 pCutoff = 0.05,
                 pointSize = 1, selectLab=VolcanoPlot_genelist, colCustom=keyvals,
                 drawConnectors = FALSE,directionConnectors = "both",  max.overlaps=Inf, colAlpha = 0.9)
ggsave(paste0(work_dir, "Fig_S5A_edgeR_P_AAN_DEG_nodir.png"), plot = plot_P_AAN_DEG, width = 5, height = 6)


############################### WTS subype
thickness = theme(legend.title = element_text(size = 13, color = "black", face = "bold"),
      legend.text = element_text(size = 13, color = "black", face = "bold"),
      axis.text.x = element_text(size = 13, color = "black", face = "bold"),
      axis.text.y = element_text(size = 13, color = "black", face = "bold"),
      axis.title.x = element_text(size = 13, color = "black", face = "bold"),
      axis.title.y = element_text(size = 13, color = "black", face = "bold"))



SCLC_subtype_WTS_meta = SCLC_meta[SCLC_meta$WTS_QC_Result %in% c("Pass", "1") & grepl("includ", SCLC_meta$study.inclusion),]

rownames(SCLC_subtype_WTS_meta) = SCLC_subtype_WTS_meta$WTS_ID
SCLC_subtype_WTS_meta$SCLC_subtype[is.na(SCLC_subtype_WTS_meta$SCLC_subtype)] = "non-IHC"
SCLC_subtype_WTS_meta$SCLC_subtype = factor(SCLC_subtype_WTS_meta$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN", "non-IHC"))
SCLC_TPM2 = SCLC_TPM[,SCLC_subtype_WTS_meta$WTS_ID]



SCLC_TPM2_Rank <- rankGenes(SCLC_TPM2)

SingScore <- simpleScore(
  rankData = SCLC_TPM2_Rank,
  upSet = NE_25_genelist,
  centerScore = T,
  knownDirection = T
  )
SCLC_subtype_WTS_meta$NE_25_genelist = SingScore$TotalScore

###########################

for(i_name in c("ASCL1 & ND1 Shared targets", "ASCL1 high signatures cell line",
                "ND1 high signatures cell line")){
        gene_list = strsplit(Signature_df[Signature_df$Signature.name==i_name, 2], ", ")[[1]] %>% unlist() %>% as.character() %>% unique()
            SingScore <- simpleScore(
          rankData = SCLC_TPM2_Rank,
          upSet = gene_list,
          centerScore = T,
          knownDirection = T
          )
        i_name = strsplit(i_name, ";")[[1]][1]
	if(i_name == "ASCL1 & ND1 Shared targets"){i_name="ASCL1_ND1_Shared_targets"}else if(i_name=="ASCL1 high signatures cell line"){
		i_name = "ASCL1_high_signatures_cell_line"}else if(i_name == "ND1 high signatures cell line"){i_name="ND1_high_signatures_cell_line"}
        SCLC_subtype_WTS_meta[,i_name] = SingScore$TotalScore
}

####################### TUFT cell marker, from all gene cpm


SingScore <- simpleScore(
  rankData = SCLC_TPM2_Rank,
  upSet = strsplit(Signature_df[Signature_df$Signature.name=="Tuft cell marker", 2], ", ")[[1]],
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$Tuft_cell_marker = SingScore$TotalScore

##################### zscaling
SCLC_subtype_WTS_meta$NEUROD1_tpm = as.numeric(SCLC_TPM2["NEUROD1",])
SCLC_subtype_WTS_meta$NEUROD1_tpm_z = (SCLC_subtype_WTS_meta$NEUROD1_tpm-mean(SCLC_subtype_WTS_meta$NEUROD1_tpm))/sd(SCLC_subtype_WTS_meta$NEUROD1_tpm)
SCLC_subtype_WTS_meta$ASCL1_tpm = as.numeric(SCLC_TPM2["ASCL1",])
SCLC_subtype_WTS_meta$ASCL1_tpm_z = (SCLC_subtype_WTS_meta$ASCL1_tpm-mean(SCLC_subtype_WTS_meta$ASCL1_tpm))/sd(SCLC_subtype_WTS_meta$ASCL1_tpm)
SCLC_subtype_WTS_meta$POU2F3_tpm = as.numeric(SCLC_TPM2["POU2F3",])
SCLC_subtype_WTS_meta$POU2F3_tpm_z = (SCLC_subtype_WTS_meta$POU2F3_tpm-mean(SCLC_subtype_WTS_meta$POU2F3_tpm))/sd(SCLC_subtype_WTS_meta$POU2F3_tpm)


SCLC_subtype_WTS_meta$Tuft_cell_marker_z = (SCLC_subtype_WTS_meta$Tuft_cell_marker-mean(SCLC_subtype_WTS_meta$Tuft_cell_marker))/sd(SCLC_subtype_WTS_meta$Tuft_cell_marker)
SCLC_subtype_WTS_meta$ASCL1_high_signatures_z = (SCLC_subtype_WTS_meta$ASCL1_high_signatures_cell_line-mean(SCLC_subtype_WTS_meta$ASCL1_high_signatures_cell_line))/sd(SCLC_subtype_WTS_meta$ASCL1_high_signatures_cell_line)

SCLC_subtype_WTS_meta$ND1_high_signatures_z = (SCLC_subtype_WTS_meta$ND1_high_signatures_cell_line-mean(SCLC_subtype_WTS_meta$ND1_high_signatures_cell_line))/sd(SCLC_subtype_WTS_meta$ND1_high_signatures_cell_line)

SCLC_subtype_WTS_meta$NE_25_genelist_z = (SCLC_subtype_WTS_meta$NE_25_genelist-mean(SCLC_subtype_WTS_meta$NE_25_genelist))/sd(SCLC_subtype_WTS_meta$NE_25_genelist)

##################### call WTS subtype

P_subtype = SCLC_subtype_WTS_meta$WTS_ID[SCLC_subtype_WTS_meta$Tuft_cell_marker_z > 1]# & SCLC_subtype_WTS_meta$POU2F3_tpm_z > 0.5]
TN_subtype = SCLC_subtype_WTS_meta$WTS_ID[SCLC_subtype_WTS_meta$Tuft_cell_marker_z < -0.25 &
                                          SCLC_subtype_WTS_meta$ASCL1_high_signatures_z < -0.25 &
                                          SCLC_subtype_WTS_meta$ND1_high_signatures_z < -0.25 &
                                          !(SCLC_subtype_WTS_meta$WTS_ID%in%P_subtype)]
N_subtype = SCLC_subtype_WTS_meta$WTS_ID[SCLC_subtype_WTS_meta$ND1_high_signatures_z > 0.25 & SCLC_subtype_WTS_meta$ASCL1_high_signatures_z < -0.25 &
                                         SCLC_subtype_WTS_meta$NEUROD1_tpm_z > 0.25 & SCLC_subtype_WTS_meta$ASCL1_tpm_z < -0.25]
N_subtype = N_subtype[!(N_subtype %in% c(P_subtype, TN_subtype))]
A_subtype = SCLC_subtype_WTS_meta$WTS_ID[SCLC_subtype_WTS_meta$ND1_high_signatures_z < 0]#&#quantile(SCLC_subtype_WTS_meta$ND1_high_signatures_z, 0.75) &
A_subtype = A_subtype[!(A_subtype %in% c(P_subtype, TN_subtype, N_subtype))]


################################ Fig 2B
require(survminer)
require(survival)



WTS_jjg_subtype = rep("AN", nrow(SCLC_subtype_WTS_meta))
WTS_jjg_subtype[SCLC_subtype_WTS_meta$WTS_ID %in% P_subtype] = "P"
WTS_jjg_subtype[SCLC_subtype_WTS_meta$WTS_ID %in% TN_subtype] = "TN"
WTS_jjg_subtype[SCLC_subtype_WTS_meta$WTS_ID %in% N_subtype] = "N"
WTS_jjg_subtype[SCLC_subtype_WTS_meta$WTS_ID %in% A_subtype] = "A"

identical(WTS_jjg_subtype, SCLC_subtype_WTS_meta$WTS_subtype)
SCLC_subtype_WTS_meta$WTS_subtype = factor(SCLC_subtype_WTS_meta$WTS_subtype, levels = c("A", "AN", "N", "P", "TN"))
SCLC_subtype_WTS_meta$WTS_subtype = factor(WTS_jjg_subtype, levels = c("A", "AN", "N", "P", "TN"))


ggdat <-  SCLC_subtype_WTS_meta %>% group_by (WTS_subtype) %>% summarise(count = n()) %>% mutate(prop = round(count/sum(count)*100, digits=1))
rownames(ggdat) = ggdat$WTS_subtype; ggdat=ggdat[rev(c("A", "AN", "N", "P", "TN")),]
ggdat <- ggdat %>% mutate(lab.ypos = cumsum(prop) - 0.5*prop)
ggdat$WTS_subtype = factor(ggdat$WTS_subtype, levels = c("A", "AN", "N", "P", "TN"))

p = ggplot(ggdat, aes(x = "", y = prop, fill = WTS_subtype)) +
  geom_bar(width = 2, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white", size = 3)+
  scale_fill_manual(values = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1")) +
  theme_void()

ggsave(paste0(work_dir,"Fig1C_right.WTS_subtype.png"), plot = p, width = 5, height = 5)



library(caret)
library(org.Hs.eg.db)
library(scales)
ggplotConfusionMatrix <- function(m){
  mytitle <- paste("Accuracy = ", percent_format()(m$overall[1]),
                   ", Kappa agreement = ", percent_format()(m$overall[2]))
  p <-
    ggplot(data = as.data.frame(m$table) ,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_text(aes(x = Reference, y = Prediction, label = Freq, fontface=2), size = 6) +
    theme(legend.position = "none") +
    ggtitle(mytitle) +
    theme(plot.title = element_text(size = 20, color = "black", face = "bold"),
                axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                axis.text.y = element_text(size = 20, color = "black", face = "bold"))
  return(p)
}


##### Calculate kappa
ggconfusion_P = SCLC_subtype_WTS_meta[, c("WTS_subtype", "SCLC_subtype")]
ggconfusion_P = ggconfusion_P[ggconfusion_P$SCLC_subtype!="non-IHC",]
ggconfusion_P = ggconfusion_P[ggconfusion_P$WTS_subtype!="TN" & ggconfusion_P$SCLC_subtype!="TN",]
kappa_df = ggconfusion_P
kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
kappa_df$WTS_subtype[kappa_df$WTS_subtype!="P"] = "non-P"
kappa_df$SCLC_subtype = as.character(kappa_df$SCLC_subtype)
kappa_df$SCLC_subtype[kappa_df$SCLC_subtype!="P"] = "non-P"
cfm3 <- confusionMatrix(factor(kappa_df$SCLC_subtype, levels = c("P", "non-P")),
			factor(kappa_df$WTS_subtype, levels = c("P", "non-P")))
print("Fig S1F P and nonP subtype confusion matrix\n")
print(cfm3)


ggconfusion_P = SCLC_subtype_WTS_meta[, c("WTS_subtype", "SCLC_subtype")]
ggconfusion_P = ggconfusion_P[ggconfusion_P$SCLC_subtype!="non-IHC",]
ggconfusion_P = ggconfusion_P[!ggconfusion_P$WTS_subtype %in% c("TN", "P") & !ggconfusion_P$SCLC_subtype %in% c("TN", "P"),]
kappa_df = ggconfusion_P
cfm3 <- confusionMatrix(factor(kappa_df$SCLC_subtype, levels = c("A", "AN", "N")),
                        factor(kappa_df$WTS_subtype, levels = c("A", "AN", "N")))
print("Fig S1F A, AN, and N subtype confusion matrix\n")
print(cfm3)


kappa_df = ggconfusion_P
kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
kappa_df$WTS_subtype[kappa_df$WTS_subtype %in% c("A", "AN")] = "A&AN"
kappa_df$SCLC_subtype = as.character(kappa_df$SCLC_subtype)
kappa_df$SCLC_subtype[kappa_df$SCLC_subtype %in% c("A", "AN")] = "A&AN"
cfm3 <- confusionMatrix(factor(kappa_df$SCLC_subtype, levels = c("N", "A&AN")),
                        factor(kappa_df$WTS_subtype, levels = c("N", "A&AN")))

print("Fig S1F, N and A&AN subtype confusion matrix")
print(cfm3)

kappa_df = ggconfusion_P
kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
kappa_df$WTS_subtype[kappa_df$WTS_subtype %in% c("N", "AN")] = "N&AN"
kappa_df$SCLC_subtype = as.character(kappa_df$SCLC_subtype)
kappa_df$SCLC_subtype[kappa_df$SCLC_subtype %in% c("N", "AN")] = "N&AN"
cfm3 <- confusionMatrix(factor(kappa_df$SCLC_subtype, levels = c("N&AN", "A")),
                        factor(kappa_df$WTS_subtype, levels = c("N&AN", "A")))

print("Fig S1F, N&AN and A subtype confusion matrix")
print(cfm3)




#####################
##### Fig S3


OS_fit <- survfit(Surv(OS_time, OS_event==1)~WTS_subtype, data = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$ini_stage=="ED",])
cox_fit = coxph(Surv(OS_time, OS_event==1)~WTS_subtype, data = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$ini_stage=="ED",])

pdf(paste0(work_dir, "Fig_S3D_OS_Analysis_WTS_subtype_ED.pdf"), width = 8, height = 6.2)
p = ggsurvplot(
  fit = OS_fit,
  size = 1.5,
  xlab = "Months since Treatment",
  ylab = "OS",
  surv.median.line = "hv",
  palette = "nejm",
  xlim=c(0,30),
  ylim=c(0,1),
  break.time.by=6,
  title="OS Analysis in ED-patients: WTS subtype",
  font.title=c("bold"), font.subtitle=c("italic"),
  legend.title="Strata",
  pval = T, pval.coord = c(1, 0.25),
  surv.scale = "percent",
  risk.table = TRUE, risk.table.height = 0.3, conf.int = F)
p$plot <- p$plot + thickness
print(p)
dev.off()



OS_fit_summary = data.frame(summary(OS_fit)$table)
rownames(OS_fit_summary) = gsub("=", "", rownames(OS_fit_summary))
OS_fit_summary$rownames = rownames(OS_fit_summary)


cox_summary = cbind(data.frame(summary(cox_fit)$coefficients), data.frame(summary(cox_fit)$conf.int))
cox_summary$rownames = rownames(cox_summary)


overall_summary = merge(OS_fit_summary, cox_summary, by="rownames", all=TRUE)

overall_summary = (data.frame(overall_summary$rownames,
                 paste0(round(overall_summary$median, 1), " (", round(overall_summary$X0.95LCL, 1), "-", round(overall_summary$X0.95UCL, 1), ")"),
                 paste0(round(overall_summary$exp.coef., 1), " (", round(overall_summary$lower..95, 1), "-", round(overall_summary$upper..95, 1),")"),
                 overall_summary$Pr...z..))
colnames(overall_summary) = c("group", "mOS", "HR", "P")
overall_summary$group = gsub(" ", "", overall_summary$group)
overall_summary$P = round(overall_summary$P, 3)

apply(overall_summary, 1, function(x) {paste0(x[seq(1,4)], collapse = ";")})



###############
##### Fig S3F

i_distance_method = "manhattan"
i_cluster_medhod = "average"
pt.matrix = t(SCLC_subtype_WTS_meta[, c("Tuft_cell_marker", "POU2F3_tpm", "NE_25_genelist")])

pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))

set.seed(100)

k_cluster = 2
cluster_heatmap = Heatmap(pt.matrix, column_split = k_cluster, clustering_distance_columns = i_distance_method, show_column_names = FALSE, clustering_method_columns = i_cluster_medhod, column_gap = unit(5, "mm"), cluster_rows=FALSE)
cluster_heatmap_ht = draw(cluster_heatmap); cluster_heatmap_tr_order = column_order(cluster_heatmap_ht)
pdf(paste0(work_dir, "Fig_S1D_WTS_clustering_P_classification_",k_cluster,"_",i_cluster_medhod,"_",i_distance_method,".pdf"), width = 6, height = 2)
print(cluster_heatmap_ht)
dev.off()

POU2F3_index = which.max(c(mean(pt.matrix["Tuft_cell_marker",cluster_heatmap_tr_order[[1]]]),
			    mean(pt.matrix["Tuft_cell_marker",cluster_heatmap_tr_order[[2]]])))
tmp = rep("nonP", ncol(pt.matrix))
tmp[seq(1,ncol(pt.matrix)) %in% cluster_heatmap_tr_order[[POU2F3_index]]] = "P"
print(table(tmp))
column_split = tmp

rownames(pt.matrix) = c(1,2,3)
pdf(paste0(work_dir, "Fig_S1D_WTS_clustering_P_classification_",k_cluster,"_",i_cluster_medhod,"_",i_distance_method,".pdf"), width = 6, height = 2)
Heatmap(pt.matrix, column_split = factor(column_split, levels = c("nonP", "P")), clustering_distance_columns = i_distance_method, show_column_names = FALSE, clustering_method_columns = i_cluster_medhod, column_gap = unit(5, "mm"), cluster_rows=FALSE)
dev.off()


SCLC_TPM2 = rbind(SCLC_TPM[c("NEUROD1", "ASCL1"),SCLC_subtype_WTS_meta$WTS_ID], t(SCLC_subtype_WTS_meta[, c("ASCL1_high_signatures_cell_line", "ND1_high_signatures_cell_line", "NE_25_genelist")]))
SCLC_TPM2 = SCLC_TPM2[c("ASCL1_high_signatures_cell_line", "ASCL1", "ND1_high_signatures_cell_line", "NEUROD1"), column_split!="P"]

pt.matrix <- t(apply(SCLC_TPM2,1,function(x){(x-mean(x))/sd(x)}))

i_distance_method = "maximum"
i_cluster_medhod = "complete"


set.seed(100)
#pt.matrix2 = pt.matrix
rownames(pt.matrix) = c(1,2,3,4)
cluster_heatmap = Heatmap(pt.matrix, column_split = 3, clustering_distance_columns = i_distance_method, show_column_names = FALSE, clustering_method_columns = i_cluster_medhod, column_gap = unit(5, "mm"), cluster_rows=FALSE)
cluster_heatmap_ht = draw(cluster_heatmap); cluster_heatmap_tr_order = column_order(cluster_heatmap_ht)

dev.off()
pdf(paste0(work_dir, "Fig_S1D_WTS_clustering_NE_classification_",k_cluster,"_",i_cluster_medhod,"_",i_distance_method,".pdf"), width = 6, height = 3)
print(cluster_heatmap_ht)
dev.off()
tmp = rep("cluster", ncol(pt.matrix))
tmp[cluster_heatmap_tr_order[[1]]] = "1"; tmp[cluster_heatmap_tr_order[[2]]] = "2"; tmp[cluster_heatmap_tr_order[[3]]] = "3";
surv_df_SCLC_subtype_WTS_meta = SCLC_subtype_WTS_meta
surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster = rep("nonAnno", nrow(surv_df_SCLC_subtype_WTS_meta))
surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster[!surv_df_SCLC_subtype_WTS_meta$WTS_ID %in% colnames(pt.matrix)] = "P"

rownames(pt.matrix) = c("ASCL1_high_signatures_cell_line", "ASCL1", "ND1_high_signatures_cell_line", "NEUROD1")

NEUROD1_index = which.max(c(mean(pt.matrix["ND1_high_signatures_cell_line",cluster_heatmap_tr_order[[1]]]),
                            mean(pt.matrix["ND1_high_signatures_cell_line",cluster_heatmap_tr_order[[2]]]),
                            mean(pt.matrix["ND1_high_signatures_cell_line",cluster_heatmap_tr_order[[3]]])))
ASCL1_index = which.max(c(mean(pt.matrix["ASCL1_high_signatures_cell_line",cluster_heatmap_tr_order[[1]]]),
                          mean(pt.matrix["ASCL1_high_signatures_cell_line",cluster_heatmap_tr_order[[2]]]),
                          mean(pt.matrix["ASCL1_high_signatures_cell_line",cluster_heatmap_tr_order[[3]]])))


surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster[surv_df_SCLC_subtype_WTS_meta$WTS_ID %in% colnames(pt.matrix)[seq(1, ncol(pt.matrix)) %in% cluster_heatmap_tr_order[[NEUROD1_index]]]] = "N"
surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster[surv_df_SCLC_subtype_WTS_meta$WTS_ID %in% colnames(pt.matrix)[seq(1, ncol(pt.matrix)) %in% cluster_heatmap_tr_order[[ASCL1_index]]]] = "A"

surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster[!surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster %in% c("A", "N", "P")] = "AN"
surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster = factor(surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster, levels = c("A", "AN", "N", "P"))


OS_fit <- survfit(Surv(OS_time, OS_event==1)~Heatmap_cluster, data = surv_df_SCLC_subtype_WTS_meta[surv_df_SCLC_subtype_WTS_meta$ini_stage=="ED",])
cox_fit = coxph(Surv(OS_time, OS_event==1)~Heatmap_cluster, data = surv_df_SCLC_subtype_WTS_meta[surv_df_SCLC_subtype_WTS_meta$ini_stage=="ED",])
pdf(paste0(work_dir, "Fig_S3E_WTS_clustering_classification_survplot.pdf"), width = 7, height = 5.5)
p = ggsurvplot(
	fit = OS_fit,
	size = 1.5,
	xlab = "Months since Treatment",
	ylab = "OS",
	surv.median.line = "hv",
	palette = "nejm",
	xlim=c(0,30),
	ylim=c(0,1),
	break.time.by=3,
	font.title=c("bold"), font.subtitle=c("italic"),
	legend.title="Strata",
	pval = F, pval.coord = c(1, 0.25),
	surv.scale = "percent",
	risk.table = TRUE, risk.table.height = 0.3, conf.int = F)
p$plot <- p$plot + thickness
print(p)
dev.off()



OS_fit_summary = data.frame(summary(OS_fit)$table)
rownames(OS_fit_summary) = gsub("=", "", rownames(OS_fit_summary))
OS_fit_summary$rownames = rownames(OS_fit_summary)


cox_summary = cbind(data.frame(summary(cox_fit)$coefficients), data.frame(summary(cox_fit)$conf.int))
cox_summary$rownames = rownames(cox_summary)


overall_summary = merge(OS_fit_summary, cox_summary, by="rownames", all=TRUE)

overall_summary = (data.frame(overall_summary$rownames,
                 paste0(round(overall_summary$median, 1), " (", round(overall_summary$X0.95LCL, 1), "-", round(overall_summary$X0.95UCL, 1), ")"),
                 paste0(round(overall_summary$exp.coef., 1), " (", round(overall_summary$lower..95, 1), "-", round(overall_summary$upper..95, 1),")"),
                 overall_summary$Pr...z..))
colnames(overall_summary) = c("group", "mOS", "HR", "P")
overall_summary$group = gsub(" ", "", overall_summary$group)
overall_summary$P = round(overall_summary$P, 3)

apply(overall_summary, 1, function(x) {paste0(x[seq(1,4)], collapse = ";")})








