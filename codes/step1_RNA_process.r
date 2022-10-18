getwd()
.libPaths("/home/jjg/tools/RprofileLibpath/")
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
git_path = "/BiO2/users/jjg/SCLC/I_SCLC/Paper/221011_Github_code/Github/"

Exp.count = read.table("./expr_matrix/Exp.count.txt", header = TRUE, sep = '\t')


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
library("xlsx")
work_dir = "//BiO2/users/jjg/SCLC/I_SCLC/Paper/221011_Github_code/Fig1//"
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

############################################
# write cpm matrix
############################################

gp <- rep(1,ncol(Exp.count))
Exp.count2 <- Exp.count
rownames(Exp.count2) <- toTable(org.Hs.egSYMBOL)$gene_id[match(rownames(Exp.count), toTable(org.Hs.egSYMBOL)$symbol)]
library(edgeR)
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
write.table(log2.CPM.count, "log2_CPM_n226.txt", sep = "\t", quote = FALSE)

# Not filtered CPM
y <- DGEList(counts=Exp.count,
             gene=rownames(Exp.count),
             group=gp)
y$genes$ids <- rownames(Exp.count2)
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ 0 + gp, data = y$samples)
y <- estimateDisp(y)
log2_CPM_n226_nonGeneFilter <- cpm(y, prior.count=2, log=TRUE)
write.table(log2_CPM_n226_nonGeneFilter, "log2_CPM_n226_nonGeneFilter.txt", sep = "\t", quote = FALSE)

SCLC_cpm = log2.CPM.count[,SCLC_subtype_WTS_meta$WTS_ID]
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
                 title = "P versus Neuro, control is Neuro",
                 labSize = 5.0)
ggsave(paste0(work_dir, "plot_N_versus_A_EnhancedVolcanoPlot.png"), plot_N_versus_A, width = 7, height = 7)
write.table(N_versus_A$table, paste0(work_dir, "plot_N_versus_A_DEG.txt"), col.names = TRUE, row.names = TRUE, sep = '\t', quote = FALSE)
saveRDS(N_versus_A, paste0(work_dir, "plot_N_versus_A_DEG.Rds"))


N_versus_A = N_versus_A$table
Filt_DEG = N_versus_A[N_versus_A$PValue<0.001, ]
Filt_DEG = Filt_DEG[rownames(Filt_DEG)%in%rownames(log2.CPM.count),]



SCLC_subtype_WTS_meta_NE = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$NE_subtype=="NE", ]

NE_log2.CPM.count2 = log2.CPM.count[rownames(Filt_DEG),SCLC_subtype_WTS_meta_NE$WTS_ID]


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



pca_df <- prcomp(t(NE_log2.CPM.count2),
                 center = T,
                 scale. = T)


Fig1D_PCA_df = data.frame(pca_df$x[,"PC1"], pca_df$x[, "PC2"], pca_df$x[,"PC3"], SCLC_subtype_WTS_meta_NE$SCLC_subtype)
colnames(Fig1D_PCA_df) = c("PC1", "PC2", "PC3", "SCLC_subtype")


p = ggscatter(Fig1D_PCA_df[Fig1D_PCA_df$PC1<50,], x="PC1", y = "PC2", shape = "SCLC_subtype", color = "SCLC_subtype", palette = c("#BC3C29", "#0072B5", "#E18727"), ellipse=TRUE,
         conf.int = FALSE, cor.coef = FALSE, ellipse.type = "convex")
ggsave(paste0(work_dir, "Fig1D_PCA_plot.pdf"), plot = p, width = 5, height = 5)

calculate_group_centroid(Fig1D_PCA_df[, "PC1"], Fig1D_PCA_df[,"PC2"], SCLC_subtype_WTS_meta_NE$SCLC_subtype)


##### Fig 1E right


DEG_dataframe = data.frame()

for(i_subtype in c("A", "trans", "N")){
        sample_group_info = as.character(SCLC_subtype_WTS_meta$SCLC_subtype)
        set.seed(100)
        sample_group_info[!(seq(1, length(sample_group_info)) %in% c(sample(seq(1,length(sample_group_info))[sample_group_info=="A"], 15),
                          sample(seq(1,length(sample_group_info))[sample_group_info=="trans"], 15),
                          sample(seq(1,length(sample_group_info))[sample_group_info=="N"], 15)))] = "except"

        sample_group_info[!(sample_group_info%in%c("A", "trans", "N"))] = "except"
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


identical(colnames(SCLC_cpm), SCLC_subtype_WTS_meta$WTS_ID)
sample_group_info = as.character(SCLC_subtype_WTS_meta$SCLC_subtype)
print(table(sample_group_info))
SCLC_cpm2 = SCLC_cpm

######### return value

AverageExpr_DEG = lapply(DEG_dataframe$Gene[DEG_dataframe$Gene %in% rownames(SCLC_cpm2)], function(x) {
                x = as.character(x[1])
                expr_vector = SCLC_cpm2[x,] %>% unlist() %>% as.numeric()
               mean_A_type = mean(expr_vector[sample_group_info=="A"])
               mean_trans_type = mean(expr_vector[sample_group_info=="trans"])
               mean_N_type = mean(expr_vector[sample_group_info=="N"])
               return(c(mean_A_type, mean_trans_type, mean_N_type))
}) %>% data.frame() %>% t()
colnames(AverageExpr_DEG) = c("A", "trans", "N")

rownames(AverageExpr_DEG) = DEG_dataframe$Gene[DEG_dataframe$Gene %in% rownames(SCLC_cpm2)]
dend = as.dendrogram(hclust(dist(cor(AverageExpr_DEG))))
dd.reorder <- reorder(dend, c(3,1,2))

col_fun = colorRamp2(c(0.8, 0.9,1), c("#2E4A89", "#56CA58", "#FEFA2B"))
print(cor(AverageExpr_DEG))
#pdf(paste0(work_dir, "Fig1E_right_NeuroEndocrine_AverageExpr_DEG_correlation_heatmap.pdf"), width = 8, height = 7)
png(paste0(work_dir, "Fig1E_right_NeuroEndocrine_AverageExpr_DEG_correlation_heatmap.png"), width = 800, height = 700)
Heatmap(cor(AverageExpr_DEG), cluster_columns=rev(dd.reorder), cluster_rows=rev(dd.reorder), col = col_fun, name = "FC>1.25 & P value top50")
#Heatmap(SCLC_cpm[,sample_group_info!="except"])
dev.off()




DEG_dataframe = data.frame()

for(i_subtype in c("A", "trans","N", "P")){
        sample_group_info = as.character(SCLC_subtype_WTS_meta$SCLC_subtype)
        set.seed(100)
        sample_group_info[!(seq(1, length(sample_group_info)) %in% c(sample(seq(1,length(sample_group_info))[sample_group_info=="A"], 15),
                          sample(seq(1,length(sample_group_info))[sample_group_info=="trans"], 15),
                          sample(seq(1,length(sample_group_info))[sample_group_info=="N"], 15),
                          sample(seq(1,length(sample_group_info))[sample_group_info=="P"], 15)))] = "except"

        sample_group_info[!(sample_group_info%in%c("A", "trans", "N", "P"))] = "except"
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




identical(colnames(SCLC_cpm), SCLC_subtype_WTS_meta$WTS_ID)
sample_group_info = as.character(SCLC_subtype_WTS_meta$SCLC_subtype)
print(table(sample_group_info))
SCLC_cpm2 = SCLC_cpm

######### return value

AverageExpr_DEG = lapply(DEG_dataframe$Gene[DEG_dataframe$Gene %in% rownames(SCLC_cpm2)], function(x) {
                x = as.character(x[1])
                expr_vector = SCLC_cpm2[x,] %>% unlist() %>% as.numeric()
               mean_A_type = mean(expr_vector[sample_group_info=="A"])
               mean_trans_type = mean(expr_vector[sample_group_info=="trans"])
               mean_N_type = mean(expr_vector[sample_group_info=="N"])
               mean_P_type = mean(expr_vector[sample_group_info=="P"])
               mean_TN_type = mean(expr_vector[sample_group_info=="TN"])
               return(c(mean_A_type, mean_trans_type, mean_N_type, mean_P_type, mean_TN_type))
}) %>% data.frame() %>% t()
colnames(AverageExpr_DEG) = c("A", "trans", "N", "P", "TN")

rownames(AverageExpr_DEG) = DEG_dataframe$Gene[DEG_dataframe$Gene %in% rownames(SCLC_cpm2)]
dend = as.dendrogram(hclust(dist(cor(AverageExpr_DEG))))
dd.reorder <- reorder(dend, c(2, 1 , 3, 4))

col_fun = colorRamp2(c(0.6, 0.8,1), c("#2E4A89", "#56CA58", "#FEFA2B"))
#pdf(paste0(work_dir, "Fig1E_left_NeuroEndo_and_P_AverageExpr_DEG_correlation_heatmap.pdf"), width = 8, height = 8)
set.seed(0)
png(paste0(work_dir, "Fig1E_left_NeuroEndo_and_P_AverageExpr_DEG_correlation_heatmap.png"), width = 800, height = 800)
Heatmap(cor(AverageExpr_DEG), cluster_columns=TRUE, cluster_rows=TRUE, col = col_fun, name = "FC > 1.5")
#Heatmap(SCLC_cpm[,sample_group_info!="except"])
dev.off()

########################################
## N versus A&trans
############################

sample_group_info = lapply(as.character(SCLC_subtype_WTS_meta$SCLC_subtype), function(x) {if(x=="A" | x=="trans"){
               return("A&trans")
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

N_versus_AandTrans <- exactTest(y, pair=c("A&trans","N")) # compare groups 1 and 2
topTags(N_versus_AandTrans, n=10)

target_gene_of_N_vs_Atrans = c("ASCL1", "CD74", "ISG20", "CD38", "STAT4", "CXCL10", "CXCL9", "IL2RA", "CIITA", "CCL18", "CXCL13",
                                            "CLKNK", "CLECL1", "MCTP2", "COTL1", "ETS2", "SEC11C",
                                              "NOTCH2", "COL9A1", "NEUROD1", "MMP8", "NOTCH1", "FGFR3","SERPINA10", "ITGB3", "TIMP3", "MMP3")
N_Atrans_DEG = N_versus_AandTrans$table
N_Atrans_DEG$gene = rownames(N_Atrans_DEG)

keyvals = apply(N_Atrans_DEG, 1, function(x) {
                        if(as.numeric(x[3])< 0.05 & abs(as.numeric(x[1])) >= 1){
                                if(as.numeric(x[1]) >= 0.1){
                                        if(sum(x[4]%in%target_gene_of_N_vs_Atrans)>0){
                                                return("#1e2f97")
                                        }else{return("#949494")}
                                }else if(as.numeric(x[1]) <= -0.1){
                                        if(sum(x[4]%in%target_gene_of_N_vs_Atrans)>0){
                                                return("#800000")
                                        }else{return("#949494")}
                                }
                        }else{return("#949494")}
                                              }) %>% unlist() %>% as.character()

names(keyvals)[keyvals == '#800000'] <- 'A&trans Immune'
#names(keyvals)[keyvals == '#c38e63'] <- 'nonhigh'# "#f6bdc0"
names(keyvals)[keyvals == "#1e2f97"] <- 'N Notch'
#names(keyvals)[keyvals == '#59788e'] <- 'nonlow' #"#b0dbf1"
names(keyvals)[keyvals == '#949494'] <- 'non-sig'


plot_N_Atrans_DEG = EnhancedVolcano(N_Atrans_DEG,
                lab = rownames(N_Atrans_DEG),
                 x = 'logFC', labSize=3,
                 y = 'PValue', xlim = c(-5.6, 5.6), ylim = c(0, 20),
                 pCutoff = 0.05, colCustom=keyvals,
                 pointSize = 1, selectLab=target_gene_of_N_vs_Atrans,
                 drawConnectors = TRUE,directionConnectors = "both",  max.overlaps=Inf, colAlpha = 0.9)
ggsave(paste0(work_dir, "Fig_S5B_edgeR_N_Atrans_DEG.png"), plot = plot_N_Atrans_DEG, width = 5, height = 6)

saveRDS(N_versus_AandTrans, paste0(work_dir, "Fig_S5B_edgeR_N_Atrans_DEG.Rds"))
########################################
## P versus A&trans
############################

sample_group_info = lapply(as.character(SCLC_subtype_WTS_meta$SCLC_subtype), function(x) {if(x=="A" | x=="trans"){
               return("A&trans")
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

P_versus_AandTrans <- exactTest(y, pair=c("A&trans","P")) # compare groups 1 and 2
topTags(P_versus_AandTrans, n=10)


target_gene_of_P_vs_Atrans = c("ASCL1", "POU2F3", "NOTCH1", "NOTCH2",
                                            "HLA-DQA1", "IL7R", "IRF8", "CD4", "LYN", "IL17RB", "CD8A", "IL1R2", "PTGER4", "CCR3", "CXCL3", "CXCL2",
                                            "SAMHD1", "PIK3R5", "RNF166", "SPIC", "SPIB")

NE_25_genelist = c("BEX1","ASCL1","INSM1","CHGA","TAGLN3","KIF5C","CRMP1","SCG3","SYT4","RTN1","MYT1","SYP","KIF1A","TMSB15A","SYN1","SYT11","RUNDC3A","TFF3","CHGB","FAM57B","SH3GL2","BSN","SEZ6","TMSB15B","CELF3")

P_Atrans_DEG = P_versus_AandTrans$table
P_Atrans_DEG$gene = rownames(P_Atrans_DEG)
keyvals = apply(P_Atrans_DEG, 1, function(x) {
                        if(as.numeric(x[3])< 0.05 & abs(as.numeric(x[1])) >= 1){
                                if(as.numeric(x[1]) >= 0.1){
                                        if(sum(x[4]%in%target_gene_of_P_vs_Atrans)>0){
                                                return('#800000')
                                        }else{return("#949494")}
                                }else if(as.numeric(x[1]) <= -0.1){
                                        if(sum(x[4]%in%c(target_gene_of_P_vs_Atrans, NE_25_genelist))>0){
                                                return('black')
                                        }else{return("#949494")}
                                }
                        }else{return("#949494")}
                                              }) %>% unlist() %>% as.character()

names(keyvals)[keyvals == 'black'] <- 'A&trans NE'
#names(keyvals)[keyvals == '#c38e63'] <- 'nonhigh'# "#f6bdc0"
names(keyvals)[keyvals == '#800000'] <- 'P immune'
#names(keyvals)[keyvals == '#59788e'] <- 'nonlow' #"#b0dbf1"
names(keyvals)[keyvals == '#949494'] <- 'non-sig'

VolcanoPlot_genelist = c(target_gene_of_P_vs_Atrans, NE_25_genelist)
VolcanoPlot_genelist = VolcanoPlot_genelist[VolcanoPlot_genelist%in%rownames(P_Atrans_DEG)]
plot_P_Atrans_DEG = EnhancedVolcano(P_Atrans_DEG,
                lab = rownames(P_Atrans_DEG),
                 x = 'logFC', labSize=3,
                 y = 'PValue', ylim = c(0, 40),#xlim = c(-5.6, 5.6), ylim = c(0, 20),
                 pCutoff = 0.05,
                 pointSize = 1, selectLab=VolcanoPlot_genelist, colCustom=keyvals,
                 drawConnectors = TRUE,directionConnectors = "both",  max.overlaps=Inf, colAlpha = 0.9)
ggsave(paste0(work_dir, "Fig_S5A_edgeR_P_Atrans_DEG.png"), plot = plot_P_Atrans_DEG, width = 5, height = 6)
saveRDS(P_versus_AandTrans, paste0(work_dir, "Fig_S5A_edgeR_P_Atrans_DEG.Rds"))


############################### WTS subype
thickness = theme(legend.title = element_text(size = 13, color = "black", face = "bold"),
      legend.text = element_text(size = 13, color = "black", face = "bold"),
      axis.text.x = element_text(size = 13, color = "black", face = "bold"),
      axis.text.y = element_text(size = 13, color = "black", face = "bold"),
      axis.title.x = element_text(size = 13, color = "black", face = "bold"),
      axis.title.y = element_text(size = 13, color = "black", face = "bold"))



SCLC_subtype_WTS_meta = SCLC_meta[SCLC_meta$WTS_QC_Result %in% c("Pass", "1") & grepl("includ", SCLC_meta$study.inclusion),]# & !is.na(SCLC_meta$SCLC_subtype),]

rownames(SCLC_subtype_WTS_meta) = SCLC_subtype_WTS_meta$WTS_ID
SCLC_subtype_WTS_meta$SCLC_subtype[is.na(SCLC_subtype_WTS_meta$SCLC_subtype)] = "non-IHC"
SCLC_subtype_WTS_meta$SCLC_subtype = factor(SCLC_subtype_WTS_meta$SCLC_subtype, levels = c("A", "trans", "N", "P", "TN", "non-IHC"))
SCLC_cpm = log2.CPM.count[,SCLC_subtype_WTS_meta$WTS_ID]


git_path



library(singscore)

SCLC_cpm_Rank <- rankGenes(SCLC_cpm)

SingScore <- simpleScore(
  rankData = SCLC_cpm_Rank,
  upSet = NE_25_genelist,
  centerScore = T,
  knownDirection = T
  )
SCLC_subtype_WTS_meta$NE_25_genelist = SingScore$TotalScore

###########################
sig_dir = "/BiO2/users/jjg/SCLC/I_SCLC/220126_Borromeo_et_al/Signature_list_adpval_0.001/"
sig_file_list = list.files(sig_dir, pattern = ".txt")


pathway_list = list()
for(i_file in sig_file_list){
        if(strsplit(i_file, "_")[[1]][1]=="Pathset"){
                signature_df = read.table(paste0(sig_dir, i_file), sep = '\t', header = TRUE)
                for(i_row in seq(1,nrow(signature_df))){
                        signature_name = signature_df[i_row,"pathway"]
                        signature_name = gsub(" ","_",signature_name)
                        pathway_list[[paste0(gsub("\\.txt", "", i_file), ";", signature_name)]] = strsplit(gsub(" ", "", signature_df[i_row,"Genes"]), ";")[[1]]
                }

        }else{
                signature_df = read.table(paste0(sig_dir, i_file), sep = '\t', header = TRUE)
                signature_name = strsplit(i_file, "\\.")[[1]][1]
                pathway_list[[paste0(gsub("\\.txt", "", i_file),";", signature_name)]] = signature_df$Symbol
        }
}

for(i_name in c("ASCL1_ND1_Shared_targets;ASCL1_ND1_Shared_targets", "ASCL1_high_signatures_cell_line;ASCL1_high_signatures_cell_line",
                "ND1_high_signatures_cell_line;ND1_high_signatures_cell_line")){
        gene_list = pathway_list[i_name] %>% unlist() %>% as.character() %>% unique()
            SingScore <- simpleScore(
          rankData = SCLC_cpm_Rank,
          upSet = gene_list,
          centerScore = T,
          knownDirection = T
          )
        i_name = strsplit(i_name, ";")[[1]][1]
        SCLC_subtype_WTS_meta[,i_name] = SingScore$TotalScore
}

####################### TUFT cell marker, from all gene cpm


SCLC_cpm_raw = read.table("./log2_CPM_n226_nonGeneFilter.txt")[,SCLC_subtype_WTS_meta$WTS_ID]
SCLC_cpm_Rank_raw <- rankGenes(SCLC_cpm_raw)

SingScore <- simpleScore(
  rankData = SCLC_cpm_Rank_raw,
  upSet = c("POU2F3", "TRPM5", "SOX9", "GFI1B", "CHAT", "ASCL2", "AVIL"),
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$Tuft_cell_marker = SingScore$TotalScore

##################### zscaling
SCLC_subtype_WTS_meta$NEUROD1_cpm = as.numeric(SCLC_cpm["NEUROD1",])
SCLC_subtype_WTS_meta$NEUROD1_cpm_z = (SCLC_subtype_WTS_meta$NEUROD1_cpm-mean(SCLC_subtype_WTS_meta$NEUROD1_cpm))/sd(SCLC_subtype_WTS_meta$NEUROD1_cpm)
SCLC_subtype_WTS_meta$ASCL1_cpm = as.numeric(SCLC_cpm["ASCL1",])
SCLC_subtype_WTS_meta$ASCL1_cpm_z = (SCLC_subtype_WTS_meta$ASCL1_cpm-mean(SCLC_subtype_WTS_meta$ASCL1_cpm))/sd(SCLC_subtype_WTS_meta$ASCL1_cpm)
SCLC_subtype_WTS_meta$POU2F3_cpm = as.numeric(SCLC_cpm["POU2F3",])
SCLC_subtype_WTS_meta$POU2F3_cpm_z = (SCLC_subtype_WTS_meta$POU2F3_cpm-mean(SCLC_subtype_WTS_meta$POU2F3_cpm))/sd(SCLC_subtype_WTS_meta$POU2F3_cpm)


SCLC_subtype_WTS_meta$Tuft_cell_marker_z = (SCLC_subtype_WTS_meta$Tuft_cell_marker-mean(SCLC_subtype_WTS_meta$Tuft_cell_marker))/sd(SCLC_subtype_WTS_meta$Tuft_cell_marker)
SCLC_subtype_WTS_meta$ASCL1_high_signatures_z = (SCLC_subtype_WTS_meta$ASCL1_high_signatures_cell_line-mean(SCLC_subtype_WTS_meta$ASCL1_high_signatures_cell_line))/sd(SCLC_subtype_WTS_meta$ASCL1_high_signatures_cell_line)

SCLC_subtype_WTS_meta$ND1_high_signatures_z = (SCLC_subtype_WTS_meta$ND1_high_signatures_cell_line-mean(SCLC_subtype_WTS_meta$ND1_high_signatures_cell_line))/sd(SCLC_subtype_WTS_meta$ND1_high_signatures_cell_line)

SCLC_subtype_WTS_meta$NE_25_genelist_z = (SCLC_subtype_WTS_meta$NE_25_genelist-mean(SCLC_subtype_WTS_meta$NE_25_genelist))/sd(SCLC_subtype_WTS_meta$NE_25_genelist)

##################### call WTS subtype

P_subtype = SCLC_subtype_WTS_meta$WTS_ID[SCLC_subtype_WTS_meta$Tuft_cell_marker_z > 1]# & SCLC_subtype_WTS_meta$POU2F3_cpm_z > 0.5]
TN_subtype = SCLC_subtype_WTS_meta$WTS_ID[SCLC_subtype_WTS_meta$Tuft_cell_marker_z < -0.25 &
                                          SCLC_subtype_WTS_meta$ASCL1_high_signatures_z < -0.25 &
                                          SCLC_subtype_WTS_meta$ND1_high_signatures_z < -0.25 &
                                          !(SCLC_subtype_WTS_meta$WTS_ID%in%P_subtype)]
N_subtype = SCLC_subtype_WTS_meta$WTS_ID[SCLC_subtype_WTS_meta$ND1_high_signatures_z > 0.5 & SCLC_subtype_WTS_meta$ASCL1_high_signatures_z < -0.5 &
                                         SCLC_subtype_WTS_meta$NEUROD1_cpm_z > 0.5 & SCLC_subtype_WTS_meta$ASCL1_cpm_z < -0.5]
N_subtype = N_subtype[!(N_subtype %in% c(P_subtype, TN_subtype))]
A_subtype = SCLC_subtype_WTS_meta$WTS_ID[SCLC_subtype_WTS_meta$ND1_high_signatures_z < 0]#&#quantile(SCLC_subtype_WTS_meta$ND1_high_signatures_z, 0.75) &
A_subtype = A_subtype[!(A_subtype %in% c(P_subtype, TN_subtype, N_subtype))]



################################ Fig 2B
require(survminer)
require(survival)



WTS_jjg_subtype = rep("trans", nrow(SCLC_subtype_WTS_meta))
WTS_jjg_subtype[SCLC_subtype_WTS_meta$WTS_ID %in% P_subtype] = "P"
WTS_jjg_subtype[SCLC_subtype_WTS_meta$WTS_ID %in% TN_subtype] = "TN"
WTS_jjg_subtype[SCLC_subtype_WTS_meta$WTS_ID %in% N_subtype] = "N"
WTS_jjg_subtype[SCLC_subtype_WTS_meta$WTS_ID %in% A_subtype] = "A"

identical(WTS_jjg_subtype, SCLC_subtype_WTS_meta$WTS_subtype)
SCLC_subtype_WTS_meta$WTS_subtype = factor(SCLC_subtype_WTS_meta$WTS_subtype, levels = c("A", "trans", "N", "P", "TN"))


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


##### 4 * 4, Fig S1D Up left
ggconfusion_P = SCLC_subtype_WTS_meta[, c("WTS_subtype", "SCLC_subtype")]
ggconfusion_P = ggconfusion_P[ggconfusion_P$SCLC_subtype!="non-IHC",]

ggconfusion_P$WTS_subtype = as.character(ggconfusion_P$WTS_subtype)
ggconfusion_P$WTS_subtype[ggconfusion_P$WTS_subtype %in% c("A", "trans")] = "A&trans"
ggconfusion_P$SCLC_subtype = as.character(ggconfusion_P$SCLC_subtype)
ggconfusion_P$SCLC_subtype[ggconfusion_P$SCLC_subtype %in% c("A", "trans")] = "A&trans"

cfm3 <- confusionMatrix(factor(ggconfusion_P$WTS_subtype, levels = c("A&trans", "N", "P", "TN")),
                        factor(ggconfusion_P$SCLC_subtype, levels = c("A&trans", "N", "P", "TN")))
ggcfm3 <- ggplotConfusionMatrix(cfm3)
ggsave(paste0(work_dir, "Fig_S1D_Upleft_WTS_subtype_4.4_cfm.png"), plot = ggcfm3, width = 9, height = 7)

###### 3 * 3, Fig S1E Down left
ggconfusion_no_TN = SCLC_subtype_WTS_meta[, c("WTS_subtype", "SCLC_subtype")]
ggconfusion_no_TN = ggconfusion_no_TN[ggconfusion_no_TN$SCLC_subtype!="non-IHC",]
ggconfusion_no_TN = ggconfusion_no_TN[ggconfusion_no_TN$WTS_subtype!="TN",]

ggconfusion_no_TN$WTS_subtype = as.character(ggconfusion_no_TN$WTS_subtype)
ggconfusion_no_TN$WTS_subtype[ggconfusion_no_TN$WTS_subtype %in% c("A", "trans")] = "A&trans"
ggconfusion_no_TN$SCLC_subtype = as.character(ggconfusion_no_TN$SCLC_subtype)
ggconfusion_no_TN$SCLC_subtype[ggconfusion_no_TN$SCLC_subtype %in% c("A", "trans")] = "A&trans"

cfm3 <- confusionMatrix(factor(ggconfusion_no_TN$WTS_subtype, levels = c("A&trans", "N", "P")),
                        factor(ggconfusion_no_TN$SCLC_subtype, levels = c("A&trans", "N", "P")))
ggcfm3 <- ggplotConfusionMatrix(cfm3)
ggsave(paste0(work_dir, "Fig_S1D_Downleft_WTS_subtype_3.3cfm.png"), plot = ggcfm3, width = 9, height = 7)


##### Calculate kappa
ggconfusion_P = SCLC_subtype_WTS_meta[, c("WTS_subtype", "SCLC_subtype")]
ggconfusion_P = ggconfusion_P[ggconfusion_P$SCLC_subtype!="non-IHC",]
cfm3 <- confusionMatrix(factor(ggconfusion_P$WTS_subtype, levels = c("A", "trans", "N", "P", "TN")),
                        factor(ggconfusion_P$SCLC_subtype, levels = c("A", "trans", "N", "P", "TN")))
ggcfm3 <- ggplotConfusionMatrix(cfm3)
ggconfusion_P = ggconfusion_P[ggconfusion_P$WTS_subtype!="TN",]
kappa_df = ggconfusion_P
kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
kappa_df$WTS_subtype[kappa_df$WTS_subtype!="P"] = "non-P"
kappa_df$SCLC_subtype = as.character(kappa_df$SCLC_subtype)
kappa_df$SCLC_subtype[kappa_df$SCLC_subtype!="P"] = "non-P"
cfm3 <- confusionMatrix(factor(kappa_df$WTS_subtype, levels = c("P", "non-P")),
                        factor(kappa_df$SCLC_subtype, levels = c("P", "non-P")))
paste("P Accuracy", percent_format()(cfm3$overall[1]), "Kappa", percent_format()(cfm3$overall[2]))

kappa_df = ggconfusion_P[ggconfusion_P$WTS_subtype!="P", ]
kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
kappa_df$WTS_subtype[kappa_df$WTS_subtype!="N"] = "non-N"
kappa_df$SCLC_subtype = as.character(kappa_df$SCLC_subtype)
kappa_df$SCLC_subtype[kappa_df$SCLC_subtype!="N"] = "non-N"
cfm3 <- confusionMatrix(factor(kappa_df$WTS_subtype, levels = c("N", "non-N")),
                        factor(kappa_df$SCLC_subtype, levels = c("N", "non-N")))
paste("N Accuracy", percent_format()(cfm3$overall[1]), "Kappa", percent_format()(cfm3$overall[2]))


kappa_df = ggconfusion_P
kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
kappa_df$WTS_subtype[kappa_df$WTS_subtype %in% c("A", "trans")] = "A&trans"
kappa_df$WTS_subtype[kappa_df$WTS_subtype %in% c("N", "P")] = "N&P"
kappa_df$SCLC_subtype = as.character(kappa_df$SCLC_subtype)
kappa_df$SCLC_subtype[kappa_df$SCLC_subtype %in% c("A", "trans")] = "A&trans"
kappa_df$SCLC_subtype[kappa_df$SCLC_subtype %in% c("N", "P")] = "N&P"


kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
cfm3 <- confusionMatrix(factor(kappa_df$WTS_subtype, levels = c("A&trans", "N&P")),
                        factor(kappa_df$SCLC_subtype, levels = c("A&trans", "N&P")))
paste("A trans Accuracy", percent_format()(cfm3$overall[1]), "Kappa", percent_format()(cfm3$overall[2]))


###################
##### Fig S1D right
ggconfusion_P = SCLC_subtype_WTS_meta[, c("WTS_subtype", "SCLC_subtype")]
ggconfusion_P = ggconfusion_P[ggconfusion_P$SCLC_subtype!="non-IHC",]

ggconfusion_P$WTS_subtype = as.character(ggconfusion_P$WTS_subtype)
ggconfusion_P$WTS_subtype[ggconfusion_P$WTS_subtype %in% c("N", "trans")] = "N&trans"
ggconfusion_P$SCLC_subtype = as.character(ggconfusion_P$SCLC_subtype)
ggconfusion_P$SCLC_subtype[ggconfusion_P$SCLC_subtype %in% c("N", "trans")] = "N&trans"

cfm3 <- confusionMatrix(factor(ggconfusion_P$WTS_subtype, levels = c("A", "N&trans", "P", "TN")),
                        factor(ggconfusion_P$SCLC_subtype, levels = c("A", "N&trans", "P", "TN")))
ggcfm3 <- ggplotConfusionMatrix(cfm3)
ggsave(paste0(work_dir, "Fig_S1D_Upright_WTS_subtype_4.4_cfm.png"), plot = ggcfm3, width = 7, height = 7)


ggconfusion_no_TN = SCLC_subtype_WTS_meta[, c("WTS_subtype", "SCLC_subtype")]
ggconfusion_no_TN = ggconfusion_no_TN[ggconfusion_no_TN$SCLC_subtype!="non-IHC",]
ggconfusion_no_TN = ggconfusion_no_TN[ggconfusion_no_TN$WTS_subtype!="TN",]

ggconfusion_no_TN$WTS_subtype = as.character(ggconfusion_no_TN$WTS_subtype)
ggconfusion_no_TN$WTS_subtype[ggconfusion_no_TN$WTS_subtype %in% c("N", "trans")] = "N&trans"
ggconfusion_no_TN$SCLC_subtype = as.character(ggconfusion_no_TN$SCLC_subtype)
ggconfusion_no_TN$SCLC_subtype[ggconfusion_no_TN$SCLC_subtype %in% c("N", "trans")] = "N&trans"

cfm3 <- confusionMatrix(factor(ggconfusion_no_TN$WTS_subtype, levels = c("A", "N&trans", "P")),
                        factor(ggconfusion_no_TN$SCLC_subtype, levels = c("A", "N&trans", "P")))
ggcfm3 <- ggplotConfusionMatrix(cfm3)
ggsave(paste0(work_dir, "Fig_S1D_Downright_WTS_subtype_3.3cfm.png"), plot = ggcfm3, width = 7, height = 7)




ggconfusion_P = SCLC_subtype_WTS_meta[, c("WTS_subtype", "SCLC_subtype")]
ggconfusion_P = ggconfusion_P[ggconfusion_P$SCLC_subtype!="non-IHC",]
ggconfusion_P = ggconfusion_P[ggconfusion_P$WTS_subtype!="TN",]

kappa_df = ggconfusion_P
kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
kappa_df$WTS_subtype[kappa_df$WTS_subtype!="P"] = "non-P"
kappa_df$SCLC_subtype = as.character(kappa_df$SCLC_subtype)
kappa_df$SCLC_subtype[kappa_df$SCLC_subtype!="P"] = "non-P"
cfm3 <- confusionMatrix(factor(kappa_df$WTS_subtype, levels = c("P", "non-P")),
                        factor(kappa_df$SCLC_subtype, levels = c("P", "non-P")))
paste("P Accuracy", percent_format()(cfm3$overall[1]), "Kappa", percent_format()(cfm3$overall[2]))

kappa_df = ggconfusion_P[ggconfusion_P$WTS_subtype!="P", ]
kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
kappa_df$WTS_subtype[kappa_df$WTS_subtype %in% c("N", "trans")] = "N&trans"
kappa_df$SCLC_subtype = as.character(kappa_df$SCLC_subtype)
kappa_df$SCLC_subtype[kappa_df$SCLC_subtype %in% c("N", "trans")] = "N&trans"


cfm3 <- confusionMatrix(factor(kappa_df$WTS_subtype, levels = c("N&trans", "A")),
                        factor(kappa_df$SCLC_subtype, levels = c("N&trans", "A")))
paste("N&trans Accuracy", percent_format()(cfm3$overall[1]), "Kappa", percent_format()(cfm3$overall[2]))


kappa_df = ggconfusion_P
kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
kappa_df$WTS_subtype[kappa_df$WTS_subtype %in% c("trans", "N", "P")] = "non-A"

kappa_df$SCLC_subtype = as.character(kappa_df$SCLC_subtype)
kappa_df$SCLC_subtype[kappa_df$SCLC_subtype %in% c("trans", "N", "P")] = "non-A"


kappa_df$WTS_subtype = as.character(kappa_df$WTS_subtype)
cfm3 <- confusionMatrix(factor(kappa_df$WTS_subtype, levels = c("non-A", "A")),
                        factor(kappa_df$SCLC_subtype, levels = c("non-A", "A")))
paste("A Accuracy", percent_format()(cfm3$overall[1]), "Kappa", percent_format()(cfm3$overall[2]))


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
colnames(overall_summary) = c("group", "mPFS", "HR", "P")
overall_summary$group = gsub(" ", "", overall_summary$group)
overall_summary$P = round(overall_summary$P, 3)

apply(overall_summary, 1, function(x) {paste0(x[seq(1,4)], collapse = ";")})

##### Fig S3E

SCLC_subtype_WTS_meta$ctx_1st_pfs = as.numeric(SCLC_subtype_WTS_meta$ctx_1st_pfs)
SCLC_subtype_WTS_meta$ctx_1st_pfs_event = as.numeric(SCLC_subtype_WTS_meta$ctx_1st_pfs_event)

Pfs_fit <- survfit(Surv(ctx_1st_pfs, ctx_1st_pfs_event==1)~WTS_subtype, data = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$ini_stage=="ED",])
cox_fit = coxph(Surv(ctx_1st_pfs, ctx_1st_pfs_event==1)~WTS_subtype, data = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$ini_stage=="ED",])

pdf(paste0(work_dir, "Fig_S3E.ctx_1st_pfs_WTS_subtype_ED.pdf"), width = 8, height = 6)
p = ggsurvplot(
  fit = Pfs_fit,
  size = 1.5,
  xlab = "Months since Treatment",
  ylab = "PFS",
  surv.median.line = "hv",
  palette = "nejm",
  xlim=c(0,8),
  ylim=c(0,1),
  break.time.by=2,
  title="Pfs Analysis in ED-patients: WTS subtype",
  font.title=c("bold"), font.subtitle=c("italic"), font.size=13,
  legend.title="Strata",
  pval = T, pval.coord = c(1, 0.25),
  surv.scale = "percent",
  risk.table = TRUE, risk.table.height = 0.3, conf.int = F, newpage = FALSE)
p$plot <- p$plot + thickness
print(p)
dev.off()


###############
##### Fig S3F

NE_SCLC_subtype_WTS_meta = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype %in% c("A", "trans", "N"),]
NE_SCLC_cpm = rbind(log2.CPM.count[c("NEUROD1", "ASCL1"),NE_SCLC_subtype_WTS_meta$WTS_ID], t(NE_SCLC_subtype_WTS_meta[, c("ASCL1_high_signatures_cell_line", "ND1_high_signatures_cell_line")]))

pt.matrix <- t(apply(NE_SCLC_cpm,1,function(x){(x-mean(x))/sd(x)}))
#pt.matrix <- apply(pt.matrix,2,function(x){(x-mean(x))/sd(x)})
#rownames(pt.matrix) <- rownames(NE_SCLC_cpm); colnames(pt.matrix) = colnames(NE_SCLC_cpm)
i_cluster_medhod = "ward.D2"
i_distance_method = "maximum"

k_cluster = 3
set.seed(100)
cpm_heatmap = Heatmap(pt.matrix, column_split = k_cluster, clustering_distance_columns = i_distance_method, show_column_names = FALSE, clustering_method_columns = i_cluster_medhod, column_gap = unit(5, "mm"))
cpm_heatmap_ht = draw(cpm_heatmap); cpm_heatmap_tr_order = column_order(cpm_heatmap_ht)
pdf(paste0(work_dir, "Fig_S3F_WTS_clustering_class_",k_cluster,"_",i_cluster_medhod,"_",i_distance_method,".pdf"), width = 12, height = 4)
print(cpm_heatmap_ht)
dev.off()
tmp = rep("cluster", ncol(NE_SCLC_cpm))
tmp[cpm_heatmap_tr_order[[1]]] = "1"; tmp[cpm_heatmap_tr_order[[2]]] = "2"; tmp[cpm_heatmap_tr_order[[3]]] = "3"
surv_df_SCLC_subtype_WTS_meta = NE_SCLC_subtype_WTS_meta
surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster = factor(tmp, levels = c("1", "2", "3"))

NEUROD1_index = which.max(c(mean(pt.matrix["NEUROD1",cpm_heatmap_tr_order[[1]]]), mean(pt.matrix["NEUROD1",cpm_heatmap_tr_order[[2]]]), mean(pt.matrix["NEUROD1",cpm_heatmap_tr_order[[3]]])))
ASCL1_index = which.max(c(mean(pt.matrix["ASCL1",cpm_heatmap_tr_order[[1]]]), mean(pt.matrix["ASCL1",cpm_heatmap_tr_order[[2]]]), mean(pt.matrix["ASCL1",cpm_heatmap_tr_order[[3]]])))

surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster = as.character(surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster)
surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster[surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster==as.character(NEUROD1_index)] = "N"
surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster[surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster==as.character(ASCL1_index)] = "A"
surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster[!surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster %in% c("A", "N")] = "trans"
surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster = factor(surv_df_SCLC_subtype_WTS_meta$Heatmap_cluster, levels = c("A", "trans", "N"))

OS_fit <- survfit(Surv(OS_time, OS_event==1)~Heatmap_cluster, data = surv_df_SCLC_subtype_WTS_meta[surv_df_SCLC_subtype_WTS_meta$ini_stage=="ED",])
cox_fit = coxph(Surv(OS_time, OS_event==1)~Heatmap_cluster, data = surv_df_SCLC_subtype_WTS_meta[surv_df_SCLC_subtype_WTS_meta$ini_stage=="ED",])

pdf(paste0(work_dir, "Fig_S3G_WTS_clustering_class_", i_cluster_medhod,"_",i_distance_method,"_K,",k_cluster, "_", pvalue, "_cluster_cpm_matrix_Survplot.pdf"), width = 8, height = 6.2)
        print(
                ggsurvplot(
                        fit = OS_fit,
                        size = 1.5,
                        xlab = "Months since Treatment",
                        ylab = "OS",
                        surv.median.line = "hv",
                        palette = "nejm",
                        xlim=c(0,30),
                        ylim=c(0,1),
                        break.time.by=3,
                        title=paste0("OS: By Heatmap Cluster:", k_cluster, " ", i_distance_method),
                        font.title=c("bold"), font.subtitle=c("italic"),
                        legend.title="Strata",
                        pval = T, pval.coord = c(1, 0.25),
                        surv.scale = "percent",
                        risk.table = TRUE, risk.table.height = 0.3, conf.int = F),
        surv.plot.height = NULL,
        risk.table.height = NULL,
        ncensor.plot.height = NULL,
        newpage = TRUE)
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
colnames(overall_summary) = c("group", "mPFS", "HR", "P")
overall_summary$group = gsub(" ", "", overall_summary$group)
overall_summary$P = round(overall_summary$P, 3)

apply(overall_summary, 1, function(x) {paste0(x[seq(1,4)], collapse = ";")})

