LibPath = paste0(getwd(), "/LibPath/")
.libPaths(LibPath)
library(maftools)
library(tidyverse)
library(data.table)
library(dplyr)
library(magrittr)
library(Hmisc)
library(EnhancedVolcano)
library(ggpubr)
library('xlsx')
source('./functions/custom_plotCE_enrichment.r')
source('./functions/Custom_gap_oncoplot.R')
source("./functions/Custom_oncomatrix.R")
source("./functions/Custom_summarizeMaf.R")
source("./functions/Custom_subsetMaf.R")



plot_dir = "../Figures/"

SCLC_meta = read.xlsx("../resource/Essential_check_Edited.xlsx", header = TRUE, sheetName = "Sheet1")

CNV = read.table("../resource/CS_CNV_matrix.txt", sep = "\t", header = TRUE)
maf.filename = "../resource/annovar_input.maf"
temp.maf<-read.maf(maf=maf.filename, cnTable =CNV )
matched_meta_info = SCLC_meta[as.numeric(lapply(unique(as.character(temp.maf@data$Tumor_Sample_Barcode)), function(x) {grep(x, SCLC_meta$SGI_ID)})),]
matched_meta_info = data.frame(unique(as.character(temp.maf@data$Tumor_Sample_Barcode)), matched_meta_info)
colnames(matched_meta_info) = c("Tumor_Sample_Barcode", colnames(SCLC_meta))
tmp = as.character(lapply(matched_meta_info$SGI_ID, function(x) {strsplit(gsub("^ ", "", x), " ")[[1]][1]}))
matched_meta_info$Tumor_Sample_Barcode = tmp
rownames(matched_meta_info) = tmp
matched_meta_info$AAN_N_P = lapply(matched_meta_info$SCLC_subtype, function(x) {if(is.na(x)){return("NA")}else if(x=="A" | x=="AN"){return("A&AN")}else{return(x)}}) %>% unlist() %>% as.character()

temp.maf<-read.maf(maf=maf.filename, cnTable =CNV, clinicalData = matched_meta_info)


matched_meta_info = matched_meta_info[!is.na(matched_meta_info$SCLC_subtype) & matched_meta_info$CS_QC == "Pass" & grepl("include", matched_meta_info$study.inclusion),]

laml = subsetMaf(temp.maf, query = 'Tumor_Sample_Barcode %in% matched_meta_info$Tumor_Sample_Barcode', mafObj = TRUE, fields = colnames(data.frame(temp.maf@data))[seq(1,99)])

#####################################################
axis_theme = theme(axis.text.x = element_text(size=10, face="bold", color = "black"),
              axis.title.x = element_text(size=10, face="bold", color = "black"),
        axis.text.y = element_text(size=10, face="bold", color = "black"),
        axis.title.y = element_text(size=10, face="bold", color = "black"), panel.border = element_rect(size = 1))
#####################################################


oncogenic_pathway = data.table::fread(system.file("extdata", "oncogenic_sig_patwhays.tsv",
            package = "maftools")) %>% data.frame()


SGI_results_of_pathways = list()
plot_maf_subset = data.frame(laml@data)
for(search_signature in unique(oncogenic_pathway$Pathway)){
        search_signature_sample_match = plot_maf_subset$Hugo_Symbol %in% as.character(oncogenic_pathway[oncogenic_pathway$Pathway==search_signature, "Gene"]) %>%
                plot_maf_subset[., "Tumor_Sample_Barcode"] %>% unique() %>% as.character()
        search_signature_gene_match = as.character(oncogenic_pathway[oncogenic_pathway$Pathway==search_signature, "Gene"]) %in% plot_maf_subset$Hugo_Symbol
        SGI_results_of_pathways[[search_signature]] = c(length(search_signature_sample_match), length(unique(laml@variant.type.summary$Tumor_Sample_Barcode)),
                                                        sum(search_signature_gene_match), length(search_signature_gene_match))
}

SGI_results_of_pathways = t(data.frame(SGI_results_of_pathways)); colnames(SGI_results_of_pathways) = c("matched_samples", "total_samples", "matched_genes", "total_pathway_genes")

signature_order_by_sample = rownames(SGI_results_of_pathways)[rev(order(apply(SGI_results_of_pathways, 1, function(x) {x[1]/x[2]})))] %>% rev()
signature_order_by_gene = rownames(SGI_results_of_pathways)[rev(order(apply(SGI_results_of_pathways, 1, function(x) {x[3]/x[4]})))] %>% rev()

SGI_results_of_pathways = data.frame(SGI_results_of_pathways, rownames(SGI_results_of_pathways),
                                          SGI_results_of_pathways[,1]/SGI_results_of_pathways[,2],
                                          SGI_results_of_pathways[,3]/SGI_results_of_pathways[,4])
colnames(SGI_results_of_pathways) = c("matched_samples", "total_samples", "matched_genes", "total_pathway_genes",
                                      "Pathway", "Sample_proportion", "Gene_proportion")

SGI_results_of_pathways = SGI_results_of_pathways[signature_order_by_sample,]
SGI_results_of_pathways$Pathway_with_prop = paste0(SGI_results_of_pathways$Pathway, "\n", SGI_results_of_pathways$matched_samples,"/",SGI_results_of_pathways$total_samples)
SGI_results_of_pathways$Pathway_with_prop = factor(SGI_results_of_pathways$Pathway_with_prop, levels = SGI_results_of_pathways$Pathway_with_prop)


######################### plot change
SGI_results_of_pathways = SGI_results_of_pathways[signature_order_by_gene,]
SGI_results_of_pathways$Pathway = factor(SGI_results_of_pathways$Pathway, levels = rev(signature_order_by_gene))
SGI_results_of_pathways$Pathway_with_prop = paste0(SGI_results_of_pathways$Pathway, "\n", SGI_results_of_pathways$matched_genes,"/",SGI_results_of_pathways$total_pathway_genes)
SGI_results_of_pathways$Pathway_with_prop = factor(SGI_results_of_pathways$Pathway_with_prop, levels = SGI_results_of_pathways$Pathway_with_prop)

plot_maf_subset = data.frame(laml@data)
pathway_and_gene_ordering = list()

for(search_signature in gsub("\\.", "-", signature_order_by_sample)){
        tmp = as.character(oncogenic_pathway[oncogenic_pathway$Pathway==search_signature, "Gene"])
        tmp = tmp[tmp %in% plot_maf_subset$Hugo_Symbol]
        tmp_order = lapply(tmp, function(x) {Tumor_Sample_Barcode = plot_maf_subset[plot_maf_subset$Hugo_Symbol==x,"Tumor_Sample_Barcode"] %>% as.character() %>% unique();
                                              length(Tumor_Sample_Barcode)/length(unique(as.character(plot_maf_subset$Tumor_Sample_Barcode)))}) %>% as.numeric()

        tmp = tmp[tmp_order>0.05]
        tmp_order = tmp_order[tmp_order>0.05]
        tmp = tmp[rev(order(tmp_order))]
        pathway_and_gene_ordering[["gene list"]] = c(pathway_and_gene_ordering[["gene list"]], rev(tmp))
        pathway_and_gene_ordering[["annotation"]] = c(pathway_and_gene_ordering[["annotation"]], rep(search_signature, length(tmp)))
}

pathway_and_gene_ordering[["gene list"]] = rev(pathway_and_gene_ordering[["gene list"]])
pathway_and_gene_ordering[["annotation"]] = rev(pathway_and_gene_ordering[["annotation"]])

pathway_and_gene_ordering[["gene list"]] = c(pathway_and_gene_ordering[["gene list"]][seq(1,34)], "MYCN", pathway_and_gene_ordering[["gene list"]][seq(35, length(pathway_and_gene_ordering[["gene list"]]))])
pathway_and_gene_ordering[["annotation"]] = c(pathway_and_gene_ordering[["annotation"]][seq(1,34)], "MYCN", pathway_and_gene_ordering[["annotation"]][seq(35, length(pathway_and_gene_ordering[["annotation"]]))])

pdf(paste0(plot_dir, "Fig3/Fig3A_All_type_oncoplot_by_pathway_geneSet.pdf"), width = 16, height = 10)
color_anno = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1"); names(color_anno) = c("A", "AN", "N", "P", "TN")
oncoplot(maf = laml, sortByAnnotation = TRUE, genes = c(pathway_and_gene_ordering[["gene list"]]), keepGeneOrder = TRUE, clinicalFeatures = c('SCLC_subtype'), annotationColor=list("SCLC_subtype"=color_anno))
dev.off()




###############
##### Fig




damaging_laml = CustomsubsetMaf(laml, query = '(Hugo_Symbol %in% c("TP53", "RB1")) | (Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Amp", "Del"))', mafObj = TRUE, fields = colnames(data.frame(laml@data))[seq(1,99)], clinQuery = 'SCLC_subtype != "TN"')

#enrichment plot condition
#cond 1: A&AN vs N
#cond 2: A&AN vs N

Cond1_maf = CustomsubsetMaf(damaging_laml,  mafObj = TRUE, fields = colnames(data.frame(damaging_laml@data))[seq(1,99)], clinQuery = 'AAN_N_P %in% c("A&AN", "N")')
Cond2_maf = CustomsubsetMaf(damaging_laml,  mafObj = TRUE, fields = colnames(data.frame(damaging_laml@data))[seq(1,99)], clinQuery = 'AAN_N_P %in% c("A&AN", "P")')

Cond1_maf_CE = clinicalEnrichment(maf = Cond1_maf, clinicalFeature = 'AAN_N_P')
pdf(paste0(plot_dir, "Fig3/Fig3C_EnrichmentPlot_condition1_AAN_vs_N_default_function_example.pdf"), width = 7, height = 4)
custom_plotEnrichmentResults(enrich_res = Cond1_maf_CE, pVal = 0.1, geneFontSize = 1, annoFontSize = 1)
dev.off()


Cond2_maf_CE = clinicalEnrichment(maf = Cond2_maf, clinicalFeature = 'AAN_N_P')
pdf(paste0(plot_dir, "Fig3/Fig3B_EnrichmentPlot_condition2_A&AN_vs_P_custom_function_example.pdf"), width = 8, height = 6)
custom_plotEnrichmentResults(enrich_res = Cond2_maf_CE, pVal = 0.1, geneFontSize = 1, annoFontSize = 1)
dev.off()



###########
##### OS
################################ Fig 2B
require(survminer)
require(survival)


thickness = theme(legend.title = element_text(size = 13, color = "black", face = "bold"),
      legend.text = element_text(size = 13, color = "black", face = "bold"),
      axis.text.x = element_text(size = 13, color = "black", face = "bold"),
      axis.text.y = element_text(size = 13, color = "black", face = "bold"),
      axis.title.x = element_text(size = 13, color = "black", face = "bold"),
      axis.title.y = element_text(size = 13, color = "black", face = "bold"))

plot_maf_subset = data.frame(damaging_laml@data)
clinical_maf_subset = data.frame(damaging_laml@clinical.data)
plotting_maf = damaging_laml
print(clinical_maf_subset$SCLC_subtype %>% table())
print(plot_maf_subset$Variant_Classification[!plot_maf_subset$Hugo_Symbol%in%c("TP53", "RB1")] %>% table())
clinical_maf_subset$OS_time = as.numeric(clinical_maf_subset$OS_time)
clinical_maf_subset$OS_event = as.numeric(clinical_maf_subset$OS_event)

maf_subset_table = data.frame(table(clinical_maf_subset$SCLC_subtype))
maf_subset_table$proportion = maf_subset_table$Freq/sum(maf_subset_table$Freq)
rownames(maf_subset_table) = maf_subset_table[,1]
maf_subset_table = maf_subset_table[c("A", "AN", "N", "P"),]
orig_maf_subset_table = maf_subset_table[!is.na(maf_subset_table[,1]),]
surv_results = list()
for(i_gene in c(unique(plot_maf_subset$Hugo_Symbol), "ADGRA2_comb")){
        tryCatch(
                 {variant_exist = rep("WT", nrow(clinical_maf_subset))
		if(i_gene == "ADGRA2_comb"){i_gene = c("ADGRA2", "ZNF703")}
                variant_exist[clinical_maf_subset$Tumor_Sample_Barcode%in%plot_maf_subset$Tumor_Sample_Barcode[plot_maf_subset$Hugo_Symbol%in%i_gene]] = "Mutant"
                clinical_maf_subset$surv_classification = variant_exist
                surv_classification_count = lapply(c("WT", "Mutant"),
                                                        function(x) {
                                                                tmp = data.frame(table(clinical_maf_subset$SCLC_subtype[clinical_maf_subset$surv_classification==x]))
                                                                rownames(tmp) = tmp[,1]
                                                                tmp = tmp[c("A", "AN", "N", "P"),]
                                                                tmp$Var1 = c("A", "AN", "N", "P")
                                                                tmp$Freq[is.na(tmp$Freq)] = 0
                                                                return(tmp)
                                                        }) %>% data.frame()
                surv_classification_count = surv_classification_count[,grepl("Freq", colnames(surv_classification_count))]
                colnames(surv_classification_count) = c("WT", "Mutant")
		rownames(surv_classification_count) = c("A", "AN", "N", "P")
                surv_classification_count = surv_classification_count[!is.na(surv_classification_count[,1]),]
                proportion_plot_df = data.frame(rownames(orig_maf_subset_table), orig_maf_subset_table$proportion, apply(surv_classification_count, 2, function(x) {x/sum(x)}))
                colnames(proportion_plot_df) = c("SCLC_subtype", "subset_proportion", "WT", "Mutant")

		tmp = apply(surv_classification_count[c("A", "AN"),], 2, function(x) {sum(x, na.rm = TRUE)})
		tmp = rbind(tmp, apply(surv_classification_count[c("N", "P"),], 2, function(x) {sum(x, na.rm = TRUE)}))
		surv_classification_count = tmp
		rownames(surv_classification_count) = c("A&AN", "N&P")
		#rownames(surv_classification_count) = c("A&AN", "N", "P")

                surv_classification_proportion = t(data.frame(apply(surv_classification_count, 1, function(x) {x/sum(x)})))

                surv.diff = survival::survdiff(formula = survival::Surv(time = OS_time, event = OS_event) ~ surv_classification, data = clinical_maf_subset)
                surv.diff.pval = signif(1 - pchisq(surv.diff$chisq, length(surv.diff$n) - 1), digits = 3)
                ### Hazard ratio
                surv.cox = survival::coxph(formula = survival::Surv(time = OS_time, event = OS_event) ~ surv_classification, data = clinical_maf_subset)
                hr = signif(1/exp(stats::coef(surv.cox)), digits = 3)

                i_gene = gsub("\\.", "-", i_gene)
                Total_samples = length(unique(as.character(clinical_maf_subset$Tumor_Sample_Barcode)))
                mutant_samples = plot_maf_subset[as.character(plot_maf_subset$Hugo_Symbol) %in% i_gene,"Tumor_Sample_Barcode"] %>% as.character() %>% unique()
                label_group_name = c(paste0(paste0(i_gene,collapse = ";"), " (Mutant: n=",length(mutant_samples),"/",Total_samples,")"),
                                                             paste0(paste0(i_gene, collapse = ";"), " (WT: n=",Total_samples-length(mutant_samples),"/",Total_samples,")"))
                plot_title_info = c(paste0(paste0(rep(" ", 120), collapse = ""),label_group_name[1]),
                paste0(label_group_name[2], " P-value:", round(surv.diff.pval,5), " HR:",round(hr,4), "(",round(fisher.test(surv_classification_count)$p.value,5),")", collapse = " "))
		print_classification = paste0(rownames(surv_classification_proportion)[which.max(surv_classification_proportion[,"Mutant"])])

                if(sum(c("MYC", "MYCL", "CREBBP", "PTEN", "IL7R") %in% i_gene)>0 | identical(i_gene, c("ADGRA2", "ZNF703"))){

                        OS_survfit = survival::survfit(formula = survival::Surv(time = OS_time, event = OS_event) ~ surv_classification, data = clinical_maf_subset)
                        print_majorType = apply(data.frame(plot_maf_subset$Variant_Type[plot_maf_subset$Hugo_Symbol%in%i_gene] %>% table()), 1, function(x) {paste0(x[1], ":", x[2])}) %>% paste0(., collapse = ",")

                        plot_title_info = paste0("Major: ", print_majorType, "\n", label_group_name[1],label_group_name[2])
			p = ggsurvplot(
				  fit = OS_survfit,
				  size = 1.5,
				  xlab = "Months since Treatment",
				  ylab = "OS",
				  surv.median.line = "hv",
				  palette = "nejm",
				  xlim=c(0,30),
				  ylim=c(0,1),
				  break.time.by=6,
				  title=plot_title_info,
				  font.title=c("bold"), font.subtitle=c("italic"),
				  legend.title="Strata",
				  pval = T, pval.coord = c(1, 0.25),
				  surv.scale = "percent",
				  risk.table = TRUE, risk.table.height = 0.3, conf.int = F)
                        p$plot <- p$plot + thickness
                        pdf(paste0(plot_dir, "Fig3/Fig3F_", print_classification, "_enriched_",paste0(i_gene, collapse = ","),"_", length(mutant_samples), "_", Total_samples,".pdf"), width = 5.5, height = 5)
                        print(p)
                        dev.off()

                }
		if(!identical(i_gene, c("ADGRA2", "ZNF703"))){
			surv_results[[paste0(i_gene, collapse = ";")]] = c(surv.diff.pval, hr, fisher.test(surv_classification_count)$p.value, print_classification, length(mutant_samples))}},
		 
        error = function(e) print(clinical_maf_subset$Update_OS_event[clinical_maf_subset$surv_classification=="Mutant"]))
}



Volcano_df = data.frame(t(data.frame(surv_results)))

colnames(Volcano_df) = c("pvalue", "HR", "Fisher", "Enriched", "Mutant_samples")
Volcano_df$HR = log2(as.numeric(Volcano_df$HR))
Volcano_df$pvalue = as.numeric(Volcano_df$pvalue)
Volcano_df$Fisher = as.numeric(Volcano_df$Fisher)
Volcano_df$Mutant_samples  = as.numeric(Volcano_df$Mutant_samples)
Volcano_df$Enriched_group = Volcano_df$Enriched
Volcano_df$Enriched_group[Volcano_df$Enriched_group=="A.AN"] = "Enriched in A or AN"
####################
####################
Volcano_df = Volcano_df[Volcano_df$Mutant_samples>=3,]

####################

Volcano_df$Enriched_group[Volcano_df$Enriched_group=="N.P"] = "Enriched in N or P"

keyvals.colour = rep("red", nrow(Volcano_df))
keyvals.colour[Volcano_df$Enriched==Volcano_df$Enriched[1]] = "blue"
names(keyvals.colour) = Volcano_df$Enriched_group
keyvals.colour[Volcano_df$pvalue>0.1 & abs(Volcano_df$HR)< 3] = "black"
names(keyvals.colour)[Volcano_df$pvalue>0.1 & abs(Volcano_df$HR)< 3] = "non-sig"

sample_count = data.frame(clinical_maf_subset$SCLC_subtype %>% table() %>%
		  reshape2::melt()) %>%
apply(., 1, function(x) {paste0(x[1],":",x[2])}) %>%
paste0(., collapse = "\t")
#p = EnhancedVolcano(Volcano_df, title = paste0(i_type," ED", "Total: ",Total_samples, "\n",sample_count),
p = EnhancedVolcano(Volcano_df, title = paste0("All samples, Prognostic Mutations (OS)"),
    #lab = rownames(Volcano_df), selectLab = rownames(Volcano_df)[Volcano_df$Fisher <= 0.05 & Volcano_df$pvalue <= 0.05],
    lab = rownames(Volcano_df), selectLab = rownames(Volcano_df)[Volcano_df$Mutant_samples>=3 & Volcano_df$pvalue <= 0.1],
    x = 'HR', y = 'pvalue', xlim = c(-3.5, 3.5), xlab = "log2 HR",drawConnectors = TRUE, cutoffLineWidth=0,
    pCutoff = 0.1, colCustom = keyvals.colour, ylim = c(0, 3), colAlpha=1)
ggsave(paste0(plot_dir, "Fig3/Fig3E_VolcanoPlot_SurvDiff_Pvalue_and_HR_Mutant_sample_labeled.png"), plot = p, width = 7, height = 7)


####################################################
#####




damaging_laml = CustomsubsetMaf(laml, query = '(Hugo_Symbol %in% c("TP53", "RB1")) | (Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site"))', mafObj = TRUE, fields = colnames(data.frame(laml@data))[seq(1,99)])

######################################
##################### CNV plotting
######################################

damaging_N_type_maf_subset = CustomsubsetMaf(maf = damaging_laml, clinQuery = 'SCLC_subtype == "N"', mafObj = TRUE, fields = colnames(data.frame(damaging_laml@data))[seq(1,99)])
damaging_A_type_maf_subset = CustomsubsetMaf(maf = damaging_laml, clinQuery = 'SCLC_subtype == "A"', mafObj = TRUE, fields = colnames(data.frame(damaging_laml@data))[seq(1,99)])
damaging_AN_type_maf_subset = CustomsubsetMaf(maf = damaging_laml, clinQuery = 'SCLC_subtype == "AN"', mafObj = TRUE, fields = colnames(data.frame(damaging_laml@data))[seq(1,99)])

CNV_A_type_maf_subset = CustomsubsetMaf(maf = laml, clinQuery = 'SCLC_subtype == "A"', query = '(Hugo_Symbol %in% c("TP53", "RB1")) | (Variant_Classification %in% c("Amp", "Del"))',
                                  mafObj = TRUE, fields = colnames(data.frame(laml@data))[seq(1,99)])
CNV_AN_type_maf_subset = CustomsubsetMaf(maf = laml, clinQuery = 'SCLC_subtype == "AN"', query = '(Hugo_Symbol %in% c("TP53", "RB1")) | (Variant_Classification %in% c("Amp", "Del"))',
                                  mafObj = TRUE, fields = colnames(data.frame(laml@data))[seq(1,99)])
CNV_N_type_maf_subset = CustomsubsetMaf(maf = laml, clinQuery = 'SCLC_subtype == "N"', query = '(Hugo_Symbol %in% c("TP53", "RB1")) | (Variant_Classification %in% c("Amp", "Del"))',
                                  mafObj = TRUE, fields = colnames(data.frame(laml@data))[seq(1,99)])
CNV_P_type_maf_subset = CustomsubsetMaf(maf = laml, clinQuery = 'SCLC_subtype == "P"', query = '(Hugo_Symbol %in% c("TP53", "RB1")) | (Variant_Classification %in% c("Amp", "Del"))',
                                  mafObj = TRUE, fields = colnames(data.frame(laml@data))[seq(1,99)])



SGI_results_of_pathways_results = data.frame()
for(i_type in c("All","damaging_A", "damaging_AN", "damaging_N", "CNV_A", "CNV_AN", "CNV_N")){
        SGI_results_of_pathways = list()
        if(i_type == "All"){
                plot_maf_subset = data.frame(laml@data)
                clinical_maf_subset = data.frame(laml@clinical.data)
        }else{
                plot_maf_subset = data.frame(get(paste0(i_type, "_type_maf_subset"))@data)
                clinical_maf_subset = data.frame(get(paste0(i_type, "_type_maf_subset"))@clinical.data)
        }
        print(table(plot_maf_subset[!(plot_maf_subset$Hugo_Symbol%in%c("TP53", "RB1")), "Variant_Classification"]))

        for(search_signature in unique(oncogenic_pathway$Pathway)){
                search_signature_sample_match = plot_maf_subset$Hugo_Symbol %in% as.character(oncogenic_pathway[oncogenic_pathway$Pathway==search_signature, "Gene"]) %>%
                        plot_maf_subset[., "Tumor_Sample_Barcode"] %>% unique() %>% as.character()
                search_signature_gene_match = as.character(oncogenic_pathway[oncogenic_pathway$Pathway==search_signature, "Gene"]) %in% plot_maf_subset$Hugo_Symbol

                #SGI_results_of_pathways[[search_signature]] = c(matched_samples, total_samples, match_genes, total_pathway_genes)
                Amp = (plot_maf_subset$Hugo_Symbol %in% as.character(oncogenic_pathway[oncogenic_pathway$Pathway==search_signature, "Gene"]) & (plot_maf_subset$Variant_Classification=="Amp")) %>%
                        plot_maf_subset[., "Tumor_Sample_Barcode"] %>% unique() %>% as.character()
                Del = (plot_maf_subset$Hugo_Symbol %in% as.character(oncogenic_pathway[oncogenic_pathway$Pathway==search_signature, "Gene"]) & (plot_maf_subset$Variant_Classification=="Del")) %>%
                        plot_maf_subset[., "Tumor_Sample_Barcode"] %>% unique() %>% as.character()
                SGI_results_of_pathways[[search_signature]] = c(length(search_signature_sample_match), length(unique(clinical_maf_subset$Tumor_Sample_Barcode)),
                                                                sum(search_signature_gene_match), length(search_signature_gene_match), length(Amp), length(Del))
        }
        SGI_results_of_pathways = t(data.frame(SGI_results_of_pathways)); colnames(SGI_results_of_pathways) = c("matched_samples", "total_samples", "matched_genes", "total_pathway_genes", "Amp", "Del")
        if(i_type == "All"){
                signature_order_by_sample = rownames(SGI_results_of_pathways)[rev(order(apply(SGI_results_of_pathways, 1, function(x) {x[1]/x[2]})))] %>% rev()
                signature_order_by_gene = rownames(SGI_results_of_pathways)[rev(order(apply(SGI_results_of_pathways, 1, function(x) {x[3]/x[4]})))] %>% rev()
        }
        SGI_results_of_pathways = data.frame(SGI_results_of_pathways, rownames(SGI_results_of_pathways),
                                                  SGI_results_of_pathways[,1]/SGI_results_of_pathways[,2],
                                                  SGI_results_of_pathways[,3]/SGI_results_of_pathways[,4])
        colnames(SGI_results_of_pathways) = c("matched_samples", "total_samples", "matched_genes", "total_pathway_genes",  "Amp", "Del",
                                              "Pathway", "Sample_proportion", "Gene_proportion")
        SGI_results_of_pathways = SGI_results_of_pathways[signature_order_by_sample,]
        SGI_results_of_pathways$Pathway_with_prop = paste0(SGI_results_of_pathways$Pathway, "\n", SGI_results_of_pathways$matched_samples,"/",SGI_results_of_pathways$total_samples)
        SGI_results_of_pathways$Pathway_with_prop = factor(SGI_results_of_pathways$Pathway_with_prop, levels = SGI_results_of_pathways$Pathway_with_prop)

        SGI_results_of_pathways = SGI_results_of_pathways[signature_order_by_gene,]
        SGI_results_of_pathways$Pathway = factor(SGI_results_of_pathways$Pathway, levels = rev(signature_order_by_gene))
        SGI_results_of_pathways$Pathway_with_prop = paste0(SGI_results_of_pathways$Pathway, "\n", SGI_results_of_pathways$matched_genes,"/",SGI_results_of_pathways$total_pathway_genes)
        SGI_results_of_pathways$Pathway_with_prop = factor(SGI_results_of_pathways$Pathway_with_prop, levels = SGI_results_of_pathways$Pathway_with_prop)
        rownames(SGI_results_of_pathways) = paste0(rownames(SGI_results_of_pathways), "_", i_type)

        if(nrow(SGI_results_of_pathways_results)>0){
                SGI_results_of_pathways_results = rbind(SGI_results_of_pathways_results, SGI_results_of_pathways)
        }else{SGI_results_of_pathways_results = SGI_results_of_pathways}
}


SGI_results_of_pathways_results$SCLC_subtype = lapply(rownames(SGI_results_of_pathways_results), function(x) {x = strsplit(x, "_")[[1]];
       x = x[length(x)]
       if(x=="A"){
               return("SCLC-A")
       }else if(x=="AN"){
               return("SCLC-AN")
       }else if(x=="N"){
               return("SCLC-N")
       }else if(x=="P"){
               return("SCLC-P")
       }else if(x=="TN"){
               return("SCLC-TN")
       }else{return(x)}
                                              }) %>% unlist %>% as.character()

SGI_results_of_pathways_results$SCLC_subtype = factor(SGI_results_of_pathways_results$SCLC_subtype, levels = c("SCLC-A", "SCLC-AN", "SCLC-N", "SCLC-P", "SCLC-TN"))



NOTCH_df = SGI_results_of_pathways_results[grepl("NOTCH_damaging", rownames(SGI_results_of_pathways_results)) & !grepl("All", rownames(SGI_results_of_pathways_results)),]

RTK.RAS_df = SGI_results_of_pathways_results[grepl("RTK.RAS_damaging", rownames(SGI_results_of_pathways_results)) & !grepl("All", rownames(SGI_results_of_pathways_results)),]


PI3K_df = SGI_results_of_pathways_results[grepl("PI3K_damaging", rownames(SGI_results_of_pathways_results)) & !grepl("All", rownames(SGI_results_of_pathways_results)),]

TP53_df = SGI_results_of_pathways_results[grepl("TP53_damaging", rownames(SGI_results_of_pathways_results)) & !grepl("All", rownames(SGI_results_of_pathways_results)),]

CellCycle_df = SGI_results_of_pathways_results[grepl("Cell_Cycle_damaging", rownames(SGI_results_of_pathways_results)) & !grepl("All", rownames(SGI_results_of_pathways_results)),]




Fig3D_plot_df = data.frame(colSums(RTK.RAS_df[c("RTK.RAS_damaging_A", "RTK.RAS_damaging_AN"),c("matched_samples", "total_samples")]),
as.numeric(RTK.RAS_df["RTK.RAS_damaging_N",c("matched_samples", "total_samples")]))
colnames(Fig3D_plot_df) = c("A&AN", "N")
fisher.test(Fig3D_plot_df)
p = ggplot(RTK.RAS_df, aes(x=SCLC_subtype, y = Sample_proportion*100, fill = SCLC_subtype)) + geom_bar(stat="identity") + ylab("% damaging mutations") +
        ggtitle(paste0("RTK-RAS: ", round(fisher.test(Fig3D_plot_df)$p.value,3)))+ xlab("") + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 30))

p = p + theme_bw() + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727")) + axis_theme + theme(legend.position="none")
ggsave(paste0(plot_dir, "Fig3/Fig3D_RTK.RAS_by_samples.png"), plot = p, width = 2.8, height = 4)


Fig3D_plot_df = data.frame(colSums(NOTCH_df[c("NOTCH_damaging_A", "NOTCH_damaging_AN"),c("matched_samples", "total_samples")]),
as.numeric(NOTCH_df["NOTCH_damaging_N",c("matched_samples", "total_samples")]))
colnames(Fig3D_plot_df) = c("A&AN", "N")
fisher.test(Fig3D_plot_df)

p = ggplot(NOTCH_df, aes(x=SCLC_subtype, y = Sample_proportion*100, fill = SCLC_subtype)) + geom_bar(stat="identity") + ylab("% damaging mutations") +
        ggtitle(paste0("NOTCH: ", round(fisher.test(Fig3D_plot_df)$p.value,3))) + xlab("") + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 30))

p = p + theme_bw() + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727")) + axis_theme + theme(legend.position="none")
ggsave(paste0(plot_dir, "Fig3/Fig3D_NOTCH_by_samples.png"), plot = p, width = 2.8, height = 4)


Fig3D_plot_df = data.frame(colSums(PI3K_df[c("PI3K_damaging_A", "PI3K_damaging_AN"),c("matched_samples", "total_samples")]),
as.numeric(PI3K_df["PI3K_damaging_N",c("matched_samples", "total_samples")]))
colnames(Fig3D_plot_df) = c("A&AN", "N")
fisher.test(Fig3D_plot_df)

p = ggplot(PI3K_df, aes(x=SCLC_subtype, y = Sample_proportion*100, fill = SCLC_subtype)) + geom_bar(stat="identity") + ylab("% damaging mutations") +
        ggtitle(paste0("PI3K: ", round(fisher.test(Fig3D_plot_df)$p.value,3))) + xlab("") + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 30))

p = p + theme_bw() + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727")) + axis_theme + theme(legend.position="none")
ggsave(paste0(plot_dir, "Fig3/Fig3D_PI3K_by_samples.png"), plot = p, width = 2.8, height = 4)



Fig3D_plot_df = data.frame(colSums(TP53_df[c("TP53_damaging_A", "TP53_damaging_AN"),c("matched_samples", "total_samples")]),
as.numeric(TP53_df["TP53_damaging_N",c("matched_samples", "total_samples")]))
colnames(Fig3D_plot_df) = c("A&AN", "N")
fisher.test(Fig3D_plot_df)

p = ggplot(TP53_df, aes(x=SCLC_subtype, y = Sample_proportion*100, fill = SCLC_subtype)) + geom_bar(stat="identity") + ylab("% damaging mutations") +
        ggtitle(paste0("TP53_df: ", round(fisher.test(Fig3D_plot_df)$p.value,3))) + xlab("") + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100))

p = p + theme_bw() + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727")) + axis_theme + theme(legend.position="none")
ggsave(paste0(plot_dir, "Fig3/Fig3D_TP53_by_samples.png"), plot = p, width = 2.8, height = 4)


Fig3D_plot_df = data.frame(colSums(CellCycle_df[c("Cell_Cycle_damaging_A", "Cell_Cycle_damaging_AN"),c("matched_samples", "total_samples")]),
as.numeric(CellCycle_df["Cell_Cycle_damaging_N",c("matched_samples", "total_samples")]))
colnames(Fig3D_plot_df) = c("A&AN", "N")
fisher.test(Fig3D_plot_df)

p = ggplot(CellCycle_df, aes(x=SCLC_subtype, y = Sample_proportion*100, fill = SCLC_subtype)) + geom_bar(stat="identity") + ylab("% damaging mutations") +
        ggtitle(paste0("CellCycle_df: ", round(fisher.test(Fig3D_plot_df)$p.value,3))) + xlab("") + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100))

p = p + theme_bw() + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727")) + axis_theme + theme(legend.position="none")
ggsave(paste0(plot_dir, "Fig3/Fig3D_CellCycle_by_samples.png"), plot = p, width = 2.8, height = 4)

################################################

cnv_alteration_by_gene = data.frame(t(data.frame(lapply(c("RICTOR", "SKP2","IL7R", "MYC", "MYCN", "MYCL", "SMAD2", "SMAD4"), function(x) {
               tmp = CNV_A_type_maf_subset@data$Tumor_Sample_Barcode[CNV_A_type_maf_subset@data$Hugo_Symbol==x] %>% unique() %>% length()
               return(c(tmp/length(unique(CNV_A_type_maf_subset@clinical.data$Tumor_Sample_Barcode)), tmp, length(unique(CNV_A_type_maf_subset@clinical.data$Tumor_Sample_Barcode))))
        }))),
t(data.frame(lapply(c("RICTOR", "SKP2","IL7R", "MYC", "MYCN", "MYCL", "SMAD2", "SMAD4"), function(x) {
               tmp = CNV_AN_type_maf_subset@data$Tumor_Sample_Barcode[CNV_AN_type_maf_subset@data$Hugo_Symbol==x] %>% unique() %>% length()
               return(c(tmp/length(unique(CNV_AN_type_maf_subset@clinical.data$Tumor_Sample_Barcode)), tmp, length(unique(CNV_AN_type_maf_subset@clinical.data$Tumor_Sample_Barcode))))
        }))),
t(data.frame(lapply(c("RICTOR", "SKP2","IL7R", "MYC", "MYCN", "MYCL", "SMAD2", "SMAD4"), function(x) {
               tmp = CNV_N_type_maf_subset@data$Tumor_Sample_Barcode[CNV_N_type_maf_subset@data$Hugo_Symbol==x] %>% unique() %>% length()
               return(c(tmp/length(unique(CNV_N_type_maf_subset@clinical.data$Tumor_Sample_Barcode)), tmp, length(unique(CNV_N_type_maf_subset@clinical.data$Tumor_Sample_Barcode))))
        }))),
t(data.frame(lapply(c("RICTOR", "SKP2","IL7R", "MYC", "MYCN", "MYCL", "SMAD2", "SMAD4"), function(x) {
               tmp = CNV_P_type_maf_subset@data$Tumor_Sample_Barcode[CNV_P_type_maf_subset@data$Hugo_Symbol==x] %>% unique() %>% length()
               return(c(tmp/length(unique(CNV_P_type_maf_subset@clinical.data$Tumor_Sample_Barcode)), tmp, length(unique(CNV_P_type_maf_subset@clinical.data$Tumor_Sample_Barcode))))
        }))))

rownames(cnv_alteration_by_gene) = c("RICTOR", "SKP2","IL7R", "MYC", "MYCN", "MYCL", "SMAD2", "SMAD4")
colnames(cnv_alteration_by_gene) = lapply(c("SCLC-A", "SCLC-AN", "SCLC-N", "SCLC-P"), function(x) {paste0(x,"_",c("Percentage", "match", "samples"))}) %>% unlist() %>% as.character()

cnv_alteration_by_gene$gene_name = rownames(cnv_alteration_by_gene)




for(i_gene in c("RICTOR", "SKP2","IL7R")){
        plot_df = cnv_alteration_by_gene[i_gene,grepl("Percentage", colnames(cnv_alteration_by_gene))]# & colnames(cnv_alteration_by_gene)!= "SCLC-P_Percentage"]
        plot_df = reshape::melt(plot_df)
        plot_df$variable = factor(gsub("_Percentage", "", plot_df$variable), levels = c("SCLC-A", "SCLC-AN", "SCLC-N", "SCLC-P"))

        Fisher_test = data.frame(as.numeric(cnv_alteration_by_gene[i_gene, c("SCLC-A_match", "SCLC-AN_match", "SCLC-N_match")]),
                                 as.numeric(cnv_alteration_by_gene[i_gene, c("SCLC-A_samples", "SCLC-AN_samples", "SCLC-N_samples")]))
        Fisher_test = rbind(colSums(Fisher_test[c(1,2),]), Fisher_test[3,])
        colnames(Fisher_test) = c("match", "notMatch")
        rownames(Fisher_test) = c("A&AN", "N")
        p = ggplot(plot_df, aes(x=variable, y = value*100, fill = variable)) + geom_bar(stat="identity") + ylab("% CNV alterations") +
                xlab("") + ggtitle(paste0(i_gene, ":", round(fisher.test(Fisher_test)$p.value,3))) + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 30))
        p = p + theme_bw() + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E")) + axis_theme + theme(legend.position="none")
        ggsave(paste0(plot_dir, "FigS5/FigS5C_CNV_alterated_Samples_", i_gene, ".png"), plot = p, width = 3.5, height = 4)
}


for(i_gene in c("MYC", "MYCL")){
        plot_df = cnv_alteration_by_gene[i_gene,grepl("Percentage", colnames(cnv_alteration_by_gene))]
        plot_df = reshape::melt(plot_df)
        plot_df$variable = factor(gsub("_Percentage", "", plot_df$variable), levels = c("SCLC-A", "SCLC-AN", "SCLC-N", "SCLC-P"))

        Fisher_test = data.frame(as.numeric(cnv_alteration_by_gene[i_gene, c("SCLC-A_match", "SCLC-AN_match", "SCLC-N_match", "SCLC-P_match")]),
                                 as.numeric(cnv_alteration_by_gene[i_gene, c("SCLC-A_samples", "SCLC-AN_samples", "SCLC-N_samples", "SCLC-P_samples")]))
        Fisher_test = rbind(colSums(Fisher_test[c(1,2, 3),]), Fisher_test[4,])
        colnames(Fisher_test) = c("match", "notMatch")
        rownames(Fisher_test) = c("NE", "P")
        p = ggplot(plot_df, aes(x=variable, y = value*100, fill = variable)) + geom_bar(stat="identity") + ylab("% CNV alterations") +
                xlab("") + ggtitle(paste0(i_gene, ":", round(fisher.test(Fisher_test)$p.value,3))) + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 50))
        p = p + theme_bw() + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E")) + axis_theme + theme(legend.position="none")
        ggsave(paste0(plot_dir, "FigS5/FigS5C_CNV_alterated_Samples_", i_gene, ".png"), plot = p, width = 3.5, height = 4)
}




################### Supple Figures




TMB_df = data.frame("SCLC_subtype"=laml@clinical.data$SCLC_subtype, "TMB"=laml@clinical.data$CS_TMB)
TMB_df$TMB = as.numeric(TMB_df$TMB)
plot_df = data.frame(as.character(unique(TMB_df$SCLC_subtype)))
colnames(plot_df) = c("SCLC_subtype")
plot_df$SCLC_subtype = factor(plot_df$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))
plot_df$median = as.numeric(lapply(as.character(unique(TMB_df$SCLC_subtype)),
                 function(x) {median(TMB_df[TMB_df$SCLC_subtype==x,2])}))
plot_df$Q3 = as.numeric(lapply(as.character(unique(TMB_df$SCLC_subtype)),
                 function(x) {quantile(TMB_df[TMB_df$SCLC_subtype==x,2], 0.75)}))
plot_df$Q1 = as.numeric(lapply(as.character(unique(TMB_df$SCLC_subtype)),
                 function(x) {quantile(TMB_df[TMB_df$SCLC_subtype==x,2], 0.25)}))

p = ggplot(plot_df,aes(x=SCLC_subtype, y=median, fill = SCLC_subtype)) + geom_bar(stat="identity") +
        geom_errorbar(mapping = aes(SCLC_subtype,ymin=Q1, ymax=Q3, width = 0.2)) + theme_bw() + guides(fill="none") +
        labs(title=paste0("TMB_SUKSES median by subtype")) + xlab("SCLC subtype") + ylab("Median TMB")
p = p + axis_theme + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E", "#7876B1"))


ggsave(paste0(plot_dir,"FigS5/FigS5A_TMB_SUKSES_median_by_subtype.png"), plot=p, width = 5, height=5)




purity_df = data.frame("SCLC_subtype"=laml@clinical.data$SCLC_subtype, "Purity"=laml@clinical.data$CS_Purity)
purity_df$Purity = as.numeric(purity_df$Purity)

plot_df = data.frame(as.character(unique(matched_meta_info$SCLC_subtype)))
colnames(plot_df) = c("SCLC_subtype")
plot_df$SCLC_subtype = factor(plot_df$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))
plot_df$mean = as.numeric(lapply(as.character(unique(purity_df$SCLC_subtype)),
                 function(x) {mean(purity_df[purity_df$SCLC_subtype==x,2])}))

plot_df$Q3 = as.numeric(lapply(as.character(unique(purity_df$SCLC_subtype)),
                 function(x) {quantile(purity_df[purity_df$SCLC_subtype==x,2], 0.75)}))

plot_df$Q1 = as.numeric(lapply(as.character(unique(purity_df$SCLC_subtype)),
                 function(x) {quantile(purity_df[purity_df$SCLC_subtype==x,2], 0.25)}))



p = ggplot(plot_df,aes(x=SCLC_subtype, y=mean*100, fill = SCLC_subtype)) + geom_bar(stat="identity") +
        geom_errorbar(mapping = aes(SCLC_subtype,ymin=Q1*100, ymax=Q3*100, width = 0.2)) + theme_bw() + guides(fill="none") +
        labs(title=paste0("purity SUKSES mean by subtype")) + xlab("SCLC subtype") + ylab("Mean purity")
p = p + axis_theme + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E", "#7876B1"))
p = p + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 105))
p = p + coord_cartesian(ylim=c(50, 105))

ggsave(paste0(plot_dir, "FigS5/FigS5B_purity_SUKSES_mean_by_subtype.png"), plot=p, width = 5, height=5)







#######################
#####




######################################
##################### damaging with CNV plotting
######################################


damaging_laml = CustomsubsetMaf(laml, query = '(Hugo_Symbol %in% c("TP53", "RB1")) | (Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site" , "Amp", "Del"))', mafObj = TRUE, fields = colnames(data.frame(laml@data))[seq(1,99)])


damaging_N_type_maf_subset = CustomsubsetMaf(maf = damaging_laml, clinQuery = 'SCLC_subtype == "N"', mafObj = TRUE, fields = colnames(data.frame(damaging_laml@data))[seq(1,99)])
damaging_A_type_maf_subset = CustomsubsetMaf(maf = damaging_laml, clinQuery = 'SCLC_subtype == "A"', mafObj = TRUE, fields = colnames(data.frame(damaging_laml@data))[seq(1,99)])
damaging_AN_type_maf_subset = CustomsubsetMaf(maf = damaging_laml, clinQuery = 'SCLC_subtype == "AN"', mafObj = TRUE, fields = colnames(data.frame(damaging_laml@data))[seq(1,99)])
damaging_P_type_maf_subset = CustomsubsetMaf(maf = damaging_laml, clinQuery = 'SCLC_subtype == "P"', mafObj = TRUE, fields = colnames(data.frame(damaging_laml@data))[seq(1,99)])


damaging_alteration_by_gene = data.frame(t(data.frame(lapply(c("MYC", "PTEN", "CDKN2A", "MCL1", "H3F3A", "NOTCH1"), function(x) {
               tmp = damaging_A_type_maf_subset@data$Tumor_Sample_Barcode[damaging_A_type_maf_subset@data$Hugo_Symbol==x] %>% unique() %>% length()
               return(c(tmp/length(unique(damaging_A_type_maf_subset@clinical.data$Tumor_Sample_Barcode)), tmp, length(unique(damaging_A_type_maf_subset@clinical.data$Tumor_Sample_Barcode))))
        }))),
                                    t(data.frame(lapply(c("MYC", "PTEN", "CDKN2A", "MCL1", "H3F3A", "NOTCH1"), function(x) {
               tmp = damaging_AN_type_maf_subset@data$Tumor_Sample_Barcode[damaging_AN_type_maf_subset@data$Hugo_Symbol==x] %>% unique() %>% length()
               return(c(tmp/length(unique(damaging_AN_type_maf_subset@clinical.data$Tumor_Sample_Barcode)), tmp, length(unique(damaging_AN_type_maf_subset@clinical.data$Tumor_Sample_Barcode))))
        }))),
                                    t(data.frame(lapply(c("MYC", "PTEN", "CDKN2A", "MCL1", "H3F3A", "NOTCH1"), function(x) {
               tmp = damaging_N_type_maf_subset@data$Tumor_Sample_Barcode[damaging_N_type_maf_subset@data$Hugo_Symbol==x] %>% unique() %>% length()
               return(c(tmp/length(unique(damaging_N_type_maf_subset@clinical.data$Tumor_Sample_Barcode)), tmp, length(unique(damaging_N_type_maf_subset@clinical.data$Tumor_Sample_Barcode))))
        }))),
t(data.frame(lapply(c("MYC", "PTEN", "CDKN2A", "MCL1", "H3F3A", "NOTCH1"), function(x) {
               tmp = damaging_P_type_maf_subset@data$Tumor_Sample_Barcode[damaging_P_type_maf_subset@data$Hugo_Symbol==x] %>% unique() %>% length()
               return(c(tmp/length(unique(damaging_P_type_maf_subset@clinical.data$Tumor_Sample_Barcode)), tmp, length(unique(damaging_P_type_maf_subset@clinical.data$Tumor_Sample_Barcode))))
        }))))

rownames(damaging_alteration_by_gene) = c("MYC", "PTEN", "CDKN2A", "MCL1", "H3F3A", "NOTCH1")
colnames(damaging_alteration_by_gene) = lapply(c("SCLC-A", "SCLC-AN", "SCLC-N", "SCLC-P"), function(x) {paste0(x,"_",c("Percentage", "match", "samples"))}) %>% unlist() %>% as.character()

damaging_alteration_by_gene$gene_name = rownames(damaging_alteration_by_gene)




for(i_gene in c("MYC", "PTEN", "CDKN2A", "MCL1", "H3F3A", "NOTCH1")){
        plot_df = damaging_alteration_by_gene[i_gene,grepl("Percentage", colnames(damaging_alteration_by_gene))]
        plot_df = reshape::melt(plot_df)
        plot_df$variable = factor(gsub("_Percentage", "", plot_df$variable), levels = c("SCLC-A", "SCLC-AN", "SCLC-N", "SCLC-P"))

        Fisher_test1 = data.frame(as.numeric(damaging_alteration_by_gene[i_gene, c("SCLC-A_match", "SCLC-AN_match")]),
                                 as.numeric(damaging_alteration_by_gene[i_gene, c("SCLC-A_samples", "SCLC-AN_samples")]))
        colnames(Fisher_test1) = c("match", "notMatch")
        rownames(Fisher_test1) = c("A", "AN")

        Fisher_test2 = data.frame(as.numeric(damaging_alteration_by_gene[i_gene, c("SCLC-A_match", "SCLC-AN_match", "SCLC-N_match")]),
                                 as.numeric(damaging_alteration_by_gene[i_gene, c("SCLC-A_samples", "SCLC-AN_samples", "SCLC-N_samples")]))
        Fisher_test2 = rbind(colSums(Fisher_test2[c(1,2),]), Fisher_test2[3,])
        colnames(Fisher_test2) = c("match", "notMatch")
        rownames(Fisher_test2) = c("A&AN", "N")

        Fisher_test3 = data.frame(as.numeric(damaging_alteration_by_gene[i_gene, c("SCLC-A_match", "SCLC-AN_match", "SCLC-P_match")]),
                                 as.numeric(damaging_alteration_by_gene[i_gene, c("SCLC-A_samples", "SCLC-AN_samples", "SCLC-P_samples")]))
        Fisher_test3 = rbind(colSums(Fisher_test3[c(1,2),]), Fisher_test3[3,])
        colnames(Fisher_test3) = c("match", "notMatch")
        rownames(Fisher_test3) = c("A&AN", "p")

        p = ggplot(plot_df, aes(x=variable, y = value*100, fill = variable)) + geom_bar(stat="identity") + ylab("% Damaging & CNV alteration") +
                xlab("") + ggtitle(paste0(i_gene, ":",
                                                 round(fisher.test(Fisher_test1)$p.value,3), ";",
                                                 round(fisher.test(Fisher_test2)$p.value,3), ";",
                                                 round(fisher.test(Fisher_test3)$p.value,3))) + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 20))
        p = p + theme_bw() + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E")) + axis_theme + theme(legend.position="none")
	if(i_gene == "MYC"){
		p = ggplot(plot_df, aes(x=variable, y = value*100, fill = variable)) + geom_bar(stat="identity") + ylab("% Damaging & CNV alteration") +
                xlab("") + ggtitle(paste0(i_gene, ":",
                                                 round(fisher.test(Fisher_test1)$p.value,3), ";",
                                                 round(fisher.test(Fisher_test2)$p.value,3), ";",
						 round(fisher.test(Fisher_test3)$p.value,3))) + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 50))
		p = p + theme_bw() + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E")) + axis_theme + theme(legend.position="none")
	}
        ggsave(paste0(plot_dir, "FigS5/FigS5D_Alterated_Samples_", i_gene, ".png"), plot = p, width = 3.8, height = 4)
}
















