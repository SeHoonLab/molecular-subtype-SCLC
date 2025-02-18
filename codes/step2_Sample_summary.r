LibPath = paste0(getwd(), "/LibPath/")
.libPaths(LibPath)
require(ggplot2)
library(reshape2)
library(ggsci)
library("xlsx")
library(tidyverse)
library(ggpubr)
library(ggsci)
library(cowplot)


plot_dir = "../Figures/"
Signature_df = read.xlsx("../ref_signature/Signature gene list.xlsx", header = TRUE, sheetName="Sheet1")
NE_25_genelist = strsplit(Signature_df[Signature_df$Signature.name=="Neuroendocrine",2], ", ")[[1]]
nonNE_25_genelist = strsplit(Signature_df[Signature_df$Signature.name=="Non-Neuroendocrine",2], ", ")[[1]]
Tuft_cell_marker = strsplit(Signature_df[Signature_df$Signature.name=="Tuft cell marker", 2], ", ")[[1]]



SCLC_meta = read.xlsx("../resource/Essential_check_Edited.xlsx", header = TRUE, sheetName = "Sheet1")
SCLC_subtype_WTS_meta = SCLC_meta[SCLC_meta$WTS_QC_Result %in% c("Pass", "1") & grepl("includ", SCLC_meta$study.inclusion),]
SCLC_TPM =read.table("../resource/log2_TPM_n226.txt", sep = "\t")

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

ggsave(paste0(plot_dir,"Fig1/Fig1C_WTS_subtype.png"), plot = p, width = 5, height = 5)



SCLC_subtype_WTS_meta = SCLC_meta[!is.na(SCLC_meta$SCLC_subtype) & grepl("includ", SCLC_meta$study.inclusion),]
##### Main Fig 1C
ggdat <-  SCLC_subtype_WTS_meta %>% group_by (SCLC_subtype) %>% summarise(count = n()) %>% mutate(prop = round(count/sum(count)*100, digits=1))
rownames(ggdat) = ggdat$SCLC_subtype; ggdat=ggdat[rev(c("A", "AN", "N", "P", "TN")),]
ggdat <- ggdat %>% mutate(lab.ypos = cumsum(prop) - 0.5*prop)
ggdat$SCLC_subtype = factor(ggdat$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))
p = ggplot(ggdat, aes(x = "", y = prop, fill = SCLC_subtype)) +
  geom_bar(width = 2, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = prop), color = "white", size = 3)+
  scale_fill_manual(values = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1")) +
  theme_void()
ggsave(paste0(plot_dir,"Fig1/Fig1C_IHC_subtype.png"), plot = p, width = 5, height = 5)



##################
SCLC_subtype_WTS_meta$SCLC_subtype = factor(SCLC_subtype_WTS_meta$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))
axis_theme = theme(axis.text.x = element_text(size=10, face="bold", color = "black"),
              axis.title.x = element_text(size=10, face="bold", color = "black"),
        axis.text.y = element_text(size=10, face="bold", color = "black"),
        axis.title.y = element_text(size=10, face="bold", color = "black"), panel.border = element_rect(size = 1))


##### Fig S2A, Related to Figure 1


smoking_history = SCLC_subtype_WTS_meta$smoking
smoking_history[!grepl("Never", smoking_history) & !is.na(smoking_history)] = "Smoker"
smoking_history[smoking_history=="Never"] = "Never-\nSmoker"
SCLC_subtype_WTS_meta$smoking2 = factor(smoking_history, levels = c("Smoker", "Never-\nSmoker"))


Fig_S2A_left = ggplot(SCLC_subtype_WTS_meta[!is.na(SCLC_subtype_WTS_meta$smoking2),], aes(x = smoking2, fill = SCLC_subtype)) + geom_bar() + guides(fill="none") +
       	ylab("Patient Number") + xlab("") + theme_bw()
Fig_S2A_left = Fig_S2A_left + ylim(0, 230) + axis_theme + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E", "#7876B1"))


#ggsave(paste0(plot_dir, "FigS2/FigS2D_left_Smoking_distribution.pdf"), plot = Fig_S2A_left, width = 2, height = 3)


smoking_history_df = SCLC_subtype_WTS_meta[!is.na(SCLC_subtype_WTS_meta$smoking2),c("smoking2", "SCLC_subtype")]

smoking_history_df = rbind(lapply(c("Smoker"), function(x) {prop_df = data.frame(table(smoking_history_df$SCLC_subtype[smoking_history_df$smoking2==x]));
	prop_df$Prop = prop_df$Freq/sum(prop_df$Freq); prop_df$subtype2 =rep(x, nrow(prop_df));return(prop_df)}) %>% data.frame(), 
lapply(c("Never-\nSmoker"), function(x) {prop_df = data.frame(table(smoking_history_df$SCLC_subtype[smoking_history_df$smoking2==x]));
        prop_df$Prop = prop_df$Freq/sum(prop_df$Freq); prop_df$subtype2 =rep(x, nrow(prop_df));return(prop_df)}) %>% data.frame())
smoking_history_df$Var1 = factor(smoking_history_df$Var1, levels = c("A", "AN", "N", "P", "TN"))
smoking_history_df$subtype2 = factor(smoking_history_df$subtype2, levels = c("Smoker", "Never-\nSmoker"))

fisher_label = fisher.test(data.frame(smoking_history_df[smoking_history_df$subtype2=="Smoker","Freq"], smoking_history_df[smoking_history_df$subtype2=="Never-\nSmoker","Freq"]))$p.value


Fig_S2A_right = ggplot(smoking_history_df, aes(x = subtype2, y=Prop*100,fill = Var1)) + geom_bar(stat="identity") + guides(fill="none") +
        ylab("Patient proportion(%)") + xlab("") + theme_bw() + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 105))
Fig_S2A_right = Fig_S2A_right  + axis_theme + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E", "#7876B1"))
ggsave(paste0(plot_dir, "FigS2/FigS2D_Smoking_distribution_Prop.pdf"), plot = Fig_S2A_right, width = 2.2, height = 3.5)



##### Fig S2B, Related to Figure 1



histology_arrange = SCLC_subtype_WTS_meta$histology_site
histology_arrange[histology_arrange%in%c("lung", "bronchus")] = "Bronchus & Lung"
histology_arrange[histology_arrange == "brain"] = "Brain"
histology_arrange[!(histology_arrange%in%c("Bronchus & Lung", "LN", "Brain")) & !is.na(histology_arrange)] = "Others"
histology_arrange[histology_arrange=="Bronchus & Lung"] = "Bronchus\n/ Lung"
SCLC_subtype_WTS_meta$histology_site2 = factor(histology_arrange, levels = c("Bronchus\n/ Lung", "LN","Brain", "Others"))


Fig_S2B_left = ggplot(SCLC_subtype_WTS_meta[!is.na(SCLC_subtype_WTS_meta$histology_site2),], aes(x = histology_site2, fill = SCLC_subtype)) + geom_bar()  + ylab("Patient Number") + guides(fill="none") + theme_bw()
Fig_S2B_left = Fig_S2B_left + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E", "#7876B1")) + xlab("") + ylim(0, 230) + axis_theme


#ggsave(paste0(plot_dir, "FigS2/FigS2D_Biopsy_Sites_distribution.pdf"), plot = Fig_S2B_left, width = 3.3, height = 3)




histology_df = SCLC_subtype_WTS_meta[!is.na(SCLC_subtype_WTS_meta$histology_site2),c("histology_site2", "SCLC_subtype")]

histology_df = rbind(lapply(c("LN"), function(x) {prop_df = data.frame(table(histology_df$SCLC_subtype[histology_df$histology_site2==x]));
        prop_df$Prop = prop_df$Freq/sum(prop_df$Freq); prop_df$subtype2 =rep(x, nrow(prop_df));return(prop_df)}) %>% data.frame(),
lapply(c("Bronchus\n/ Lung"), function(x) {prop_df = data.frame(table(histology_df$SCLC_subtype[histology_df$histology_site2==x]));
        prop_df$Prop = prop_df$Freq/sum(prop_df$Freq); prop_df$subtype2 =rep(x, nrow(prop_df));return(prop_df)}) %>% data.frame(),
	    lapply(c("Brain"), function(x) {prop_df = data.frame(table(histology_df$SCLC_subtype[histology_df$histology_site2==x]));
        prop_df$Prop = prop_df$Freq/sum(prop_df$Freq); prop_df$subtype2 =rep(x, nrow(prop_df));return(prop_df)}) %>% data.frame(),
	    lapply(c("Others"), function(x) {prop_df = data.frame(table(histology_df$SCLC_subtype[histology_df$histology_site2==x]));
        prop_df$Prop = prop_df$Freq/sum(prop_df$Freq); prop_df$subtype2 =rep(x, nrow(prop_df));return(prop_df)}) %>% data.frame())

fisher.test(data.frame(histology_df[histology_df$subtype2=="LN","Freq"],
				      histology_df[histology_df$subtype2=="Bronchus\n/ Lung","Freq"]))$p.value %>% print()
fisher.test(data.frame(histology_df[histology_df$subtype2=="Bronchus\n/ Lung","Freq"],
				      histology_df[histology_df$subtype2=="Brain","Freq"]))$p.value %>% print()

histology_df$Var1 = factor(histology_df$Var1, levels = c("A", "AN", "N", "P", "TN"))
histology_df$subtype2 = factor(histology_df$subtype2, levels = c("Bronchus\n/ Lung", "LN", "Brain", "Others"))
Fig_S2B_right = ggplot(histology_df, aes(x = subtype2, y=Prop*100, fill = Var1)) + geom_bar(stat="identity") + guides(fill="none") +
        ylab("Patient proportion(%)") + xlab("") + theme_bw() + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 105))
Fig_S2B_right = Fig_S2B_right + axis_theme + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E", "#7876B1"))
ggsave(paste0(plot_dir, "FigS2/FigS2D_Biopsy_Sites_distribution_Prop.pdf"), plot = Fig_S2B_right, width = 3.5, height = 3)




##### Fig S2C, Related to Figure 1


ini_stage_df = SCLC_subtype_WTS_meta[,c("ini_stage", "SCLC_subtype")]

ini_stage_df = rbind(lapply(c("LD"), function(x) {prop_df = data.frame(table(ini_stage_df$SCLC_subtype[ini_stage_df$ini_stage==x]));
        prop_df$Prop = prop_df$Freq/sum(prop_df$Freq); prop_df$subtype2 =rep(x, nrow(prop_df));return(prop_df)}) %>% data.frame(),
lapply(c("ED"), function(x) {prop_df = data.frame(table(ini_stage_df$SCLC_subtype[ini_stage_df$ini_stage==x]));
        prop_df$Prop = prop_df$Freq/sum(prop_df$Freq); prop_df$subtype2 =rep(x, nrow(prop_df));return(prop_df)}) %>% data.frame())



fisher_label = fisher.test(data.frame(ini_stage_df[ini_stage_df$subtype2=="LD","Freq"],
                                      ini_stage_df[ini_stage_df$subtype2=="ED","Freq"]))$p.value


ini_stage_df$Var1 = factor(ini_stage_df$Var1, levels = c("A", "AN", "N", "P", "TN"))
ini_stage_df$subtype2 = factor(ini_stage_df$subtype2, levels = c("LD", "ED"))
Fig_S2C_right = ggplot(ini_stage_df, aes(x = subtype2, y=Prop*100, fill = Var1)) + geom_bar(stat="identity") + guides(fill="none") +
        ylab("Patient proportion(%)") + xlab("") + theme_bw() + scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 105))
Fig_S2C_right = Fig_S2C_right + axis_theme + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E", "#7876B1"))
ggsave(paste0(plot_dir, "FigS2/FigS2D_LD_ED_IHC_distribution_Prop", round(fisher_label,4),".pdf"), plot = Fig_S2C_right, width = 2.2, height = 3.5)

Fig_S2C_left = ggplot(ini_stage_df, aes(x = subtype2, y=Freq, fill = Var1)) + geom_bar(stat="identity") + guides(fill="none") +
        ylab("Patient Number") + xlab("") + theme_bw() + ylim(0, 200)
Fig_S2C_left = Fig_S2C_left + axis_theme + scale_fill_manual(values=c("#BC3C29", "#0072B5","#E18727", "#20854E", "#7876B1"))
#ggsave(paste0(plot_dir, "FigS2/FigS2D_LD_ED_IHC_distribution.pdf"), plot = Fig_S2C_left, width = 1.8, height = 3)




######## Fig2


clinic.info <- SCLC_meta
#View(clinic.info)
clinic.info.OS <- clinic.info %>% filter(!is.na(SCLC_subtype)) %>% filter(grepl("included", study.inclusion))
nrow(clinic.info.OS) #n=252 for OS analysis
table(clinic.info.OS$SCLC_subtype)
clinic.info.OS$NE <- "NE(A/N/AN)"
clinic.info.OS$NE[clinic.info.OS$SCLC_subtype %in% c("P","TN")] <- "Non-NE(P/TN)"
table(clinic.info.OS$NE)


clinic.info.OS$IO <- as.factor(as.character(clinic.info.OS$IO))
table(clinic.info.OS$IO)

clinic.info.OS$IO_line[is.na(clinic.info.OS$IO_line)] <- "0"
clinic.info.OS$IO_line <- as.factor(clinic.info.OS$IO_line)
table(clinic.info.OS$IO_line)

clinic.info.OS$first_IO <- "first-line ICB"
clinic.info.OS$first_IO[clinic.info.OS$IO == 0] <- "No ICB"
clinic.info.OS$first_IO[clinic.info.OS$IO == 1 &
                          clinic.info.OS$IO_line %in% c(2, 3, 4, 5, 6, 7)] <- "later-line ICB"
clinic.info.OS$first_IO <- as.factor(clinic.info.OS$first_IO)
table(clinic.info.OS$first_IO)


clinic.info.OS$OS_time
clinic.info.OS$ctx_1st_pfs <- as.numeric(clinic.info.OS$ctx_1st_pfs)

clinic.info.OS$SCLC_subtype = factor(clinic.info.OS$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))

clinic.info.OS.LD <- clinic.info.OS %>% filter(grepl("LD", clinic.info.OS$ini_stage))
clinic.info.OS.ED <- clinic.info.OS %>% filter(!grepl("LD", clinic.info.OS$ini_stage))
clinic.info.OS.ED.IO <- clinic.info.OS.ED %>% filter(IO_line == 1)
clinic.info.OS.ED.CTx <- clinic.info.OS.ED %>% filter(IO_line != 1)

clinic.info.OS.ED.NE <- clinic.info.OS.ED %>% filter(NE == "NE(A/N/AN)")

clinic.info.ED.PFS <- clinic.info.OS.ED %>% filter(!is.na(ctx_1st_pfs_event))
clinic.info.ED.PFS.IO <- clinic.info.ED.PFS %>% filter(IO_line == 1)
clinic.info.ED.PFS.CTx <- clinic.info.ED.PFS %>% filter(IO_line != 1)


################################ Fig 2B
require(survminer)
require(survival)


thickness = theme(legend.title = element_text(size = 13, color = "black", face = "bold"),
      legend.text = element_text(size = 13, color = "black", face = "bold"),
      axis.text.x = element_text(size = 13, color = "black", face = "bold"),
      axis.text.y = element_text(size = 13, color = "black", face = "bold"),
      axis.title.x = element_text(size = 13, color = "black", face = "bold"),
      axis.title.y = element_text(size = 13, color = "black", face = "bold"))


##### Fig 2A
clinic.info.OS.LD$SCLC_subtype = factor(clinic.info.OS.LD$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))

OS_fit <- survfit(Surv(OS_time, OS_event==1)~SCLC_subtype, data = clinic.info.OS.LD)
cox_fit = coxph(Surv(OS_time, OS_event==1)~SCLC_subtype, data = clinic.info.OS.LD) # 220628 OS -> OS_time

pdf(paste0(plot_dir, "Fig2/Fig2A_OS_Analysis_IHC_subtype_LD.pdf"), width = 8, height = 6.2)
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
  title="OS Analysis in LD-patients: IHC subtype",
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
colnames(overall_summary) = c("group", "mPFS", "HR", "P")
overall_summary$group = gsub(" ", "", overall_summary$group)
overall_summary$P = round(overall_summary$P, 3)

apply(overall_summary, 1, function(x) {paste0(x[seq(1,4)], collapse = ";")})


##### Fig 2B

clinic.info.OS.ED$SCLC_subtype = factor(clinic.info.OS.ED$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))


OS_fit <- survfit(Surv(OS_time, OS_event==1)~SCLC_subtype, data = clinic.info.OS.ED)
cox_fit = coxph(Surv(OS_time, OS_event==1)~SCLC_subtype, data = clinic.info.OS.ED)

pdf(paste0(plot_dir, "Fig2/Fig2B_OS_Analysis_IHC_subtype_ED.pdf"), width = 8, height = 6.2)
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
  title="OS Analysis in ED-patients: IHC subtype",
  font.title=c("bold"), font.subtitle=c("italic"), font.size=13,
  legend.title="Strata",
  pval = F, pval.coord = c(1, 0.25),
  surv.scale = "percent",
  risk.table = TRUE, risk.table.height = 0.3, conf.int = F, newpage = FALSE)
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



##### Fig 2C

OS_fit <- survfit(Surv(OS_time, OS_event==1)~SCLC_subtype, data = clinic.info.OS.ED.CTx)
cox_fit = coxph(Surv(OS_time, OS_event==1)~SCLC_subtype, data = clinic.info.OS.ED.CTx) # 220628 OS -> OS_time

pdf(paste0(plot_dir, "Fig2/Fig2C_OS_Analysis_IHC_subtype_ED_CTX.pdf"), width = 8, height = 6.2)
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
  title="OS Analysis in ED-patients: CTx only",
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
colnames(overall_summary) = c("group", "mPFS", "HR", "P")
overall_summary$group = gsub(" ", "", overall_summary$group)
overall_summary$P = round(overall_summary$P, 3)

apply(overall_summary, 1, function(x) {paste0(x[seq(1,4)], collapse = ";")})



##### Fig 2D


OS_fit <- survfit(Surv(OS_time, OS_event==1)~SCLC_subtype, data = clinic.info.OS.ED.IO)
cox_fit = coxph(Surv(OS_time, OS_event==1)~SCLC_subtype, data = clinic.info.OS.ED.IO) 

pdf(paste0(plot_dir, "Fig2/Fig2D_OS_Analysis_IHC_subtype_ED_IO.pdf"), width = 8, height = 6)
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
  title="OS Analysis in ED-patients: IO",
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
colnames(overall_summary) = c("group", "mPFS", "HR", "P")
overall_summary$group = gsub(" ", "", overall_summary$group)
overall_summary$P = round(overall_summary$P, 3)

apply(overall_summary, 1, function(x) {paste0(x[seq(1,4)], collapse = ";")})




########## Figure S3
##### Fig S3A

clinic.info.ED.PFS$ctx_1st_pfs = as.numeric(clinic.info.ED.PFS$ctx_1st_pfs)
clinic.info.ED.PFS$ctx_1st_pfs_event = as.numeric(clinic.info.ED.PFS$ctx_1st_pfs_event)


Pfs_fit <- survfit(Surv(ctx_1st_pfs, ctx_1st_pfs_event==1)~SCLC_subtype, data = clinic.info.ED.PFS)
cox_fit = coxph(Surv(ctx_1st_pfs, ctx_1st_pfs_event==1)~SCLC_subtype, data = clinic.info.ED.PFS)

pdf(paste0(plot_dir, "FigS4/FigS4A_ctx_1st_pfs_IHC_subtype_ED.pdf"), width = 8, height = 6.2)
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
  title="Pfs Analysis in ED-patients: IHC subtype",
  font.title=c("bold"), font.subtitle=c("italic"), font.size=13,
  legend.title="Strata",
  pval = F, pval.coord = c(1, 0.25),
  surv.scale = "percent",
  risk.table = TRUE, risk.table.height = 0.3, conf.int = F, newpage = FALSE)
p$plot <- p$plot + thickness
print(p)
dev.off()



Pfs_fit_summary = data.frame(summary(Pfs_fit)$table)
rownames(Pfs_fit_summary) = gsub("=", "", rownames(Pfs_fit_summary))
Pfs_fit_summary$rownames = rownames(Pfs_fit_summary)


cox_summary = cbind(data.frame(summary(cox_fit)$coefficients), data.frame(summary(cox_fit)$conf.int))
cox_summary$rownames = rownames(cox_summary)


overall_summary = merge(Pfs_fit_summary, cox_summary, by="rownames", all=TRUE)

overall_summary = (data.frame(overall_summary$rownames,
                 paste0(round(overall_summary$median, 1), " (", round(overall_summary$X0.95LCL, 1), "-", round(overall_summary$X0.95UCL, 1), ")"),
                 paste0(round(overall_summary$exp.coef., 1), " (", round(overall_summary$lower..95, 1), "-", round(overall_summary$upper..95, 1),")"),
                 overall_summary$Pr...z..))
colnames(overall_summary) = c("group", "mPFS", "HR", "P")
overall_summary$group = gsub(" ", "", overall_summary$group)
overall_summary$P = round(overall_summary$P, 3)

apply(overall_summary, 1, function(x) {paste0(x[seq(1,4)], collapse = ";")})


##### Fig S3B


Pfs_fit <- survfit(Surv(ctx_1st_pfs, ctx_1st_pfs_event==1)~SCLC_subtype, data = clinic.info.ED.PFS.CTx)
cox_fit = coxph(Surv(ctx_1st_pfs, ctx_1st_pfs_event==1)~SCLC_subtype, data = clinic.info.ED.PFS.CTx)


pdf(paste0(plot_dir, "FigS4/FigS4B_PFS_IHC_subtype_ED_CTx_only.pdf"), width = 8, height = 6.2)
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
  title="Pfs Analysis in ED-patients: CTx only",
  font.title=c("bold"), font.subtitle=c("italic"), font.size=13,
  legend.title="Strata",
  pval = F, pval.coord = c(1, 0.25),
  surv.scale = "percent",
  risk.table = TRUE, risk.table.height = 0.3, conf.int = F, newpage = FALSE)
p$plot <- p$plot + thickness
print(p)
dev.off()



Pfs_fit_summary = data.frame(summary(Pfs_fit)$table)
rownames(Pfs_fit_summary) = gsub("=", "", rownames(Pfs_fit_summary))
Pfs_fit_summary$rownames = rownames(Pfs_fit_summary)


cox_summary = cbind(data.frame(summary(cox_fit)$coefficients), data.frame(summary(cox_fit)$conf.int))
cox_summary$rownames = rownames(cox_summary)


overall_summary = merge(Pfs_fit_summary, cox_summary, by="rownames", all=TRUE)

overall_summary = (data.frame(overall_summary$rownames,
                 paste0(round(overall_summary$median, 1), " (", round(overall_summary$X0.95LCL, 1), "-", round(overall_summary$X0.95UCL, 1), ")"),
                 paste0(round(overall_summary$exp.coef., 1), " (", round(overall_summary$lower..95, 1), "-", round(overall_summary$upper..95, 1),")"),
                 overall_summary$Pr...z..))
colnames(overall_summary) = c("group", "mPFS", "HR", "P")
overall_summary$group = gsub(" ", "", overall_summary$group)
overall_summary$P = round(overall_summary$P, 3)

apply(overall_summary, 1, function(x) {paste0(x[seq(1,4)], collapse = ";")})


##### Fig S3C


Pfs_fit <- survfit(Surv(ctx_1st_pfs, ctx_1st_pfs_event==1)~SCLC_subtype, data = clinic.info.ED.PFS.IO)
cox_fit = coxph(Surv(ctx_1st_pfs, ctx_1st_pfs_event==1)~SCLC_subtype, data = clinic.info.ED.PFS.IO)


pdf(paste0(plot_dir, "FigS4/FigS4C_PFS_IHC_subtype_ED_CTx_IO.pdf"), width = 8, height = 6)
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
  title="Pfs Analysis in ED-patients: IO + CTX",
  font.title=c("bold"), font.subtitle=c("italic"), font.size=13,
  legend.title="Strata",
  pval = F, pval.coord = c(1, 0.25),
  surv.scale = "percent",
  risk.table = TRUE, risk.table.height = 0.3, conf.int = F, newpage = FALSE)
p$plot <- p$plot + thickness
print(p)
dev.off()




Pfs_fit_summary = data.frame(summary(Pfs_fit)$table)
rownames(Pfs_fit_summary) = gsub("=", "", rownames(Pfs_fit_summary))
Pfs_fit_summary$rownames = rownames(Pfs_fit_summary)


cox_summary = cbind(data.frame(summary(cox_fit)$coefficients), data.frame(summary(cox_fit)$conf.int))
cox_summary$rownames = rownames(cox_summary)


overall_summary = merge(Pfs_fit_summary, cox_summary, by="rownames", all=TRUE)

overall_summary = (data.frame(overall_summary$rownames,
                 paste0(round(overall_summary$median, 1), " (", round(overall_summary$X0.95LCL, 1), "-", round(overall_summary$X0.95UCL, 1), ")"),
                 paste0(round(overall_summary$exp.coef., 1), " (", round(overall_summary$lower..95, 1), "-", round(overall_summary$upper..95, 1),")"),
                 overall_summary$Pr...z..))
colnames(overall_summary) = c("group", "mPFS", "HR", "P")
overall_summary$group = gsub(" ", "", overall_summary$group)
overall_summary$P = round(overall_summary$P, 3)

apply(overall_summary, 1, function(x) {paste0(x[seq(1,4)], collapse = ";")})



#####################
##### Fig S1C
library(ggpubr)

SCLC_subtype_WTS_meta = SCLC_meta[!is.na(SCLC_meta$SCLC_subtype) &  grepl("includ", SCLC_meta$study.inclusion),]
SCLC_subtype_WTS_meta$SCLC_subtype = factor(SCLC_subtype_WTS_meta$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))

SCLC_subtype_WTS_meta$NE_subtype = lapply(SCLC_subtype_WTS_meta$SCLC_subtype, function(x) {
                                                           if(x=="A" | x=="N" | x=="AN"){
                                                                return("NE")
                                                           }else if(x=="P" | x=="TN"){
                                                                return("Non-NE")
                                                           }
}) %>% as.character()
SCLC_subtype_WTS_meta$NE_subtype = factor(SCLC_subtype_WTS_meta$NE_subtype, levels = c("NE", "Non-NE"))

SCLC_subtype_WTS_meta = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$SCLC_subtype%in%c("A", "AN", "N"),]
SCLC_subtype_WTS_meta$ASCL1_H_Score = as.numeric(SCLC_subtype_WTS_meta$ASCL1_H_Score)
SCLC_subtype_WTS_meta$NEUROD1_H_Score = as.numeric(SCLC_subtype_WTS_meta$NEUROD1_H_Score)



p = ggscatter(SCLC_subtype_WTS_meta, x="NEUROD1_H_Score", y = "ASCL1_H_Score", shape = "SCLC_subtype", color = "SCLC_subtype", palette = c("#BC3C29", "#0072B5", "#E18727"), ellipse=TRUE)
ggsave(paste0(plot_dir, "FigS1/FigS1B_ASCL1_NEUROD1_H_score_scatterPlot.pdf"), plot = p, width = 7, height = 7)





###############
##### Fig S1F

SCLC_subtype_WTS_meta = SCLC_meta[SCLC_meta$WTS_QC_Result %in% c("Pass", "1") & grepl("includ", SCLC_meta$study.inclusion) & !is.na(SCLC_meta$SCLC_subtype),]
SCLC_subtype_WTS_meta$NE_subtype = lapply(SCLC_subtype_WTS_meta$SCLC_subtype, function(x) {
                                                           if(x=="A" | x=="N" | x=="AN"){
                                                                return("NE")
                                                           }else if(x=="P" | x=="TN"){
                                                                return("Non-NE")
                                                           }
}) %>% as.character()
SCLC_subtype_WTS_meta$SCLC_subtype = factor(SCLC_subtype_WTS_meta$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))


library(singscore)

SCLC_TPM2 = SCLC_TPM[,SCLC_subtype_WTS_meta$WTS_ID]


SCLC_subtype_WTS_meta$ASCL1_expr = as.numeric(SCLC_TPM2["ASCL1", ])
SCLC_subtype_WTS_meta$NEUROD1_expr = as.numeric(SCLC_TPM2["NEUROD1", ])
SCLC_subtype_WTS_meta$POU2F3_expr = as.numeric(SCLC_TPM2["POU2F3", ])
SCLC_subtype_WTS_meta$YAP1_expr = as.numeric(SCLC_TPM2["YAP1", ])






SCLC_TPM_Rank <- rankGenes(SCLC_TPM2)

P_vs_NE = c()
N_vs_AAN = c()
results_rownames = c()

my_comparisons <- list(c("A", "AN"), c("AN", "N"))
axis_theme = theme(axis.text.x = element_text(size=10, face="bold", color = "black"),
              axis.title.x = element_text(size=10, face="bold", color = "black"),
        axis.text.y = element_text(size=10, face="bold", color = "black"),
        axis.title.y = element_text(size=10, face="bold", color = "black"), panel.border = element_rect(size = 1))

for(i_name in c("ASCL1 & ND1 Shared targets", "ASCL1 high signatures cell line",
                "ND1 high signatures cell line")){
        gene_list = strsplit(Signature_df[Signature_df$Signature.name==i_name, 2], ", ")[[1]] %>% unlist() %>% as.character() %>% unique()
            SingScore <- simpleScore(
          rankData = SCLC_TPM_Rank,
          upSet = gene_list,
          centerScore = T,
          knownDirection = T
          )
	if(i_name == "ASCL1 & ND1 Shared targets"){i_name="ASCL1_ND1_Shared_targets"}else if(i_name=="ASCL1 high signatures cell line"){
                i_name = "ASCL1_high_signatures_cell_line"}else if(i_name == "ND1 high signatures cell line"){i_name="ND1_high_signatures_cell_line"}
            if(nrow(SingScore)>0){
                    identical(rownames(SingScore), SCLC_subtype_WTS_meta$WTS_ID) %>% print()

                    plot_df = data.frame(SingScore$TotalScore, SCLC_subtype_WTS_meta$SCLC_subtype, SCLC_subtype_WTS_meta$NE_subtype)
                    colnames(plot_df) = c("Singscore", "SCLC_subtype", "NE_subtype")
                    plot_df$SCLC_subtype = factor(plot_df$SCLC_subtype, levels = c("A", "AN", "N", "P", "TN"))
                    pvalue = c(wilcox.test(plot_df$Singscore[plot_df$NE_subtype=="NE"], plot_df$Singscore[plot_df$SCLC_subtype=="P"])$p.value,
                               wilcox.test(plot_df$Singscore[plot_df$SCLC_subtype=="N"], plot_df$Singscore[plot_df$SCLC_subtype=="A" | plot_df$SCLC_subtype=="AN"])$p.value)
                    rawFC = c(mean(plot_df$Singscore[plot_df$SCLC_subtype=="P"])>mean(plot_df$Singscore[plot_df$NE_subtype=="NE"]),
                              mean(plot_df$Singscore[plot_df$SCLC_subtype=="N"])>mean(plot_df$Singscore[plot_df$SCLC_subtype=="A" | plot_df$SCLC_subtype=="AN"]))
                    if(pvalue[1] <= 0.05 & rawFC[1]==TRUE){
                            P_vs_NE = c(P_vs_NE, "P Up")
                    }else if(pvalue[1] <= 0.05 & rawFC[1]==FALSE){
                            P_vs_NE = c(P_vs_NE, "P Down")
                    }else if(round(pvalue[1], 1)<=0.1 & rawFC[1]==TRUE){
                            P_vs_NE = c(P_vs_NE, "P Up(lessSig)")
                    }else if(round(pvalue[1], 1)<=0.1 & rawFC[1]==FALSE){
                            P_vs_NE = c(P_vs_NE, "P Down(lessSig)")
                    }else{P_vs_NE = c(P_vs_NE, "")}

                    if(pvalue[2] <= 0.05 & rawFC[2]==TRUE){
                            N_vs_AAN = c(N_vs_AAN, "N Up")
                    }else if(pvalue[2] <= 0.05 & rawFC[2]==FALSE){
                            N_vs_AAN = c(N_vs_AAN, "N Down")
                    }else if(round(pvalue[2], 1)<=0.1 & rawFC[2]==TRUE){
                            N_vs_AAN = c(N_vs_AAN, "N Up(lessSig)")
                    }else if(round(pvalue[2], 1)<=0.1 & rawFC[2]==FALSE){
                            N_vs_AAN = c(N_vs_AAN, "N Down(lessSig)")
                    }else{N_vs_AAN = c(N_vs_AAN, "")}

                    results_rownames = c(results_rownames, i_name)
                    if(sum(pvalue<0.05)>=1){plot_dir_loop = paste0(plot_dir, "Sig/")}else if(sum(round(pvalue,1)<=0.1)>=1){plot_dir_loop = paste0(plot_dir, "lessSig/")}else{plot_dir_loop= plot_dir}
                    p = ggplot(plot_df, aes(x=SCLC_subtype, y = Singscore)) + geom_boxplot(aes(fill = SCLC_subtype)) + scale_fill_nejm() + theme_bw()
                    p = p + ggtitle(paste0(i_name, "\n", "wilcox(NE vs P):",round(pvalue[1], 5), "\n","wilcox(N vs A&AN):", round(pvalue[2], 5)))
                    p = p + stat_compare_means(comparisons = my_comparisons, method = "wilcox") + axis_theme + theme(legend.position="none") + ylim(min(plot_df$Singscore), max(plot_df$Singscore)+0.1)
                    ggsave(paste0(plot_dir, "FigS2/FigS2A_", i_name, ".pdf"), plot = p, width = 4, height = 4)
                    p = p + ggtitle("")
                    #ggsave(paste0(plot_dir, "FigS2/FigS2A", i_name, ".png"), plot = p, width = 4, height = 4)
            }
  }




SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = NE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$NE_signature = SingScore$TotalScore

SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = nonNE_25_genelist,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$`non-NE_signature` = SingScore$TotalScore


for(i_poster_gene in c("NE_signature", "non-NE_signature")){
        tmp = SCLC_subtype_WTS_meta[,c("SCLC_subtype", i_poster_gene)]
        colnames(tmp) = c("SCLC_subtype", "gene_expr")
        pvalue = wilcox.test(tmp[tmp$SCLC_subtype%in%c("A", "AN", "N"),"gene_expr"], tmp[tmp$SCLC_subtype%in%c("P"),"gene_expr"])$p.value %>% round(.,4)
        pvalue = wilcox.test(tmp[tmp$SCLC_subtype%in%c("A", "AN"),"gene_expr"], tmp[tmp$SCLC_subtype%in%c("N"),"gene_expr"])$p.value %>% round(.,4)

        print(pvalue)
        p = ggplot(tmp, aes(x=SCLC_subtype, y=gene_expr, fill = SCLC_subtype)) + geom_boxplot() #ggtitle(paste0(i_poster_gene, ": 182 WTS\nwilcox(A&AN&N vs P):", pvalue)) + ylab(i_poster_gene)
        p = p + theme_bw() + scale_fill_nejm() + theme(legend.position="none") +  ylab(gsub("_", " ", i_poster_gene)) + xlab("")
        ggsave(paste0(plot_dir, "FigS2/FigS2A_", i_poster_gene, "_CancerCell_gene.png"), plot = p + axis_theme, width = 4, height = 4)
}


##### Tuft cell marker


SingScore <- simpleScore(
  rankData = SCLC_TPM_Rank,
  upSet = Tuft_cell_marker,
  centerScore = T,
  knownDirection = T
  )

SCLC_subtype_WTS_meta$Tuft_cell_marker = SingScore$TotalScore

tmp = SCLC_subtype_WTS_meta[,c("SCLC_subtype", "Tuft_cell_marker")]
colnames(tmp) = c("SCLC_subtype", "gene_expr")
pvalue = wilcox.test(tmp[tmp$SCLC_subtype%in%c("A", "AN", "N"),"gene_expr"], tmp[tmp$SCLC_subtype%in%c("P"),"gene_expr"])$p.value %>% round(.,7)
p = ggplot(tmp, aes(x=SCLC_subtype, y=gene_expr, fill = SCLC_subtype)) + geom_boxplot() #ggtitle(paste0(i_poster_gene, ": 182 WTS\nwilcox(A&AN&N vs P):", pvalue)) + ylab(i_poster_gene)
p = p + theme_bw() + scale_fill_nejm() + theme(legend.position="none") +  ylab("Tuft cell signature") + xlab("") #+ xlab("IHC subtype") + ggtitle(paste0(paste0(c("Singscore: Tuft cell markers"), collapse = ";"), "\nwilcox (A&AN&N vs P)", pvalue))

ggsave(paste0(plot_dir, "FigS2/FigS2A_Tuft_gene.png"), plot = p + axis_theme, width = 4, height = 4)





######### NMF


SCLC_subtype_WTS_meta$NMF_subtype = factor(SCLC_subtype_WTS_meta$NMF_subtype, levels = c("A", "N", "P", "I", "equivocal"))



OS_fit <- survfit(Surv(OS_time, OS_event==1)~NMF_subtype, data = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$ini_stage=="ED" & SCLC_subtype_WTS_meta$NMF!="equivocal",])
cox_fit = coxph(Surv(OS_time, OS_event==1)~NMF_subtype, data = SCLC_subtype_WTS_meta[SCLC_subtype_WTS_meta$ini_stage=="ED" & SCLC_subtype_WTS_meta$NMF!="equivocal",])

pdf(paste0(plot_dir, "FigS4/FigS4F_OS_Analysis_NMF_subtype_ED.pdf"), width = 8, height = 6.2)
p = ggsurvplot(
  fit = OS_fit,
  size = 1.5,
  xlab = "Months since Treatment",
  ylab = "OS",
  surv.median.line = "hv",
  palette = c("#BC3C29", "#E18727", "#20854E", "#7876B1"),
  xlim=c(0,30),
  ylim=c(0,1),
  break.time.by=6,
  title="OS Analysis in ED-patients: NMF subtype",
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





#################################



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

Confusion_matrix_df = data.frame(SCLC_subtype_WTS_meta$NMF_subtype, SCLC_subtype_WTS_meta$SCLC_subtype)
Confusion_matrix_df = Confusion_matrix_df[Confusion_matrix_df[,1]!="equivocal",]
colnames(Confusion_matrix_df) = c("NMF_subtype", "IHC_subtype")

##### Calculate kappa
kappa_df = Confusion_matrix_df
kappa_df$NMF_subtype = as.character(kappa_df$NMF_subtype)
kappa_df$NMF_subtype[kappa_df$NMF_subtype!="A"] = "non-A"
kappa_df$IHC_subtype = as.character(kappa_df$IHC_subtype)
kappa_df$IHC_subtype[kappa_df$IHC_subtype!="A"] = "non-A"
cfm3 <- confusionMatrix(factor(kappa_df$IHC_subtype, levels = c("A", "non-A")),
                        factor(kappa_df$NMF_subtype, levels = c("A", "non-A")))
print("FigS3B : A subtype confusion matrix\n")
print(cfm3)



##### Calculate kappa
kappa_df = Confusion_matrix_df
kappa_df$NMF_subtype = as.character(kappa_df$NMF_subtype)
kappa_df$NMF_subtype[kappa_df$NMF_subtype!="N"] = "non-N"
kappa_df$IHC_subtype = as.character(kappa_df$IHC_subtype)
kappa_df$IHC_subtype[kappa_df$IHC_subtype!="N"] = "non-N"
cfm3 <- confusionMatrix(factor(kappa_df$IHC_subtype, levels = c("N", "non-N")),
                        factor(kappa_df$NMF_subtype, levels = c("N", "non-N")))
print("FigS3B : N subtype confusion matrix\n")
print(cfm3)



##### Calculate kappa
kappa_df = Confusion_matrix_df
kappa_df$NMF_subtype = as.character(kappa_df$NMF_subtype)
kappa_df$NMF_subtype[kappa_df$NMF_subtype!="P"] = "non-P"
kappa_df$IHC_subtype = as.character(kappa_df$IHC_subtype)
kappa_df$IHC_subtype[kappa_df$IHC_subtype!="P"] = "non-P"
cfm3 <- confusionMatrix(factor(kappa_df$IHC_subtype, levels = c("P", "non-P")),
                        factor(kappa_df$NMF_subtype, levels = c("P", "non-P")))
print("FigS3B : P subtype confusion matrix\n")
print(cfm3)


