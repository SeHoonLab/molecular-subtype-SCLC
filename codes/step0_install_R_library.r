LibPath = paste0(getwd(), "/LibPath/")
system(paste0("mkdir ", LibPath))
.libPaths(LibPath)
system("mkdir ../Figures/")
install.packages("pacman", repos = "https://cloud.r-project.org/")
library(pacman)
install.packages('BiocManager', repos = "https://cloud.r-project.org/")
install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", repos=NULL, type="source")
install.packages("devtools", repos = "https://cloud.r-project.org/")

devtools::install_github('kevinblighe/EnhancedVolcano')

p_install_version(
    c("edgeR", "reshape2", "ggsci", "xlsx", "tidyverse", "cowplot", "fgsea", "scales", "singscore", "ggpubr", "circlize", "ComplexHeatmap", "survminer", "survival", "maftools", "data.table", "dplyr", "magrittr", "Hmisc", "NMF", "sva", "RColorBrewer", "GSVA", "tweeDEseq", "pracma", "ggplot2", "caret", "org.Hs.eg.db", "tibble", "stringr"),
    c("3.32.1", "1.4.4", "2.9", "0.6.5", "1.3.1", "1.1.1", "1.16.0", "1.1.1", "1.10.0", "0.4.0", "0.4.13", "2.6.2", "0.4.9", "3.2.7", "2.6.5", "1.14.0", "1.0.8", "2.0.2", "4.5.0", "0.23.0", "3.38.0", "1.1.2", "1.38.2", "1.36.0", "2.3.3", "3.3.5", "6.0.88", "3.12.0", "3.1.6", "1.4.0")
)
