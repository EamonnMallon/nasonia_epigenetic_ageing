EnsemblBetaVMPsSIG_final <- read.csv("~/Dropbox/Projects/bernice_vmp/MiniProject_BCh/FINALData/EnsemblBetaVMPsSIG_final.csv")

EnsemblBetaVMPsSIG<-  EnsemblBetaVMPsSIG_final[,-c(7)]

EnsemblBetaVMPsSIG<-
  EnsemblBetaVMPsSIG[!duplicated(EnsemblBetaVMPsSIG[,c('ensembl_gene_id')]),]



write.csv(EnsemblBetaVMPsSIG, "/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/supplemental_table_5_vmps.csv")
