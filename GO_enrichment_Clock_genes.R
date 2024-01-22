#-----------------------------------------------
# Script created by Alun Jones, see paper Bebane et al. (2019) Neonics and bumblebees...
#-----------------------------------------------

library(GOstats)
library(GSEABase)
library(treemap)
library(readr)

#-----------------------------------------------
# Read in background GO set and make compatible with GOstats


GO_annotations <- read_csv("~/Dropbox/Projects/kristi_epigenetic_clock/backgroundGOannotations.csv")
GO_annotations <- data.frame(GO_annotations[c(2,3,1)])#need to make an old fashion dataframe otherwise GOFrame doesn't work. I think because I used read_csv



#-----------------------------------------------
# Create necessary objects

GO_frame <- GOFrame(GO_annotations,organism = "N.vitripennis")
goAllFrame <- GOAllFrame(GO_frame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.vector(unique(GO_annotations[,3]))

#-----------------------------------------------
# Read in gene's of interest 


my_genes <- read.csv("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/clock_genes.csv")

#my_genes <- as.data.frame(na.omit(my_genes$geneID)) why do I need this?
my_genes = subset(my_genes, select = c(Symbol) )
colnames(my_genes) <- "genes"
my_genes <- as.vector(my_genes[,1])

#-----------------------------------------------
# Keep only genes with annotated GOs

my_genes <- my_genes[my_genes %in% universe]
length(my_genes)


#-----------------------------------------------
# Set up paramters for hypergeometric test

Get_GO_params_all <- function(genes_of_i,universe,pvalue_cut){
  onto_terms <- c("BP","CC","MF")
  directions <- c("over","under")
  param_list <- list()
  name_1 <- list()
  for(i in 1:3){
    for(j in 1:2){
      name_1 <- c(name_1,paste(onto_terms[i],directions[j],sep = "_"))
      parameters <- GSEAGOHyperGParams(name="Hygeo params",
                                       geneSetCollection = gsc,
                                       universeGeneIds = universe,
                                       geneIds = my_genes,
                                       ontology = paste(onto_terms[i]),
                                       pvalueCutoff = pvalue_cut,
                                       conditional = T,testDirection = paste(directions[j]))
      param_list <- c(param_list,parameters)
    }
  }
  names(param_list) <- name_1
  return(param_list)
}

param_list <- Get_GO_params_all(genes_of_i = DE_Genes_A,universe = universe,
                                pvalue_cut = 0.05)


#-----------------------------------------------
# Hypergeometric test

Hyper_G_test <- function(param_list){
  Hyper_G_list <- list()
  for(i in 1:length(param_list)){
    res <- hyperGTest(param_list[[i]])
    Hyper_G_list <- c(Hyper_G_list,res)
  }
  names(Hyper_G_list) <- names(param_list)
  return(Hyper_G_list)
}


GO_enrichment <- Hyper_G_test(param_list = param_list)


#-----------------------------------------------
# Choose what you want to look for, here the choice is biological process over-represented 

Result <- summary(GO_enrichment[["BP_over"]])

#-----------------------------------------------
# FDR correction

Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

#-----------------------------------------------
# Make an output for REVIGO and write out

REVIGO <- Result_FDR[,1:2]

#write.table(REVIGO,"./Levels_FDR0.05/No_meth/A_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.05/No_meth/B_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.05/No_meth/C_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.05/No_meth/D_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.05/No_meth/E_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.05/No_meth/F_no_meth_against_all_genes_GOs.txt",row.names = F,sep = "\t",quote = F)

#write.table(REVIGO,"./Levels_FDR0.05/High_meth/A_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.05/High_meth/B_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.05/High_meth/C_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.05/High_meth/D_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./Levels_FDR0.05/High_meth/E_high_meth_against_all_meth_genes_GOs.txt",row.names = F,sep = "\t",quote = F)
write.table(REVIGO,"/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/GO_analysis/clock_GOs.txt",row.names = F,sep = "\t",quote = F)
