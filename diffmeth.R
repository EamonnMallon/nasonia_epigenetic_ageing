####Identifying the significant ones############
#################################################


library(dplyr)
library(ggplot2)
library(reshape2)


temp = list.files(path = "/Users/emb3/CpG_meth/significant/", pattern="*.csv", full.names = TRUE) # in mac


list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
         read.csv), envir = .GlobalEnv)



signif <- bind_rows(list(F016 = X.Users.emb3.CpG_meth.significant..F0_vs_F16_15._sig, 
                         F08=X.Users.emb3.CpG_meth.significant..F0_vs_F8_Hypermethylated_w15.diff, 
                         F08=X.Users.emb3.CpG_meth.significant..F0_vs_F8_Hypomethylated_w15.diff, 
                         m016=X.Users.emb3.CpG_meth.significant..M0_vs_M16_15._sig, 
                         m08=X.Users.emb3.CpG_meth.significant..M0_vs_M8_15._sig), 
                         #m0816=X.Users.emb3.CpG_meth.significant..M8_vs_M16_15._sig), 
                         .id = 'source') #f816 is empty, took m0816 out

signif$id<-paste(signif$chr,signif$start,signif$end)


# No point doing 08,816 common, there are no 816 females. Below I subsetted those sign. in 16 or 8 and then only kept those who were in both


females<-subset(signif,source == "F016" | source == "F08")
males<-subset(signif,source == "m016"| source == "m08")

commonfemales<- females[duplicated(females$id)|duplicated(females$id, fromLast=TRUE),]
commonmales<- males[duplicated(males$id)|duplicated(males$id, fromLast=TRUE),]






#################################################################
#Although required to check that hypo stayed hypo etc can ignore this bit. They all did, so fine to used dataframes commonfemales and common males.

### Making sure the pattern is the same for females
com_female_short<- subset(commonfemales, select = c(source,id,meth.diff) )
commonfemale_wide <- dcast(com_female_short, id ~ source, value.var="meth.diff")
commonfemale_wide$sign16<-sign(commonfemale_wide$F016)
commonfemale_wide$sign8<-sign(commonfemale_wide$F08)
commonfemale_wide$pattern<-commonfemale_wide$sign16+commonfemale_wide$sign8
commonfemale_wide_final <- filter(commonfemale_wide, pattern != 0)

ggplot(commonfemales, aes(source, meth.diff)) +
  geom_point() +
  facet_wrap(~id)

### Making sure the pattern is the same for males
com_male_short<- subset(commonmales, select = c(source,id,meth.diff) )
commonmale_wide <- dcast(com_male_short, id ~ source, value.var="meth.diff")
commonmale_wide$sign16<-sign(commonmale_wide$m016)
commonmale_wide$sign8<-sign(commonmale_wide$m08)
commonmale_wide$pattern<-commonmale_wide$sign16+commonmale_wide$sign8
commonmale_wide_final <- filter(commonmale_wide, pattern != 0)

ggplot(commonmales, aes(source, meth.diff)) +
  geom_point() +
  facet_wrap(~id)
##### Graph looks good but need to make sure that they stay hypo hyper.

########################################################################
########################################################################


hyperfemales<-subset(commonfemales, meth.diff > 0)
hypofemales<-subset(commonfemales, meth.diff < 0)
hypermales<-subset(commonmales, meth.diff > 0)
hypomales<-subset(commonmales, meth.diff < 0)


###################Testing for overlap#####################
hyper<-merge(hyperfemales,hypermales, by = "id")
hypo<-merge(hypofemales,hypomales, by = "id")
opp1<-merge(hypofemales,hypermales, by = "id")
opp2<-merge(hyperfemales,hypomales, by = "id")

###############Annotating the genes of interest########################
#https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html


##GO terms not directly used in hypergeometric test , but need LOCs


############################################

# Setting up Ensembl database
library(biomaRt)
ensembl = useEnsemblGenomes(biomart="metazoa_mart", dataset="nvitripennis_eg_gene")
attributes <- listAttributes(ensembl) #all the columns you can get
emsemblgenes <- getBM(attributes=c('ensembl_gene_id',
                                   "description",
                                   "chromosome_name",
                                   "start_position",
                                   "end_position",
                                   "strand",
                                   "go_id"), 
                      mart=ensembl)
ensemblgenes <- as.data.frame(emsemblgenes)
###################################################################



##Hyperfemales 
# Changing chr name so it's compatible with Ensembl data 
Ensembl_chr <- hyperfemales$chr
for(i in c(1:nrow(hyperfemales))){
  if(Ensembl_chr[i] == "NC_045757.1"){
    Ensembl_chr[i] <- "CM020934.1"
  }else if(Ensembl_chr[i] == "NC_045758.1"){
    Ensembl_chr[i] <- "CM020935.1"
  }else if(Ensembl_chr[i] == "NC_045759.1"){
    Ensembl_chr[i] <- "CM020936.1"
  }else if(Ensembl_chr[i] == "NC_045760.1"){
    Ensembl_chr[i] <- "CM020937.1"
  }else if(Ensembl_chr[i] == "NC_045761.1"){
    Ensembl_chr[i] <- "CM020938.1"
  } else {
    Ensembl_chr[i] <- "WELF01000060.1"
  }
}
hyperfemales <- cbind(hyperfemales, Ensembl_chr)
hyperfemales %>% distinct(id, .keep_all = TRUE) -> hyperfemales



# Annotating hyperfemales with Ensembl database
ensemblcompare <- hyperfemales[, c(11,4,5)]
ensemblres <- data.frame(matrix(nrow=nrow(ensemblcompare), ncol=ncol(ensemblgenes)))
colnames(ensemblres) <- c("ensembl_gene_id",
                          "description",
                          "chromosome_name",
                          "start_position",
                          "end_position",
                          "strand",
                          "go_id")

for(i in 1:nrow(hyperfemales)){ 
  for(j in 1:nrow(ensemblgenes)){
    if(ensemblcompare$Ensembl_chr[i] == ensemblgenes$chromosome_name[j] &&
       ensemblcompare$start[i] >= ensemblgenes$start_position[j] &&
       ensemblcompare$end[i] <= ensemblgenes$end_position[j]){
      ensemblres[i, ] <- ensemblgenes[j, ]
    }
  }
  # Print progress 
  cat("\r",round(i/nrow(ensemblres)*100),'% done')
}

ensemblres %>% distinct(ensembl_gene_id, .keep_all = TRUE) -> hyperfemaleensemblres
write.csv(hyperfemaleensemblres, "/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/diffmeth_GO_analysis/annotatedhyperfem.csv")

##############################################################################
##Hypofemales 
# Changing chr name so it's compatible with Ensembl data 
Ensembl_chr <- hypofemales$chr
for(i in c(1:nrow(hypofemales))){
  if(Ensembl_chr[i] == "NC_045757.1"){
    Ensembl_chr[i] <- "CM020934.1"
  }else if(Ensembl_chr[i] == "NC_045758.1"){
    Ensembl_chr[i] <- "CM020935.1"
  }else if(Ensembl_chr[i] == "NC_045759.1"){
    Ensembl_chr[i] <- "CM020936.1"
  }else if(Ensembl_chr[i] == "NC_045760.1"){
    Ensembl_chr[i] <- "CM020937.1"
  }else if(Ensembl_chr[i] == "NC_045761.1"){
    Ensembl_chr[i] <- "CM020938.1"
  } else {
    Ensembl_chr[i] <- "WELF01000060.1"
  }
}
hypofemales <- cbind(hypofemales, Ensembl_chr)
hypofemales %>% distinct(id, .keep_all = TRUE) -> hypofemales



# Annotating hypofemales with Ensembl database
ensemblcompare <- hypofemales[, c(11,4,5)]
ensemblres <- data.frame(matrix(nrow=nrow(ensemblcompare), ncol=ncol(ensemblgenes)))
colnames(ensemblres) <- c("ensembl_gene_id",
                          "description",
                          "chromosome_name",
                          "start_position",
                          "end_position",
                          "strand",
                          "go_id")

for(i in 1:nrow(hypofemales)){ 
  for(j in 1:nrow(ensemblgenes)){
    if(ensemblcompare$Ensembl_chr[i] == ensemblgenes$chromosome_name[j] &&
       ensemblcompare$start[i] >= ensemblgenes$start_position[j] &&
       ensemblcompare$end[i] <= ensemblgenes$end_position[j]){
      ensemblres[i, ] <- ensemblgenes[j, ]
    }
  }
  # Print progress 
  cat("\r",round(i/nrow(ensemblres)*100),'% done')
}

ensemblres %>% distinct(ensembl_gene_id, .keep_all = TRUE) -> hypofemaleensemblres

write.csv(hypofemaleensemblres, "/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/diffmeth_GO_analysis/annotatedhypofem.csv")
##############################################################################


##############################################################################
##Hypomales 
# Changing chr name so it's compatible with Ensembl data 
Ensembl_chr <- hypomales$chr
for(i in c(1:nrow(hypomales))){
  if(Ensembl_chr[i] == "NC_045757.1"){
    Ensembl_chr[i] <- "CM020934.1"
  }else if(Ensembl_chr[i] == "NC_045758.1"){
    Ensembl_chr[i] <- "CM020935.1"
  }else if(Ensembl_chr[i] == "NC_045759.1"){
    Ensembl_chr[i] <- "CM020936.1"
  }else if(Ensembl_chr[i] == "NC_045760.1"){
    Ensembl_chr[i] <- "CM020937.1"
  }else if(Ensembl_chr[i] == "NC_045761.1"){
    Ensembl_chr[i] <- "CM020938.1"
  } else {
    Ensembl_chr[i] <- "WELF01000060.1"
  }
}
hypomales <- cbind(hypomales, Ensembl_chr)
hypomales %>% distinct(id, .keep_all = TRUE) -> hypomales



# Annotating hypofemales with Ensembl database
ensemblcompare <- hypomales[, c(11,4,5)]
ensemblres <- data.frame(matrix(nrow=nrow(ensemblcompare), ncol=ncol(ensemblgenes)))
colnames(ensemblres) <- c("ensembl_gene_id",
                          "description",
                          "chromosome_name",
                          "start_position",
                          "end_position",
                          "strand",
                          "go_id")

for(i in 1:nrow(hypomales)){ 
  for(j in 1:nrow(ensemblgenes)){
    if(ensemblcompare$Ensembl_chr[i] == ensemblgenes$chromosome_name[j] &&
       ensemblcompare$start[i] >= ensemblgenes$start_position[j] &&
       ensemblcompare$end[i] <= ensemblgenes$end_position[j]){
      ensemblres[i, ] <- ensemblgenes[j, ]
    }
  }
  # Print progress 
  cat("\r",round(i/nrow(ensemblres)*100),'% done')
}

ensemblres %>% distinct(ensembl_gene_id, .keep_all = TRUE) -> hypomaleensemblres

write.csv(hypomaleensemblres, "/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/diffmeth_GO_analysis/annotatedhypomal.csv")
##############################################################################




##############################################################################
##Hypermales 
# Changing chr name so it's compatible with Ensembl data 
Ensembl_chr <- hypermales$chr
for(i in c(1:nrow(hypermales))){
  if(Ensembl_chr[i] == "NC_045757.1"){
    Ensembl_chr[i] <- "CM020934.1"
  }else if(Ensembl_chr[i] == "NC_045758.1"){
    Ensembl_chr[i] <- "CM020935.1"
  }else if(Ensembl_chr[i] == "NC_045759.1"){
    Ensembl_chr[i] <- "CM020936.1"
  }else if(Ensembl_chr[i] == "NC_045760.1"){
    Ensembl_chr[i] <- "CM020937.1"
  }else if(Ensembl_chr[i] == "NC_045761.1"){
    Ensembl_chr[i] <- "CM020938.1"
  } else {
    Ensembl_chr[i] <- "WELF01000060.1"
  }
}
hypermales <- cbind(hypermales, Ensembl_chr)
hypermales %>% distinct(id, .keep_all = TRUE) -> hypermales



# Annotating hyperfemales with Ensembl database
ensemblcompare <- hypermales[, c(11,4,5)]
ensemblres <- data.frame(matrix(nrow=nrow(ensemblcompare), ncol=ncol(ensemblgenes)))
colnames(ensemblres) <- c("ensembl_gene_id",
                          "description",
                          "chromosome_name",
                          "start_position",
                          "end_position",
                          "strand",
                          "go_id")

for(i in 1:nrow(hypermales)){ 
  for(j in 1:nrow(ensemblgenes)){
    if(ensemblcompare$Ensembl_chr[i] == ensemblgenes$chromosome_name[j] &&
       ensemblcompare$start[i] >= ensemblgenes$start_position[j] &&
       ensemblcompare$end[i] <= ensemblgenes$end_position[j]){
      ensemblres[i, ] <- ensemblgenes[j, ]
    }
  }
  # Print progress 
  cat("\r",round(i/nrow(ensemblres)*100),'% done')
}

ensemblres %>% distinct(ensembl_gene_id, .keep_all = TRUE) -> hypermaleensemblres

write.csv(hypermaleensemblres, "/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/diffmeth_GO_analysis/annotatedhypermal.csv")
##############################################################################





##############################################################################
##All sign diff

totalsign<- rbind(X.Users.emb3.CpG_meth.significant..F0_vs_F16_15._sig,X.Users.emb3.CpG_meth.significant..F0_vs_F8_Hypermethylated_w15.diff,X.Users.emb3.CpG_meth.significant..F0_vs_F8_Hypomethylated_w15.diff,X.Users.emb3.CpG_meth.significant..F8_vs_F16_15._sig,X.Users.emb3.CpG_meth.significant..M0_vs_M16_15._sig, X.Users.emb3.CpG_meth.significant..M0_vs_M8_15._sig, X.Users.emb3.CpG_meth.significant..M8_vs_M16_15._sig)
totalsign$id<-paste(totalsign$chr,totalsign$start,totalsign$end)


# Changing chr name so it's compatible with Ensembl data 
Ensembl_chr <- totalsign$chr
for(i in c(1:nrow(totalsign))){
  if(Ensembl_chr[i] == "NC_045757.1"){
    Ensembl_chr[i] <- "CM020934.1"
  }else if(Ensembl_chr[i] == "NC_045758.1"){
    Ensembl_chr[i] <- "CM020935.1"
  }else if(Ensembl_chr[i] == "NC_045759.1"){
    Ensembl_chr[i] <- "CM020936.1"
  }else if(Ensembl_chr[i] == "NC_045760.1"){
    Ensembl_chr[i] <- "CM020937.1"
  }else if(Ensembl_chr[i] == "NC_045761.1"){
    Ensembl_chr[i] <- "CM020938.1"
  } else {
    Ensembl_chr[i] <- "WELF01000060.1"
  }
}
totalsign <- cbind(totalsign, Ensembl_chr)
totalsign %>% distinct(id, .keep_all = TRUE) -> totalsign



# Annotating total with Ensembl database
ensemblcompare <- totalsign[, c(10,3,4)]
ensemblres <- data.frame(matrix(nrow=nrow(ensemblcompare), ncol=ncol(ensemblgenes)))
colnames(ensemblres) <- c("ensembl_gene_id",
                          "description",
                          "chromosome_name",
                          "start_position",
                          "end_position",
                          "strand",
                          "go_id")

for(i in 1:nrow(totalsign)){ 
  for(j in 1:nrow(ensemblgenes)){
    if(ensemblcompare$Ensembl_chr[i] == ensemblgenes$chromosome_name[j] &&
       ensemblcompare$start[i] >= ensemblgenes$start_position[j] &&
       ensemblcompare$end[i] <= ensemblgenes$end_position[j]){
      ensemblres[i, ] <- ensemblgenes[j, ]
    }
  }
  # Print progress 
  cat("\r",round(i/nrow(ensemblres)*100),'% done')
}

ensemblres %>% distinct(ensembl_gene_id, .keep_all = TRUE) -> totalsigneensemblres

write.csv(totalsigneensemblres, "/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/diffmeth_GO_analysis/annotatedtotal.csv")
##############################################################################




######
#Setting Up
#geneCount <- read.csv("/Users/u2240321/Documents/MiniProject1/Data/FINAL/Gene_res_final.csv")

GO_annotations <- read.table("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/Nasonia_PSR1.1_Gene_GO.txt")
GO_annotations <- GO_annotations[-1,]
colnames(GO_annotations) <- c("Gene_ID","GO_Terms")




#############################################
########Alun's script (GOStats-based)########
library(GOstats)
library(GSEABase)

# Re-sorting GO_annotations
evi <- c(rep("IEA",nrow(GO_annotations))) # Evidence code = IEA = Inferred from Electronic Annotation
GO_annotations <- cbind(GO_annotations$GO_Terms,evi,GO_annotations$Gene_ID)
colnames(GO_annotations) <- c("GOIds","evi","genes")
GO_annotations <- as.data.frame(GO_annotations)


#Only comparing to diffmeth genes############################################################################

GO_annotations %>% 
  filter(genes %in% totalsigneensemblres$ensembl_gene_id)->GO_annotations # much cooler way of what I usually do by merge


########################################################################
#Setting Up
GO_frame <- GOFrame(GO_annotations,organism = "Nasonia Vitripennis")
goAllFrame <- GOAllFrame(GO_frame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.vector(unique(GO_annotations[,3]))

# AT the minute comparing to all genes. To follow what we are doing I should subset it to the 5000 that are in any way diff methed. Not sure this makes sense to me. Talk to Hollie.---Thats what you did in line 404


########################################################################

#####Hyper females


# Genes & GO Terms only as vector
Gene <- na.omit(hyperfemaleensemblres$ensembl_gene_id)


# Count for gene
genelist <- unique(Gene)
genelist_count <- numeric(length(genelist))
for(i in 1:length(genelist)){
  for(j in 1:length(Gene)){
    cur_gene <- genelist[i]
    cur_dum <- Gene[j]
    if(cur_gene == cur_dum){
      genelist_count[i] <- genelist_count[i] + 1
    }
  }
}
Gene_res <- cbind(genelist,genelist_count)


my_genes <- as.vector(Gene_res[,1])
my_genes <- my_genes[my_genes %in% universe] #Getting rid of genes that are actually not in the annotation file
length(my_genes) #Only 1887/2451 genes are actually annotated (final)

#Setting up parameters
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

#Hypergeometry Test
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
Result <- summary(GO_enrichment[["BP_over"]])

#FDR Correction
Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

#REVIGO Write-Out
REVIGO <- Result_FDR[,1:2]
write.table(REVIGO,"/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/hyperfemale_diffmeth_revigo.txt",row.names = F,sep = "\t",quote = F)
##########################################################################



#####Hypo females


# Genes & GO Terms only as vector
Gene <- na.omit(hypofemaleensemblres$ensembl_gene_id)


# Count for gene
genelist <- unique(Gene)
genelist_count <- numeric(length(genelist))
for(i in 1:length(genelist)){
  for(j in 1:length(Gene)){
    cur_gene <- genelist[i]
    cur_dum <- Gene[j]
    if(cur_gene == cur_dum){
      genelist_count[i] <- genelist_count[i] + 1
    }
  }
}
Gene_res <- cbind(genelist,genelist_count)


my_genes <- as.vector(Gene_res[,1])
my_genes <- my_genes[my_genes %in% universe] #Getting rid of genes that are actually not in the annotation file
length(my_genes) #Only 1887/2451 genes are actually annotated (final)

#Setting up parameters
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

#Hypergeometry Test
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
Result <- summary(GO_enrichment[["BP_over"]])

#FDR Correction
Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

#REVIGO Write-Out
REVIGO <- Result_FDR[,1:2]
write.table(REVIGO,"/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/hypofemale_diffmeth_revigo.txt",row.names = F,sep = "\t",quote = F)

#####Hypo males


# Genes & GO Terms only as vector
Gene <- na.omit(hypomaleensemblres$ensembl_gene_id)


# Count for gene
genelist <- unique(Gene)
genelist_count <- numeric(length(genelist))
for(i in 1:length(genelist)){
  for(j in 1:length(Gene)){
    cur_gene <- genelist[i]
    cur_dum <- Gene[j]
    if(cur_gene == cur_dum){
      genelist_count[i] <- genelist_count[i] + 1
    }
  }
}
Gene_res <- cbind(genelist,genelist_count)


my_genes <- as.vector(Gene_res[,1])
my_genes <- my_genes[my_genes %in% universe] #Getting rid of genes that are actually not in the annotation file
length(my_genes) #Only 1887/2451 genes are actually annotated (final)

#Setting up parameters
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

#Hypergeometry Test
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
Result <- summary(GO_enrichment[["BP_over"]])

#FDR Correction
Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

#REVIGO Write-Out
REVIGO <- Result_FDR[,1:2]
write.table(REVIGO,"/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/hypomale_diffmeth_revigo.txt",row.names = F,sep = "\t",quote = F)



#####Hyper males


# Genes & GO Terms only as vector
Gene <- na.omit(hypermaleensemblres$ensembl_gene_id)


# Count for gene
genelist <- unique(Gene)
genelist_count <- numeric(length(genelist))
for(i in 1:length(genelist)){
  for(j in 1:length(Gene)){
    cur_gene <- genelist[i]
    cur_dum <- Gene[j]
    if(cur_gene == cur_dum){
      genelist_count[i] <- genelist_count[i] + 1
    }
  }
}
Gene_res <- cbind(genelist,genelist_count)


my_genes <- as.vector(Gene_res[,1])
my_genes <- my_genes[my_genes %in% universe] #Getting rid of genes that are actually not in the annotation file
length(my_genes) #Only 1887/2451 genes are actually annotated (final)

#Setting up parameters
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

#Hypergeometry Test
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
Result <- summary(GO_enrichment[["BP_over"]])

#FDR Correction
Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

#REVIGO Write-Out
REVIGO <- Result_FDR[,1:2]
write.table(REVIGO,"/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/long_form_paper/hypermale_diffmeth_revigo.txt",row.names = F,sep = "\t",quote = F)


