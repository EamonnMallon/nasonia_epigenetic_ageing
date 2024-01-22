library(GOstats)
library(GSEABase)
library(treemap)
library(dplyr)


#---------------------------------------------------
#---------------------------------------------------
#---------------------------------------------------
# 1.Creating and formatting background gene sets for GO analysis
#---------------------------------------------------
#---------------------------------------------------
#---------------------------------------------------

#---------------------------------------------------
#---------------------------------------------------
# 1.1 Getting in all methylated CpGs in ageing dataset
#---------------------------------------------------
#---------------------------------------------------
kristi_data <- read.delim("~/Dropbox/Projects/kristi_epigenetic_clock/nasonia_methylkitobject_percent.txt")
kristi_data %>% mutate(across(F01:M163, ~ .x/100))->kristi_data # turning into betas
kristi_data$cpg<- 1:375403
kristi_data$id<-paste(kristi_data$chr,kristi_data$start,kristi_data$end)

## ----------------------------------------------------------------
# 1.1.1 Annotating all methylated CpGs with gene ID (from Hollie)
## ----------------------------------------------------------------

library(sqldf)
library(readr)

# All annotations, including numbered exons and introns
annotation <- read_delim("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/GCF_009193385.2_Nvit_psr_1.1_genomic_numbered_exons.txt", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

# CpG sites (chr, cpg_position) 
select(kristi_data, chr,start)%>%rename(cpg_position = start)->cpg_sites

output <- sqldf("SELECT sample.chr,
                    sample.cpg_position,
                    annot.chr,
                    annot.feature,
                    annot.start,
                    annot.end,
                    annot.gene_id,
                    annot.number
                    FROM cpg_sites AS sample
                    LEFT JOIN annotation AS annot
                    ON sample.chr = annot.chr
                    AND (sample.cpg_position >= annot.start AND sample.cpg_position <= annot.end)")

write_csv(output, "~/Dropbox/Projects/kristi_epigenetic_clock/methylatedcpgs_with_LOC.csv") # just a temp so I don't have to rerun sqldf

temp<-data.frame(output)
background_genes <- filter(temp, feature == "gene")
background_genes<-unique(background_genes[c("gene_id")])
names(background_genes) <- c("genes")

## ---------------------------------------------------------------
# GO terms associated with geneIDs
## ----------------------------------------------------------------
GO_annotations <- read.table("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/scripts/Nasonia_PSR1.1_Gene_GO.txt")

GO_annotations[,3] <- paste("IEA")
names(GO_annotations) <- c("genes","GOIds","evi")
GO_annotations[,3] <- paste("IEA")
GO_annotations <- GO_annotations[c(2,3,1)]
GO_annotations <- GO_annotations[-1,] # remove the original


GO_annotations<-merge(GO_annotations,background_genes, by="genes")# after all that only changes it from 771163 to 699047


write_csv(GO_annotations, "~/Dropbox/Projects/kristi_epigenetic_clock/backgroundGOannotations.csv") 













