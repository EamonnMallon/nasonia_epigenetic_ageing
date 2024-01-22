library(tidyverse)
library(viridis)
library(ggpubr)
library(ggfortify)
library(emmeans)
set.seed(123)
#####################################################
#####################################################





####Getting in the cpgs
kristi_data <- read.delim("~/Dropbox/Projects/kristi_epigenetic_clock/nasonia_methylkitobject_percent.txt")
kristi_data %>% mutate(across(F01:M163, ~ .x/100))->kristi_data # turning into betas
kristi_data$cpg<- 1:375403
kristi_data$id<-paste(kristi_data$chr,kristi_data$start,kristi_data$end)
#####################################################
#####################################################

####Identifying the significant ones
temp = list.files(path = "/Users/emb3/CpG_meth/significant/", pattern="*.csv", full.names = TRUE)

list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
         read.csv), envir = .GlobalEnv)

signif<- rbind(X.Users.emb3.CpG_meth.significant..F0_vs_F16_15._sig,X.Users.emb3.CpG_meth.significant..F0_vs_F8_Hypermethylated_w15.diff,X.Users.emb3.CpG_meth.significant..F0_vs_F8_Hypomethylated_w15.diff,X.Users.emb3.CpG_meth.significant..F8_vs_F16_15._sig,X.Users.emb3.CpG_meth.significant..M0_vs_M16_15._sig, X.Users.emb3.CpG_meth.significant..M0_vs_M8_15._sig, X.Users.emb3.CpG_meth.significant..M8_vs_M16_15._sig)
signif$id<-paste(signif$chr,signif$start,signif$end)
signif<-select(signif,id)
signif<-unique(signif)


kristi_data <- merge(kristi_data,signif,by="id")
#####################################################
#####################################################



####### Wrangling the data ######################

entropy_data_loci<-select(kristi_data,c(1,6:23)) # just diff cpgs
#entropy_data_loci<-select(kristi_data,c(23,5:21,22)) # all meth cpgs

entropy_data<- gather(entropy_data_loci, sample, betavalue, F01:M163, factor_key=TRUE)
entropy_data$sex<- substr(entropy_data$sample,1,1)
entropy_data$sex<- as.factor(entropy_data$sex)
entropy_data$chrono<- substr(entropy_data$sample,2,2)
entropy_data$chrono<-as.numeric(recode(entropy_data$chrono, "1" = "16"))


beta_normalize <- function(x) {
  x_ <- ((x - min(x)) / (max(x) - min(x)))
  (x_ * (length(x_) - 1) + 0.5) / length(x_)
} #https://static1.squarespace.com/static/58a7d1e52994ca398697a621/t/5a2ebc43e4966b0fab6b02de/1513012293857/betareg_politics.pdf

entropy_data$betavalue<-beta_normalize(entropy_data$betavalue)
#####################################################
#####################################################







######Calculating entropy ###################################
entropy_factor<-(1/(length(unique(entropy_data$id))))*log2(0.5)
#entropy_factor<-0.1
resid_entropy_values = entropy_data %>% group_by(sample)  %>%
  summarise(entropy_pre = 
            sum(
              (betavalue*log2(betavalue))
              +((1-betavalue)*log2((1-betavalue))
              )),
            .groups = 'drop') # equation from https://www-sciencedirect-com.ezproxy3.lib.le.ac.uk/science/article/pii/S1097276512008933#sec3
# details on summarise https://sparkbyexamples.com/r-programming/group-by-summarise-in-r/
resid_entropy_values$entropy_measure<-entropy_factor*resid_entropy_values$entropy_pre
resid_entropy_values$chrono<-as.numeric(c(0,0,0,8,8,8,16,16,16,0,0,0,8,8,16,16,16)) #should this be numeric? " If you are just interested in fitting longitudinal data over a limited number of times of special interest, then fitting them as factors would be preferred, because then you could explain which time points had an actual effect on the response. However, if you have several time points and only want to know what the general trend is, it would be better to fit the data as numeric." I think because there are so few time points, it makes sense to treat it categorically. But the factor analyis is not very robust, results start disappearing when you look at it pairwise. 
resid_entropy_values$sex<-factor(c("Females","Females","Females","Females","Females","Females","Females","Females","Females","Males","Males","Males","Males","Males","Males","Males","Males"))
#####################################################
#####################################################







###Beta regression#############################################
library(betareg)
entropy_model<-betareg(resid_entropy_values$entropy_measure~resid_entropy_values$sex*resid_entropy_values$chrono)
summary(entropy_model)
library(car) # does this interfer with recode: yes, put it here where its needed
Anova(entropy_model)
detach("package:car", unload = TRUE) # just so it doesn't mess with recode
#####################################################
#####################################################



### Post hoc tests #####################################################
emmip(entropy_model, sex ~ chrono,CIs = TRUE, 
      cov.reduce = FALSE) 
joint_tests(entropy_model, by = "sex")
#emmM1 = emmeans(entropy_model, ~ sex * chrono)
#emmM1
#pairs(emmM1, by = "sex")
#####################################################
#####################################################




###### Producing figure for paper###########################################
#ggboxplot(resid_entropy_values, x = "chrono", y = "entropy_measure",
  #        color = "sex",  palette = "viridis",
   #       add = "jitter", ylab = "Entropy", xlab = "Chronological age (Days)") + theme(legend.title = element_blank(),legend.position = c(0.1, 0.8))  

#ggsave("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/entropy_scatter.pdf")




ggplot(resid_entropy_values, aes(x = chrono, y = entropy_measure, color =sex)) +
  geom_point()+
  geom_line(data = cbind(resid_entropy_values, pred = predict(entropy_model)), aes(y = pred)) +
  ylab("Entropy")+
  xlab("Chronological age (Days)") + 
  theme_classic(base_size = 16) +
  theme(legend.title = element_blank(),legend.position = c(0.1, 0.8))

ggsave("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/entropy_regression_scatter.pdf",width = 10.9, height = 4.5, units = "in")
 
#####################################################
#####################################################