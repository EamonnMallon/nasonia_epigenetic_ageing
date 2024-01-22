#### Script for An epigenetic clock in an insect model system
library(dplyr)
library(purrr)
library(tidyverse)
library(ggfortify)
library(glmnet)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(vegan)
library(modeest)
library(Metrics)
# Model interpretability packages
library(vip)      # for variable importance
library(viridis)
library(broom)
library(lme4)
library(emmeans)
set.seed(123)



####Getting in the cpgs##########
##################################
kristi_data <- read.delim("~/Dropbox/Projects/kristi_epigenetic_clock/nasonia_methylkitobject_percent.txt")
kristi_data %>% mutate(across(F01:M163, ~ .x/100))->kristi_data # turning into betas
kristi_data$cpg<- 1:375403
kristi_data$id<-paste(kristi_data$chr,kristi_data$start,kristi_data$end)
###############################################################################
###############################################################################

####Calculating average methylation for each sample across time##############
#############################################################################
transpose_df<-t(subset(kristi_data,select = -c(chr,start,end,strand,cpg,id)))
tr_averages<- as.data.frame(matrix(ncol=2, nrow=17))
tr_averages$means<-rowMeans(transpose_df)
tr_averages$samples<-rownames(transpose_df)
tr_averages$chrono<-c(0,0,0,8,8,8,16,16,16,0,0,0,8,8,16,16,16)
tr_averages$sex<-c("f","f","f","f","f","f","f","f","f","m","m","m","m","m","m","m","m")

####Using a mixed effect model to check if genome is going hypo or hyper with ageing####
########################################################################################

kristi_long<- subset(kristi_data,select = -c(chr,start,end,strand,cpg))
kristi_long <- gather(kristi_long, sample, beta, F01:M163, factor_key=TRUE)
kristi_long$sex<- substr(kristi_long$sample,1,1)
kristi_long$chrono<- substr(kristi_long$sample,2,2)
kristi_long %>% mutate(chrono=recode(chrono, '1'='16'))->kristi_long
kristi_long$chrono<-as.numeric(kristi_long$chrono)
mixed.lmer <- lmer(beta ~ chrono*sex + (1|id), data = kristi_long)
summary(mixed.lmer)
sex.lmer <- lmer(beta ~ chrono + (1|id), data = kristi_long)
age.lmer <- lmer(beta ~ sex + (1|id), data = kristi_long)
interaction.lmer <- lmer(beta ~ chrono+sex + (1|id), data = kristi_long)
anova(sex.lmer,mixed.lmer)
anova(age.lmer,mixed.lmer)
anova(interaction.lmer,mixed.lmer)
emmip(mixed.lmer, sex ~ chrono, CIs = TRUE, 
      cov.reduce = FALSE, xlab = "Chronological age (Days)", ylab = "Methylation proportion per CpG") +
  theme_classic(base_size = 16) +
  # Move legend into plot to fill empty space
  theme(legend.title = element_blank(),legend.position = c(.1, .8)) +
  scale_color_discrete(labels = c("Females", "Males"))#https://yuzar-blog.netlify.app/posts/2022-12-29-emmeans2interactions/
ggsave("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/methylome.pdf",width = 10.9, height = 4.5, units = "in")#,width = 10.9, height = 4, units = "in"


#ggplot(kristi_long, aes(x = chrono, y = beta, color =sex)) +
 # geom_point()+
 # geom_line(data = cbind(kristi_long, pred = predict(mixed.lmer)), aes(y = pred)) +
 # ylab("Methylation proportion per CpG")+
 # xlab("Chronological age (Days)") + 
 #theme_classic() +
  #theme(legend.title = element_blank(),legend.position = c(0.1, 0.8))

#ggsave("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/methylome2.pdf")



####Identifying the significant ones############
#################################################
temp = list.files(path = "/Users/emb3/CpG_meth/significant/", pattern="*.csv", full.names = TRUE) #mac
#temp = list.files(path = "/home/emb3/Dropbox/Projects/kristi_epigenetic_clock/CpG_meth/significant/", pattern="*.csv", full.names = TRUE) #linux

list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
         read.csv), envir = .GlobalEnv)

#signif<- rbind(X.home.emb3.Dropbox.Projects.kristi_epigenetic_clock.CpG_meth.significant..F0_vs_F16_15._sig,X.home.emb3.Dropbox.Projects.kristi_epigenetic_clock.CpG_meth.significant..F0_vs_F8_Hypermethylated_w15.diff,X.home.emb3.Dropbox.Projects.kristi_epigenetic_clock.CpG_meth.significant..F0_vs_F8_Hypomethylated_w15.diff,X.home.emb3.Dropbox.Projects.kristi_epigenetic_clock.CpG_meth.significant..F8_vs_F16_15._sig,X.home.emb3.Dropbox.Projects.kristi_epigenetic_clock.CpG_meth.significant..M0_vs_M16_15._sig, X.home.emb3.Dropbox.Projects.kristi_epigenetic_clock.CpG_meth.significant..M0_vs_M8_15._sig, X.home.emb3.Dropbox.Projects.kristi_epigenetic_clock.CpG_meth.significant..M8_vs_M16_15._sig)
signif<- rbind(X.Users.emb3.CpG_meth.significant..F0_vs_F16_15._sig,X.Users.emb3.CpG_meth.significant..F0_vs_F8_Hypermethylated_w15.diff,X.Users.emb3.CpG_meth.significant..F0_vs_F8_Hypomethylated_w15.diff,X.Users.emb3.CpG_meth.significant..F8_vs_F16_15._sig,X.Users.emb3.CpG_meth.significant..M0_vs_M16_15._sig, X.Users.emb3.CpG_meth.significant..M0_vs_M8_15._sig, X.Users.emb3.CpG_meth.significant..M8_vs_M16_15._sig)
signif$id<-paste(signif$chr,signif$start,signif$end)
signif<-select(signif,id)
signif<-unique(signif)
############################################################################################
############################################################################################

kristi_data <- merge(kristi_data,signif,by="id")


#### Measures for all differential cpgs##########################
##################################################################
methlevels = NULL
methlevels$zero <- rowMeans(kristi_data[ , c(6,7,8,15,16,17)], na.rm=TRUE)
methlevels$sixteen <- rowMeans(kristi_data[ , c(12,13,14,20,21,22)], na.rm=TRUE)
methlevels$methpattern<-methlevels$sixteen-methlevels$zero
methlevels<-as.data.frame(methlevels)
methlevels$hh<-ifelse(methlevels$methpattern< 0,"hypo","hyper")
table(methlevels$hh)
test<-filter(methlevels, methpattern==0)








play_males<-select(kristi_data, 6:23)
trans_play<-t(play_males)
as.data.frame(trans_play) %>%
  purrr::set_names(as.character(slice(., 18))) %>%
  slice(-18) -> trans_play
#trans_play$samples <- row.names(trans_play)
trans_play$chrono<-c(0,0,0,8,8,8,16,16,16,0,0,0,8,8,16,16,16)
trans_play$sex<-c("f","f","f","f","f","f","f","f","f","m","m","m","m","m","m","m","m")

########################################################
#########################################################



##Elastic net regression
###########################
###########################

x<- as.matrix(trans_play[,1:5290])
y=as.matrix(trans_play$chrono)


# running it a 100 times
lambdas = NULL
for (i in 1:100)
{
  fit <- cv.glmnet(y=y, x=x, alpha=0.5, nfolds=3, relax = TRUE) #alpha =0.5 (elastic net regression)
  errors = data.frame(fit$lambda,fit$cvm)
  lambdas <- rbind(lambdas,errors)
}
#Recall that λ is a tuning parameter that helps to control our model from over-fitting to the training data. To identify the optimal λ value we can use k-fold cross-validation (CV). glmnet::cv.glmnet() can perform k-fold CV, and by default, performs 10-fold CV.

# take mean cvm for each lambda
lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)

# select the best one
plot(fit)
bestindex = which(lambdas[2]==min(lambdas[2]))
bestlambda = lambdas[bestindex,1]




model1 = glmnet(y=y, x=x, lambda=bestlambda)

weights = data.frame(coef.name = dimnames(coef(model1))[[1]], coef.value = matrix(coef(model1)))

weights1 = weights[weights$coef.value!=0,]

a = intersect(names(trans_play),weights1$coef.name) 

data_m_2 <- as.matrix(trans_play[,a])

pred <- data_m_2 %*% weights1$coef.value[2:nrow(weights1)] 

trans_play$pred_age <- pred + weights1$coef.value[1]


names(trans_play)<-make.names(names(trans_play))



graph_data<-select(trans_play,chrono,sex,pred_age)

ggplot(graph_data,aes(y=pred_age,x=chrono, color = sex)) + 
  geom_point(alpha=0.7,size=2,position = position_jitter(width = 0.2, height = 0.2)) + #in paper jitter = 0.4
  xlab("Chronological age (Days)") +
  ylab("Epigenetic age (Days)") +
  scale_y_continuous(breaks=c(0,8,16))+
  scale_x_continuous(breaks=c(0,8,16))+ 
  theme_classic(base_size = 16) +
  # Move legend into plot to fill empty space
  theme(legend.title = element_blank(),legend.position = c(.1, .8)) +
  scale_color_discrete(labels = c("Females", "Males")) #+
  #geom_smooth(method = "lm",  se = FALSE)
 



ggsave("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/epigenetic_clock.pdf", width = 10.9, height = 4.5, units = "in")


final_lm<-lm(chrono~pred_age*sex, data = graph_data)

summary(final_lm)

vip(model1, num_features = 19, geom = "point")

corr<-cor.test(trans_play$pred_age, trans_play$chrono,
               method = "spearman",
               exact = FALSE, conf.level = 0.95, continuity = FALSE)

#corr<-cor.test(data_m$pred_age, data_m$chrono_age,
     #          method = "spearman",
     #          exact = FALSE, conf.level = 0.95, continuity = FALSE)
corr

rmse(trans_play$pred_age, trans_play$chrono)

#-------------------------------------------------------------
# Age acceleration
# ------------------------------------------------------------

# Here I consider two measures of age acceleration.
# The first acceleration meausure is based on the difference.
graph_data$AgeAccelerationDiff=graph_data$pred_age- graph_data$chrono
# The second acceleration meausure equals the residual resulting 
# from regressing DNAmAge on chronological age
graph_data$AgeAccelerationResid= residuals(lm(pred_age~chrono, data = graph_data))
boxplot(AgeAccelerationDiff ~sex,data=graph_data)
boxplot(AgeAccelerationResid ~sex,data=graph_data)



df.summary <- graph_data %>%
  group_by(sex) %>%
  summarise(
    sd = IQR(AgeAccelerationResid, na.rm = TRUE),
    AgeAccelerationResid = median(AgeAccelerationResid)
  )


ggplot(graph_data, aes(sex, AgeAccelerationResid, color = sex)) +
      geom_jitter(
       position = position_jitter(0.2)
       )  + 
     geom_pointrange(
       aes(ymin = AgeAccelerationResid-sd, ymax = AgeAccelerationResid+sd),
       data = df.summary
  ) +theme_classic(base_size = 16) +
  theme(legend.title = element_blank(),legend.position = c(.1, .8)) +
  scale_color_discrete(labels = c("Females", "Males"))+
  ylab("Age acceleration (Days)") +
  scale_x_discrete(breaks=c("f","m"),
  labels=c("Females", "Males"))

ggsave("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/supp_age_accel.pdf", width = 10.9, height = 4.5, units = "in")


wilcox.test(AgeAccelerationResid~sex,data=graph_data)
################################################################################
################################################################################


rename(weights1, cpg = coef.name)->weights1
clock_genes<- merge(weights1,kristi_data,by="cpg")
















