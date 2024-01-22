#https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(ggsci)


surviv <- read.csv("~/Dropbox/Projects/kristi_epigenetic_clock/survivourship/GGplot_data.csv")


# Kaplan Meier Survival Curve
km <- with(surviv, Surv(Death))
head(km,80)

km_fit <- survfit(Surv(Death) ~ 1, data=surviv)
summary(km_fit, times = c(1,10,30,50*(1:10)))

#plot(km_fit, xlab="Days", main = 'Kaplan Meyer Plot') #base graphics is always ready
autoplot(km_fit)


km_trt_fit <- survfit(Surv(Death) ~ Treatment, data=surviv)
autoplot(km_trt_fit, ylab = "Survival (%)", xlab = "Chronological age (Days)")+
  theme_classic(base_size = 16) + theme(legend.title = element_blank(),legend.position = c(0.8, 0.8))

ggsave("/Users/emb3/Dropbox/Projects/kristi_epigenetic_clock/manuscript/survivourship.pdf",width = 10.9, height = 4.5, units = "in")



# # Fit Cox Model. Not sure which of the below coc is correct for treating Tube as a nested fixed effect
# cox <- coxph(Surv(Death) ~ Treatment/Tube, data = surviv) #nested
# cox <- coxph(Surv(Death) ~ Treatment, cluster =Tube, data = surviv) #nested
# summary(cox)

# Much better to treat tube as a random effect, with treatment nested in tube
library(coxme)
mefit <- coxme(Surv(Death) ~ Treatment + (1 | Tube/Treatment), data = surviv)#https://cran.r-project.org/web/packages/coxme/vignettes/coxme.pdf
mefit

#For exp(𝐵)>1, the interpretation is even easier, as a value of, say exp(𝐵)=1.259 
#(as is the case for your "stenosis" variable), means that scoring "stenosis" = 1 
#will result in an increased probability (25.9%) of experiencing an end point compared to
#when "stenosis" = 0.     SO 2.4 mean a 140% chance of dying if male?


