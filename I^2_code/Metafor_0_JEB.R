
################################################################################
rm(list = ls())

library(readxl)
library(ape)
library(metafor)
library(ggplot2)
library(gplots)
library(lme4)
library(car)
library(scales)
library(emmeans)
library(arm)
library(QGglmm)
library(MCMCglmm)
library(HDInterval)
library(r2glmm)
library(dispmod)
library(coda)
library(rstan)
library(rethinking)
library(finalfit)
library(dplyr)
library(taxize)
library(rJava)
library(glmulti)
library(MuMIn)



options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


df<-read.table(file="0-analysis_JEB.txt", header= TRUE,sep = "\t")


#### checking and converting 
df$Sel_temp<-as.numeric(df$Sel_temp)
df$Assay_temp<-as.numeric(df$Assay_temp)
df$Rel_assay_temp<-as.numeric(df$Rel_assay_temp)
df$Abs_rel_assay_temp<-as.numeric(df$Abs_rel_assay_temp)

df$Sign<-as.factor(df$Sign)
df$Recombination<-as.factor(df$Recombination)
df$Type_Gen_var<-as.factor(df$Type_Gen_var)
df$Anc_vs_Control<-as.factor(df$Anc_vs_Control)
df$Study_ID<-as.factor(df$Study_ID)
df$Curve_ID<-as.factor(df$Curve_ID)
df$No_of_Gen<-as.numeric(df$No_of_Gen)

df$Sel_coeff<-as.numeric(df$Sel_coeff)
df$Error_new<-as.numeric(df$Error_new)
df$Scientific_name<-as.factor(df$Scientific_name)
df$Scientific_name<-factor(df$Scientific_name)

df$Variance<-(df$Error_new)^2



###### metafor model, used to obtain some values and matrix structure for I² calculations
full_meta <- rma.mv(Sel_coeff, Variance, random = list(~ 1 | Study_ID, ~ 1 | Scientific_name),
                     data = df, method="ML",
                     mods = ~ Type_Gen_var )
summary(full_meta)


##### obtain model matrix
W <- diag(1/full_meta$vi)
X <- model.matrix(full_meta)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

###### the version with frequentist approach in metafor can be found here:
###### https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate


###################################################################################
#I² calculation


#####Posteriors from the best model
post0<-read.table(file="posteriors_0_analysis_JEB.txt", header= TRUE,sep = "\t")


# getting total I²
# first sum the squared posteriors corresponding to the variance components of the random effect structure
##### remember when using medians from Bayesian to ^2
##### s_sigma, t_sigma, b_sigma
##### s_sigma = Study 
##### t_sigma = Species
##### b_sigma = Curve

sum_posts2 <- post0$s_sigma^2 + post0$t_sigma^2


# get I² total
mean(100 * sum_posts2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))))

# get I² for the s_sigma component
mean(100 * post0$s_sigma^2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))))

# get I² for the t_sigma component
mean(100 * post0$t_sigma^2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))))



# get I² CIs total 
quantile(100 * sum_posts2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))),p=c(0.025,0.5,0.975))

# get I² CIs for the s_sigma component
quantile(100 * post0$s_sigma^2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))),p=c(0.025,0.5,0.975))

# get I² CIs for the t_sigma component
quantile(100 * post0$t_sigma^2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))),p=c(0.025,0.5,0.975))


###################################################################################


