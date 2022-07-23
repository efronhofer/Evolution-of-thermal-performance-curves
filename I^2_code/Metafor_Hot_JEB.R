
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



df<-read.table(file="Hotter_JEB.txt", header= TRUE,sep = "\t")


#### checking and converting 
df$Sign_Max<-as.factor(df$Sign_Max)
df$Recombination_Max<-as.factor(df$Recombination_Max)
df$Type_Gen_var_Max<-as.factor(df$Type_Gen_var_Max)
df$Anc_vs_Control_Max<-as.factor(df$Anc_vs_Control_Max)
df$Study_ID_Max<-as.factor(df$Study_ID_Max)
df$Curve_ID_Max<-as.factor(df$Curve_ID_Max)

df$MaxSel_coeff<-as.numeric(df$MaxSel_coeff)
df$MaxError_new<-as.numeric(df$MaxError_new)
df$Sign_Max<-factor(df$Sign_Max)
df$Scientific_name<-as.factor(df$Scientific_name)
df$Scientific_name<-factor(df$Scientific_name)


df$Variance<-(df$MaxError_new)^2



full_meta <- rma.mv(MaxSel_coeff, Variance, random = list(~ 1 | Study_ID_Max, ~ 1 | Scientific_name),
                    data = df, method="ML")
summary(full_meta)



W <- diag(1/full_meta$vi)
X <- model.matrix(full_meta)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
###### the version with frequentist approach in metafor can be found here:
###### https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate



###################################################################################
#I² calculation 


#####Posteriors from the best model
post0<-read.table(file="posteriors_hotter_JEB.txt", header= TRUE,sep = "\t")


# getting total I²
# first sum the squared posteriors corresponding to the variance components of the random effect structure
##### remember when using medians from Bayesian to ^2
##### s_sigma, t_sigma, b_sigma
##### s_sigma = Study 
##### t_sigma = Species
sum_posts2 <- post0$s_sigma^2 + post0$t_sigma^2



# get I² total
mean(100 * sum_posts2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))))

# get I² for the s_sigma component
mean(100 * post0$s_sigma^2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))))

# get I² for the t_sigma component
mean(100 * post0$t_sigma^2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))))



# get I² CIs total total
quantile(100 * sum_posts2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))),p=c(0.025,0.5,0.975))

# get I² CIs for the s_sigma component
quantile(100 * post0$s_sigma^2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))),p=c(0.025,0.5,0.975))

# get I² CIs for the t_sigma component
quantile(100 * post0$t_sigma^2 / (sum_posts2 + (full_meta$k-full_meta$p)/sum(diag(P))),p=c(0.025,0.5,0.975))


###################################################################################


