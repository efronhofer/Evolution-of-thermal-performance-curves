
rm(list=ls())

library(coda)
library(rstan)
library(rethinking)



options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#################  
#################  Funnel plot 0-analysis 
#################  
df<-read.table(file="0-analysis_JEB.txt", header= TRUE,sep = "\t")


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



#################  Making datalist for main analyses
datalist_funnel <- list(
  Rel_fitness= df$Sel_coeff,
  Error=1/df$Error_new,
  Abs_Temp=df$Abs_rel_assay_temp,
  Sign=df$Sign,
  Recomb=df$Recombination,
  Var=df$Type_Gen_var,
  Anc=df$Anc_vs_Control,
  Taxa=df$Scientific_name,
  Study=df$Study_ID,
  Curve=df$Curve_ID,
  N=nrow(df)
)


#################  General intercept + general slope
set.seed(2)
model.funnel.0intercept<-ulam(
  alist(
    Error ~ dnorm(mu, sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + g*(Rel_fitness),
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_funnel,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  General intercept
set.seed(2)
model.funnel.0slope<-ulam(
  alist(
    Error ~ dnorm(mu, sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma ,
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_funnel,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  compare models
compare(model.funnel.0intercept,model.funnel.0slope, func=WAIC)





#################  
#################  Funnel plot all studies
#################  
df_tot<-read.table(file="All_data_JEB.txt", header= TRUE,sep = "\t")

df_tot$Sel_temp<-as.numeric(df_tot$Sel_temp)
df_tot$Assay_temp<-as.numeric(df_tot$Assay_temp)
df_tot$Rel_assay_temp<-as.numeric(df_tot$Rel_assay_temp)
df_tot$Abs_rel_assay_temp<-as.numeric(df_tot$Abs_rel_assay_temp)

df_tot$Sign<-as.factor(df_tot$Sign)
df_tot$Recombination<-as.factor(df_tot$Recombination)
df_tot$Type_Gen_var<-as.factor(df_tot$Type_Gen_var)
df_tot$Anc_vs_Control<-as.factor(df_tot$Anc_vs_Control)
df_tot$Study_ID<-as.factor(df_tot$Study_ID)
df_tot$Curve_ID<-as.factor(df_tot$Curve_ID)

df_tot$Sel_coeff<-as.numeric(df_tot$Sel_coeff)
df_tot$Error_new<-as.numeric(df_tot$Error_new)
df_tot$Scientific_name<-as.factor(df_tot$Scientific_name)
df_tot$Scientific_name<-factor(df_tot$Scientific_name)


#################  Making datalist for main analyses
datalist_funnel_tot <- list(
  Rel_fitness= standardize(df_tot$Sel_coeff),
  Error=1/df_tot$Error_new,
  Abs_Temp=df_tot$Abs_rel_assay_temp,
  Sign=df_tot$Sign,
  Recomb=df_tot$Recombination,
  Var=df_tot$Type_Gen_var,
  Anc=df_tot$Anc_vs_Control,
  Taxa=df_tot$Scientific_name,
  Study=df_tot$Study_ID,
  Curve=df_tot$Curve_ID,
  N=nrow(df_tot)
)


#################  General intercept + general slope
set.seed(2)
model.funnel.total_slope<-ulam(
  alist(
    Error ~ dnorm(mu, sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Curve]*b_sigma + g*(Rel_fitness) + k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Curve]~ dnorm(0, 1),
    k[Study]~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma ~ dexp(1),
    b_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<- z*b_sigma,
    gq> vector[Study]:s<<- k*s_sigma
    
  ), data=datalist_funnel_tot,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.9999))


#################  General intercept 
set.seed(2)
model.funnel.total_intercept<-ulam(
  alist(
    Error ~ dnorm(mu, sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Curve]*b_sigma + k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Curve]~ dnorm(0, 1),
    k[Study]~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma ~ dexp(1),
    b_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<- z*b_sigma,
    gq> vector[Study]:s<<- k*s_sigma
    
  ), data=datalist_funnel_tot,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.9999))


#################  compare models
compare(model.funnel.total_intercept,model.funnel.total_slope, func=WAIC)




#################  
#################  Funnel plot multi analysis
#################  
df<-read.table(file="Multi_analysis_JEB.txt", header= TRUE,sep = "\t")


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

df$Sel_coeff<-as.numeric(df$Sel_coeff)
df$Error_new<-as.numeric(df$Error_new)
df$Scientific_name<-as.factor(df$Scientific_name)
df$Scientific_name<-factor(df$Scientific_name)


###### Making datalist for multi analysis 
df_tot_3<-df

datalist_funnel_tot3 <- list(
  Rel_fitness= standardize(df_tot_3$Sel_coeff), #to avoid diverteng trasnsitions
  Error=(1/df_tot_3$Error_new),
  Abs_Temp=df_tot_3$Abs_rel_assay_temp,
  Sign=df_tot_3$Sign,
  Recomb=df_tot_3$Recombination,
  Var=df_tot_3$Type_Gen_var,
  Anc=df_tot_3$Anc_vs_Control,
  Taxa=df_tot_3$Scientific_name,
  Study=df_tot_3$Study_ID,
  Curve=df_tot_3$Curve_ID,
  N=nrow(df_tot_3)
)


#################  General intercept + general slope
set.seed(2)
model.funnel.Multi_slope<-ulam(
  alist(
    Error ~ dnorm(mu, sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Curve]*b_sigma + g*(Rel_fitness) + k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Curve]~ dnorm(0, 1),
    k[Study]~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma ~ dexp(1),
    b_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<- z*b_sigma,
    gq> vector[Study]:s<<- k*s_sigma
    
  ), data=datalist_funnel_tot3,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.9999))


#################  General intercept 
set.seed(2)
model.funnel.Multi_intercept<-ulam(
  alist(
    Error ~ dnorm(mu, sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Curve]*b_sigma + k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Curve]~ dnorm(0, 1),
    k[Study]~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma ~ dexp(1),
    b_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<- z*b_sigma,
    gq> vector[Study]:s<<- k*s_sigma
    
  ), data=datalist_funnel_tot3,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.9999))

#################  compare models
compare(model.funnel.Multi_intercept,model.funnel.Multi_slope, func = WAIC)






#################  
#################  Funnel plot 2-analysis
#################  
df<-read.table(file="2-analysis_JEB.txt", header= TRUE,sep = "\t")

df$Sel_temp<-as.numeric(df$Sel_temp)
df$Assay_temp<-as.numeric(df$Assay_temp)
df$Rel_assay_temp<-as.numeric(df$Rel_assay_temp)
df$Abs_rel_assay_temp<-as.numeric(df$Abs_rel_assay_temp)

df$Design<-as.factor(df$Exp_Design)
df$Side<-as.factor(df$Side)
df$Sign<-as.factor(df$Sign)
df$Recombination<-as.factor(df$Recombination)
df$Type_Gen_var<-as.factor(df$Type_Gen_var)
df$Ans_vs_Control<-as.factor(df$Ans_vs_Control)
df$Study_ID<-as.factor(df$Study_ID)
df$Curve_ID<-as.factor(df$Curve_ID)
df$No_of_Gen<-as.numeric(df$No_of_Gen)

df$Sel_coeff<-as.numeric(df$Sel_coeff)
df$Error_new<-as.numeric(df$Error_new)
df$Scientific_name<-as.factor(df$Scientific_name)
df$Scientific_name<-factor(df$Scientific_name)


###### Making datalist for funnel plot 2-analysis 
df_tot_3<-df

datalist_funnel_tot2 <- list(
  Rel_fitness= standardize(df_tot_3$Sel_coeff), #to avoid diverteng trasnsitions
  Error=(1/df_tot_3$Error_new),
  Abs_Temp=df_tot_3$Abs_rel_assay_temp,
  Sign=df_tot_3$Sign,
  Recomb=df_tot_3$Recombination,
  Var=df_tot_3$Type_Gen_var,
  Anc=df_tot_3$Anc_vs_Control,
  Taxa=df_tot_3$Scientific_name,
  Study=df_tot_3$Study_ID,
  Curve=df_tot_3$Curve_ID,
  N=nrow(df_tot_3)
)


#################  General intercept + general slope
set.seed(2)
model.funnel.2slope<-ulam(
  alist(
    Error ~ dnorm(mu, sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Curve]*b_sigma + g*(Rel_fitness) + k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Curve]~ dnorm(0, 1),
    k[Study]~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma ~ dexp(1),
    b_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<- z*b_sigma,
    gq> vector[Study]:s<<- k*s_sigma
    
  ), data=datalist_funnel_tot2,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 150, adapt_delta = 0.999999999999999,stepsize = 0.01))


#################  General intercept 
set.seed(2)
model.funnel.2intercept<-ulam(
  alist(
    Error ~ dnorm(mu, sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Curve]*b_sigma + k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Curve]~ dnorm(0, 1),
    k[Study]~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma ~ dexp(1),
    b_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<- z*b_sigma,
    gq> vector[Study]:s<<- k*s_sigma
    
  ), data=datalist_funnel_tot2,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 150, adapt_delta = 0.999999999999999,stepsize = 0.01))


#################  compare models
compare(model.funnel.2slope,model.funnel.2intercept, func = WAIC)






#################  
#################  Funnel plot 0-analysis 
#################  
df<-read.table(file="Hotter_JEB.txt", header= TRUE,sep = "\t")


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

str(df)

#################  Making datalist for main analyses
datalist <- list(
  Rel_fitness= (df$MaxSel_coeff),
  Error=(1/df$MaxError_new),
  Taxa=df$Scientific_name,
  Study=df$Study_ID_Max,
  N=nrow(df)
)



#################  General intercept + general slope
set.seed(2)
model.funnel.HOTslope<-ulam(
  alist(
    Error ~ dnorm(mu, sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + g*(Rel_fitness),
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  General intercept
set.seed(2)
model.funnel.HOTintercept<-ulam(
  alist(
    Error ~ dnorm(mu, sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma ,
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  compare models
compare(model.funnel.HOTintercept,model.funnel.HOTslope, func=WAIC)




