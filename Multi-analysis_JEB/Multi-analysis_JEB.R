
rm(list=ls())

library(coda)
library(rstan)
library(rethinking)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


df<-read.table(file="Multi_analysis_JEB.txt", header= TRUE,sep = "\t")


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

df$Sel_coeff<-as.numeric(df$Sel_coeff)
df$Error_new<-as.numeric(df$Error_new)
df$Scientific_name<-as.factor(df$Scientific_name)
df$Scientific_name<-factor(df$Scientific_name)


###### Making datalist for main analyses 
datalist <- list(
  Rel_fitness= df$Sel_coeff,
  Error=df$Error_new,
  Abs_Temp=df$Abs_rel_assay_temp,
  Rel_Temp=standardize(df$Rel_assay_temp),
  Rel_Temp_s2=standardize(df$Rel_assay_temp)^2,
  Rel_Temp_s3=standardize(df$Rel_assay_temp)^3,
  Sign=df$Sign,
  Recomb=df$Recombination,
  Var=df$Type_Gen_var,
  Anc=df$Anc_vs_Control,
  Taxa=df$Scientific_name,
  Study=df$Study_ID,
  Curve=df$Curve_ID,
  N=nrow(df)
)


#making a new datalist by excluding some studies (reduced dataset) for additional analysis
df2 <- df[complete.cases(df[ , "No_of_Gen"]), ] 
df2$Study_ID<-factor(df2$Study_ID)
df2$Curve_ID<-factor(df2$Curve_ID)
df2$Scientific_name<-as.factor(df2$Scientific_name)
df2$Scientific_name<-factor(df2$Scientific_name)

datalist_gen <- list(
  Rel_fitness= df2$Sel_coeff,
  Error=df2$Error_new,
  Abs_Temp=df2$Abs_rel_assay_temp,
  Rel_Temp=standardize(df2$Rel_assay_temp),
  Rel_Temp_s2=standardize(df2$Rel_assay_temp)^2,
  Rel_Temp_s3=standardize(df2$Rel_assay_temp)^3,
  Gen=standardize(df2$No_of_Gen),
  Sign=df2$Sign,
  Recomb=df2$Recombination,
  Var=df2$Type_Gen_var,
  Anc=df2$Anc_vs_Control,
  Taxa=df2$Scientific_name,
  Study=df2$Study_ID,
  Curve=df2$Curve_ID,
  N=nrow(df2)
)




########### Models complete dataset ###########
set.seed(1)
model.1<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-alpha + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(1)
model.2<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))




################# 
set.seed(1)
model.3<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.4<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.5<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var] + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  
set.seed(1)
model.6<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var]+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.7<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.8<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.9<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-sai[Sign] + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma,
    
    #priors
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.10<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- sai[Sign]+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.11<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.12<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.13<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var] + sai[Sign] + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.14<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + sai[Sign]+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.15<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.16<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  model comparison 
compare(model.1,model.2,model.3,model.4,
        model.5,model.6,model.7,model.8,model.9,model.10,model.11,model.12,
        model.13,model.14,model.15,model.16,
        func=WAIC)










######## Additional analysis with reduced dataset #########
set.seed(1)
model.1re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-alpha + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################
set.seed(1)
model.2re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))





#################
set.seed(1)
model.3re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################
set.seed(1)
model.4re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################
set.seed(1)
model.5re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var] + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################
set.seed(1)
model.6re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var]+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################
set.seed(1)
model.7re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################
set.seed(1)
model.8re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))




#################
set.seed(1)
model.9re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-sai[Sign] + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma,
    
    #priors
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################
set.seed(1)
model.10re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- sai[Sign]+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################
set.seed(1)
model.11re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################
set.seed(1)
model.12re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################
set.seed(1)
model.13re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var] + sai[Sign] + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################set.seed(1)
model.14re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + sai[Sign]+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################
set.seed(1)
model.15re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################
set.seed(1)
model.16re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))










######## Additional analysis with reduced dataset #########
set.seed(1)
model.17re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-alpha + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma + gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.18re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma + gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.19re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma + gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.20re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma + gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.21re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var] + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.22re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var]+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.23re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.24re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.25re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-sai[Sign] + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.26re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- sai[Sign]+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.27re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.28re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.29re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var] + sai[Sign] + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar + k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.30re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + sai[Sign]+ g*(Rel_Temp) + z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.31re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.32re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + sai[Sign]+ g*(Rel_Temp) + g2*(Rel_Temp_s2)+ g3*(Rel_Temp_s3)+ z[Curve]*b_sigma + x[Taxa]*t_sigma + t_bar+ k[Study]*s_sigma+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    sai[Sign]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    g2 ~ dnorm(0, 1),
    g3 ~ dnorm(0, 1),
    x[Taxa]~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    t_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma,
    gq> vector[Study]:s<<-  k*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))




#################  compare models with reduced dataset
compare(model.1,model.2,model.3,model.4,
        model.5,model.6,model.7,model.8,model.9,model.10,model.11,model.12,
        model.13,model.14,model.15,model.16,
        model.1re,model.2re,model.3re,model.4re,
        model.5re,model.6re,model.7re,model.8re,model.9re,model.10re,model.11re,model.12re,
        model.13re,model.14re,model.15re,model.16re,
        model.17re,model.18re,model.19re,model.20re,
        model.21re,model.22re,model.23re,model.24re,
        model.25re,model.26re,model.27re,model.28re,
        model.29re,model.30re,model.31re,model.32re,
        func=WAIC)



