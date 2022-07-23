
rm(list=ls())

library(coda)
library(rstan)
library(rethinking)



options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


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


#################  Making datalist for main analyses
datalist <- list(
  Rel_fitness= df$MaxSel_coeff,
  Error=df$MaxError_new,
  Sign=df$Sign_Max,
  Recomb=df$Recombination_Max,
  Var=df$Type_Gen_var_Max,
  Anc=df$Anc_vs_Control_Max,
  Taxa=df$Scientific_name,
  Study=df$Study_ID_Max,
  Curve=df$Curve_ID_Max,
  N=nrow(df)
)




################# making a new datalist by excluding some studies (reduced dataset) for additianal analysis
df2<-df
df2$numero<-1:length(df2$MaxSel_coeff)
df2$No_of_Gen_Max<-as.numeric(df2$No_of_Gen_Max)
df2<-df2[-1,]
df2$Study_ID_Max<- factor(df2$Study_ID_Max)
df2$Curve_ID_Max<- factor(df2$Curve_ID_Max)
df2$Scientific_name<-as.factor(df2$Scientific_name)
df2$Scientific_name<-factor(df2$Scientific_name)

datalist_gen <- list(
  Rel_fitness= df2$MaxSel_coeff,
  Error=df2$MaxError_new,
  Sign=df2$Sign_Max,
  Gen=as.numeric(standardize(df2$No_of_Gen_Max)),
  Recomb=df2$Recombination_Max,
  Var=df2$Type_Gen_var_Max,
  Anc=df2$Anc_vs_Control_Max,
  Taxa=df2$Scientific_name,
  Study=df2$Study_ID_Max,
  Curve=df2$Curve_ID_Max,
  N=nrow(df2)
)


########### Models complete dataset ###########
set.seed(2)
model.1<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  
set.seed(2)
model.2<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-sai[Sign] + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma,
    
    #priors
    sai[Sign]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(2)
model.3<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var] + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(2)
model.4<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-sai[Sign] + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + v[Var],
    
    #priors
    sai[Sign]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    v[Var]~dnorm(0,2),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))

#################  compare models
compare(model.1,model.2,model.3,model.4,
        func=WAIC)







######## Additional analysis with reduced dataset #########
set.seed(2)
model.1.1<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma,
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  
set.seed(2)
model.2.1<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-sai[Sign] + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma,
    
    #priors
    sai[Sign]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(2)
model.3.1<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var] + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma,
    
    #priors
    v[Var]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(2)
model.4.1<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-sai[Sign] + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + v[Var],
    
    #priors
    sai[Sign]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    v[Var]~dnorm(0,2),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(2)
model.5<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-alpha + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + g*(Gen),
    
    #priors
    alpha~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    g ~dnorm(0,1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(2)
model.6<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var]*(Gen) + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma ,
    
    #priors
    v[Var]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(2)
model.7<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var] + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + g*(Gen),
    
    #priors
    v[Var]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    g ~dnorm(0,1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(2)
model.8<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-alpha[Var]*(Gen) + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + v[Var],
    
    #priors
    alpha[Var]~dnorm(0,2),
    v[Var]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  
set.seed(2)
model.9<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-g*(Gen) + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + sai[Sign],
    
    #priors
    sai[Sign]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    g ~dnorm(0,1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  
set.seed(2)
model.10<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var]*(Gen) + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + sai[Sign],
    
    #priors
    v[Var]~dnorm(0,2),
    sai[Sign]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(2)
model.11<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-v[Var] +g*(Gen) + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + sai[Sign],
    
    #priors
    v[Var]~dnorm(0,2),
    sai[Sign]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    g ~dnorm(0,1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(2)
model.12<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <-alpha[Var]*(Gen) + x[Taxa]*t_sigma+ t_bar + z[Study]*s_sigma + sai[Sign] +v[Var],
    
    #priors
    alpha[Var]~dnorm(0,2),
    v[Var]~dnorm(0,2),
    sai[Sign]~dnorm(0,2),
    x[Taxa] ~ dnorm(0, 1),
    z[Study]~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,3),
    t_sigma ~ dexp(1),
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<- z*s_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))




#################  compare final models with reduced dataset
compare(model.1.1,
        model.2.1,model.3.1,model.4.1,model.5,
        model.6,model.7,model.8,model.9,model.10,model.11,model.12,
        func=WAIC)



