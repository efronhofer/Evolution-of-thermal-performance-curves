
rm(list=ls())


library(coda)
library(rstan)
library(rethinking)


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

#################  Making datalist for main analyses
datalist <- list(
  Rel_fitness= df$Sel_coeff,
  Error=df$Error_new,
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



#################  Making datalist for analysis with generation time. 
#################  Some points are excluded due to NA in generation time 
df2<-df
df2 <- df[complete.cases(df[ , "No_of_Gen"]), ] 
df2$Study_ID<- factor(df2$Study_ID)
df2$Scientific_name<- factor(df2$Scientific_name)
df2$Curve_ID<- factor(df2$Curve_ID)


datalist_gen <- list(
  Rel_fitness= df2$Sel_coeff,
  Error=df2$Error_new,
  Abs_Temp=df2$Abs_rel_assay_temp,
  Sign=df2$Sign,
  Gen=as.numeric(standardize(df2$No_of_Gen)),
  Recomb=df2$Recombination,
  Var=df2$Type_Gen_var,
  Anc=df2$Anc_vs_Control,
  Taxa=df2$Scientific_name,
  Study=df2$Study_ID,
  Curve=df2$Curve_ID,
  N=nrow(df2)
)


########### Models complete dataset ###########
#################  
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



