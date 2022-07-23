
rm(list=ls())

library(coda)
library(rstan)
library(rethinking)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


df<-read.table(file="2-analysis_JEB.txt", header= TRUE,sep = "\t")


#### checking and converting 
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

df$Exp_Design<-as.factor(df$Exp_Design)
df$Ans_vs_Control<-factor(df$Ans_vs_Control)
df$Scientific_name<-as.factor(df$Scientific_name)
df$Scientific_name<-factor(df$Scientific_name)

#################  Making datalist for main analyses
datalist <- list(
  Rel_fitness= df$Sel_coeff,
  Error=df$Error_new,
  Abs_Temp=df$Abs_rel_assay_temp,
  Side=df$Side,
  Sign=df$Sign,
  Design=as.factor(df$Exp_Design),
  Recomb=df$Recombination,
  Var=df$Type_Gen_var,
  Anc=df$Ans_vs_Control,
  Taxa=df$Scientific_name,
  Study=df$Study_ID,
  Curve=df$Curve_ID,
  N=nrow(df)
)



#################  making a new datalist by excluding some studies (reduced dataset) for additianal analyses
df2 <- df[complete.cases(df[ , "No_of_Gen"]), ] 
df2$Study_ID<-factor(df2$Study_ID)
df2$Curve_ID<-factor(df2$Curve_ID)
df2$Scientific_name<-factor(df2$Scientific_name)


datalist_gen <- list(
  Rel_fitness= df2$Sel_coeff,
  Error=df2$Error_new,
  Abs_Temp=df2$Abs_rel_assay_temp,
  Design=as.factor(df2$Exp_Design),
  Gen=standardize(df2$No_of_Gen),
  Recomb=df2$Recombination,
  Var=df2$Type_Gen_var,
  Anc=df2$Ans_vs_Control,
  Taxa=df2$Scientific_name,
  Study=df2$Study_ID,
  Curve=df2$Curve_ID,
  N=nrow(df2)
)


########### Models complete dataset ###########
#################  
set.seed(1)
model.1<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha + z[Curve]*b_sigma  + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar,
    
    #priors
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(1)
model.2<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  
set.seed(1)
model.3<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha + z[Curve]*b_sigma + g*(Abs_Temp) + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar,
    
    #priors
    alpha ~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
   
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.4<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + z[Curve]*b_sigma + g*(Abs_Temp) + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.5<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] * (Abs_Temp) + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.6<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- g[Var] * (Abs_Temp) + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar + v[Var],
    
    #priors
    g[Var]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################   
set.seed(1)
model.7<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha + z[Curve]*b_sigma  + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar +d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(1)
model.8<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  
set.seed(1)
model.9<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha + z[Curve]*b_sigma + g*(Abs_Temp) + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    alpha ~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.10<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + z[Curve]*b_sigma + g*(Abs_Temp) + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.11<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] * (Abs_Temp) + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.12<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- g[Var] * (Abs_Temp) + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar + v[Var]+d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    g[Var]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  compare models 
compare(model.1,model.2,model.3,model.4,
        model.5,model.6,model.7,model.8,
        model.9,model.10,model.11,model.12,
        func=WAIC)







######## Additional analysis with reduced dataset #########
set.seed(1)
model.1re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha + z[Curve]*b_sigma  + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar,
    
    #priors
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.2re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  
set.seed(1)
model.3re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha + z[Curve]*b_sigma + g*(Abs_Temp) + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar,
    
    #priors
    alpha ~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.4re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + z[Curve]*b_sigma + g*(Abs_Temp) + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.5re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] * (Abs_Temp) + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar,
    
    #priors
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.6re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- g[Var] * (Abs_Temp) + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar + v[Var],
    
    #priors
    g[Var]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(1)
model.7re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha + z[Curve]*b_sigma  + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar +d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(1)
model.8re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  
set.seed(1)
model.9re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha + z[Curve]*b_sigma + g*(Abs_Temp) + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    alpha ~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.10re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + z[Curve]*b_sigma + g*(Abs_Temp) + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.11re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] * (Abs_Temp) + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.12re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- g[Var] * (Abs_Temp) + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar + v[Var]+d[Design],
    
    #priors
    d[Design]~dnorm(0,1),
    g[Var]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))




#################   
set.seed(1)
model.13re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha + z[Curve]*b_sigma  + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar + gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    alpha~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    s_sigma ~ dexp(1),
    sigma~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



#################  
set.seed(1)
model.14re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


#################  
set.seed(1)
model.15re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- alpha + z[Curve]*b_sigma + g*(Abs_Temp) + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    alpha ~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1),
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.16re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] + z[Curve]*b_sigma + g*(Abs_Temp) + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    g ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))


################# 
set.seed(1)
model.17re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- v[Var] * (Abs_Temp) + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))



################# 
set.seed(1)
model.18re<-ulam(
  alist(
    Rel_fitness ~ dnorm(Rel_fitness_true, Error),
    vector[N]:Rel_fitness_true ~ dnorm(mu,sigma),
    mu <- g[Var] * (Abs_Temp) + z[Curve]*b_sigma + k[Study]*s_sigma + x[Taxa]*t_sigma + t_bar + v[Var]+ gen*(Gen),
    
    #priors
    gen~dnorm(0,1),
    g[Var]~dnorm(0,1),
    v[Var]~dnorm(0,1),
    z[Curve] ~ dnorm(0, 1),
    k[Study] ~ dnorm(0, 1),
    x[Taxa] ~ dnorm(0, 1),
    
    t_bar ~ dnorm(0,1),
    b_sigma ~ dexp(1), 
    sigma~ dexp(1),
    s_sigma ~ dexp(1),
    t_sigma ~ dexp(1),
    
    gq> vector[Taxa]:t<<- t_bar + x*t_sigma,
    gq> vector[Study]:s<<-  k*s_sigma,
    gq> vector[Curve]:b<<-  z*b_sigma
    
  ), data=datalist_gen,chains=4, cores=4, iter= 30000,warmup = 1000,log_lik=TRUE, control = list(max_treedepth = 50, adapt_delta = 0.99))




#################  compare final models with reduced dataset
compare(model.1re,model.2re,model.3re,model.4re,model.5re,model.6re,
        model.7re,model.8re,model.9re,model.10re,model.11re,model.12re,
        model.13re,model.14re,model.15re,model.16re,model.17re,model.18re,
        func=WAIC)



