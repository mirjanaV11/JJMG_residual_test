##########################################################################
## datasets come from the German weather service at nearby stations 
## Freudenstadt (fr) / Hornisgrinde (ho): distance about 22 km
##
## datasets contain wind directions measured at Wednesday 6:00 and noon,  
## from 2015-01-07 to 2023-12-27, some missing values
## dataset data.bf: n=463
##########################################################################
library(circular); library(foreach); library(doRNG); library(doFuture); 
library(cylcop); library(xtable)

######################################################################
## function definition
######################################################################

# Tn test statistic
test_sdj = function(theta, phi, K, weight, ...) {
  Rn = In = rep(0, K+1)
  for (i in 1:(K+1)) {
    Rn[i] = (mean(cos(as.vector(theta)*(i-1))) - phi^(i-1))^2
    In[i] = (mean(sin(as.vector(theta)*(i-1))))^2
  }
  ts = length(theta)*sum((Rn + In)*weight(0:K, ...))
  return(ts)
}

# Wn test statistic
test_stat_w = function(x){
  n = length(x)
  index = 1:n
  return(1/(12*n)+sum((x-(2*index-1)/(2*n)-(mean(x)-1/2))^2))
  
}

# Kn test statistic
test_stat_k = function(x){
  n = length(x)
  index = 1:n
  return(max(x - (index-1)/n)+max(index/n-x))
}

# log likelihood
logLm = function(par, theta_x, theta_y) {
  theta_0_pom = par[1];  theta_1_pom = par[2]; phi_pom = par[3]; r_pom = par[4] 
  -length(theta_x)*log(1 - phi_pom^2) +
    sum(log( 1 - 2*phi_pom*cos(theta_y - theta_0_pom - theta_x + 2*Arg(1 + r_pom*exp(1i*(theta_x - theta_1_pom)))) 
             + phi_pom^2 ))
}


# perform tests
test.result = function(theta_x, theta_y, B=10000, par.start) { #theta_x, theta_y in rad
  x = exp(1i*theta_x)
  y = exp(1i*theta_y)
  res = try(nlminb( par.start, logLm, lower=c(0.001, 0, 0.001, 0.001), 
                    upper=c(2*pi, 2*pi, 1, Inf), theta_x=theta_x, theta_y=theta_y))
  theta_0 = res$par[1]; theta_1 = res$par[2]; phi = res$par[3]; r = res$par[4]
  beta_0 = exp(1i*theta_0)
  beta_1 = r*exp(1i*theta_1)
  y_fit = beta_0*(exp(1i*theta_x) + beta_1)/(1 + Conj(beta_1)*exp(1i*theta_x))
  reziduali = y/y_fit
  
  ts1_03 = test_sdj((Arg(reziduali)+2*pi)%%(2*pi), phi, K = 10, dpois, lambda = 0.3)
  ts1_05 = test_sdj((Arg(reziduali)+2*pi)%%(2*pi), phi, K = 10, dpois, lambda = 0.5)
  ts1_1 = test_sdj((Arg(reziduali)+2*pi)%%(2*pi), phi, K = 10, dpois, lambda = 1)
  sorted.data = pwrappedcauchy(as.vector(sort((Arg(reziduali)+2*pi)%%(2*pi))), scale = -log(phi))
  ts1.k = test_stat_k(sorted.data)
  ts1.w = test_stat_w(sorted.data)    
  
  n = length(theta_x)
  Nponavljanja = B
  registerDoFuture()
  plan(multisession, workers=5)
  lista = foreach (j = 1:Nponavljanja) %dorng% {
    ok = FALSE
    while(!ok) {
      eps.0 = circular::rwrappedcauchy(n, mu = circular(0), rho = phi, control.circular=list())
      y.0 = beta_0*(x + beta_1)/(1 + Conj(beta_1)*x)*exp(1i*eps.0)
      theta_y.0 = Arg(y.0) 
      
      logL.0 = function(par) {
        theta_0_pom = par[1];  theta_1_pom = par[2]; phi_pom = par[3]; r_pom = par[4] 
        -n*log(1-phi_pom^2) +
          sum(log(1 - 2*phi_pom*cos(theta_y.0 - theta_0_pom - theta_x + 2*Arg(1 + r_pom*exp(1i*(theta_x - theta_1_pom)))) + phi_pom^2))
      }
      
      ocena.0 = try(nlminb( par.start, logL.0, lower = c(0.001, 0, 0.001, 0.001), upper = c(2*pi, 2*pi, 1, Inf)))
      if (('try-error' %in% class(ocena.0)) | any(is.na(ocena.0$par) == TRUE)) {next}
      else {ok = TRUE}
      theta_0.0 = ocena.0$par[1]
      theta_1.0 = ocena.0$par[2]
      phi.0 = ocena.0$par[3]
      r.0 = ocena.0$par[4]
      
      beta_0.0 = exp(1i*theta_0.0)
      beta_1.0 = r.0*exp(1i*theta_1.0)
      y_fit.0 = beta_0.0*(exp(1i*theta_x) + beta_1.0)/(1 + Conj(beta_1.0)*exp(1i*theta_x))
      reziduali.0 = y.0/y_fit.0
      
      ts0_03 = test_sdj((Arg(reziduali.0)+2*pi)%%(2*pi), phi.0, K = 10, dpois, lambda = 0.3)
      ts0_05 = test_sdj((Arg(reziduali.0)+2*pi)%%(2*pi), phi.0, K = 10, dpois, lambda = 0.5)
      ts0_1 = test_sdj((Arg(reziduali.0)+2*pi)%%(2*pi), phi.0, K = 10, dpois, lambda = 1)
      sorted.data.0 = pwrappedcauchy(as.vector(sort((Arg(reziduali.0)+2*pi) %% (2*pi))), scale = -log(phi.0))
      ts0.k = test_stat_k(sorted.data.0)
      ts0.w = test_stat_w(sorted.data.0)    
    }
    return(list(c(ts0_03, ts0_05, ts0_1, ts0.k, ts0.w)))
  }
  
  ts0_03 = ts0_05 = ts0_1 = ts0.k = ts0.w = rep(0, Nponavljanja)
  for(i in 1:length(lista)) {
    ts0_03[i] = lista[[i]][[1]][1]
    ts0_05[i] = lista[[i]][[1]][2]
    ts0_1[i] = lista[[i]][[1]][3]
    ts0.k[i] = lista[[i]][[1]][4]
    ts0.w[i] = lista[[i]][[1]][5]
  }
  
  # p-values
  p.03 = 1 - ecdf(ts0_03)(ts1_03)
  p.05 = 1 - ecdf(ts0_05)(ts1_05)
  p.1 = 1 - ecdf(ts0_1)(ts1_1)
  p.k = 1 - ecdf(ts0.k)(ts1.k)
  p.w = 1 - ecdf(ts0.w)(ts1.w) 
  return( c(p.03, p.05, p.1, p.k, p.w) )
}

######################################################################
## end function definition
######################################################################

## Black forest ##
#save(data.bf, file="data.bf.Rda")
load("data.bf.Rda")
attach(data.bf)

## regression models ##
## ho12 ~ fr12
theta_x = frwe12; theta_y = howe12
x = exp(1i*theta_x); y = exp(1i*theta_y)
par = c(1,3,0.5,0.5)
res = try(nlminb( par, logLm, lower=c(0.001, 0, 0.001, 0.001), 
                     upper=c(2*pi, 2*pi, 1, Inf), theta_x=theta_x, theta_y=theta_y))
list(res$par, res$objective)
theta_0 = res$par[1]; theta_1 = res$par[2]; phi = res$par[3]; r = res$par[4]
beta_0 = exp(1i*theta_0); beta_1 = r*exp(1i*theta_1)
## conditional distribution (Kato, p. 637-638)
# parameter phi is the concentration or precision parameter. 
# If phi=1, then covariates and responses are correlated without error. 
# The smaller the value of p, the less concentrated the error variables. 
# When phi=0, the variable e has a uniform distribution on the circle
# Arg(E[y|x]) = m_yx, |E[y|x]| = phi: mean resultant length, 1-phi: variance
xx = seq(0, 2*pi-0.1, pi/4)
m_yx = Arg( beta_0*(exp(1i*xx) + beta_1)/(1 + Conj(beta_1)*exp(1i*xx)) ) %% (2*pi)
rbind(xx, m_yx, m_yx-xx)
y_fit = beta_0*(exp(1i*theta_x) + beta_1)/(1 + Conj(beta_1)*exp(1i*theta_x))
hofr12.res = y/y_fit
c(1-phi, 1-rho.circular(Arg(hofr12.res)), 1-rho.circular(theta_y)) # same as var(theta_y)
res.tab = matrix(nrow=4, ncol=6)
res.tab[1,] = round( c(theta_0,theta_1, r, phi, m_yx[c(2,4)]), 2)
###############################################################################

## ho06 ~ fr06
theta_x = frwe06; theta_y = howe06
x = exp(1i*theta_x); y = exp(1i*theta_y)
par = c(1,3,0.5,0.5)
res = try(nlminb( par, logLm, lower=c(0.001, 0, 0.001, 0.001), 
                  upper=c(2*pi, 2*pi, 1, Inf), theta_x=theta_x, theta_y=theta_y))
list(res$par, res$objective)
theta_0 = res$par[1]; theta_1 = res$par[2]; phi = res$par[3]; r = res$par[4]
beta_0 = exp(1i*theta_0); beta_1 = r*exp(1i*theta_1)
## conditional distribution
xx = seq(0, 2*pi-0.1, pi/4)
m_yx = Arg( beta_0*(exp(1i*xx) + beta_1)/(1 + Conj(beta_1)*exp(1i*xx)) ) %% (2*pi)
rbind(xx, m_yx, m_yx-xx)
#
y_fit = beta_0*(exp(1i*theta_x) + beta_1)/(1 + Conj(beta_1)*exp(1i*theta_x))
hofr06.res = y/y_fit
c(1-phi, 1-rho.circular(Arg(hofr06.res)), 1-rho.circular(theta_y))
# regression model reduces variance by factor 1.35
res.tab[2,] = round( c(theta_0,theta_1, r, phi, m_yx[c(2,4)]), 2)
###############################################################################

## ho12 ~ ho06
theta_x = howe06; theta_y = howe12
x = exp(1i*theta_x); y = exp(1i*theta_y)
par = c(5,5,0.5,0.5)
res = try(nlminb( par, logLm, lower=c(0.01, 0, 0.001, 0.001), 
                  upper=c(2*pi, 2*pi, 1, Inf), theta_x=theta_x, theta_y=theta_y))
list(res$par, res$objective)
theta_0 = res$par[1]; theta_1 = res$par[2]; phi = res$par[3]; r = res$par[4]
beta_0 = exp(1i*theta_0); beta_1 = r*exp(1i*theta_1)
## conditional distribution
xx = seq(0, 2*pi-0.1, pi/4)
m_yx = Arg( beta_0*(exp(1i*xx) + beta_1)/(1 + Conj(beta_1)*exp(1i*xx)) ) %% (2*pi)
rbind(xx, m_yx, m_yx-xx)
#
y_fit = beta_0*(exp(1i*theta_x) + beta_1)/(1 + Conj(beta_1)*exp(1i*theta_x))
ho1206.res = y/y_fit
c(1-phi, 1-rho.circular(Arg(ho1206.res)), 1-rho.circular(theta_y))
# regression model reduces variance by factor 2.4/2.1
res.tab[3,] = round( c(theta_0,theta_1, r, phi, m_yx[c(2,4)]), 2)
###############################################################################

## fr12 ~ fr06
theta_x = frwe06; theta_y = frwe12
x = exp(1i*theta_x); y = exp(1i*theta_y)
par = c(1,3,0.5,0.5)
res = try(nlminb( par, logLm, lower=c(0.001, 0, 0.001, 0.001), 
                  upper=c(2*pi, 2*pi, 1, Inf), theta_x=theta_x, theta_y=theta_y))
list(res$par, res$objective)
theta_0 = res$par[1]; theta_1 = res$par[2]; phi = res$par[3]; r = res$par[4]
beta_0 = exp(1i*theta_0); beta_1 = r*exp(1i*theta_1)
## conditional distribution
xx = seq(0, 2*pi-0.1, pi/4)
m_yx = Arg( beta_0*(exp(1i*xx) + beta_1)/(1 + Conj(beta_1)*exp(1i*xx)) ) %% (2*pi)
rbind(xx, m_yx, m_yx-xx)
c(1-phi, 1-rho.circular(theta_y)) # same as var(theta_y)
#
y_fit = beta_0*(exp(1i*theta_x) + beta_1)/(1 + Conj(beta_1)*exp(1i*theta_x))
fr1206.res = y/y_fit
c(1-phi, 1-rho.circular(Arg(fr1206.res)), 1-rho.circular(theta_y))
# regression model reduces variance by factor 1.8
res.tab[4,] = round( c(theta_0,theta_1, r, phi, m_yx[c(2,4)]), 2)
xtable(res.tab)
###############################################################################

par(mfrow=c(2,2), mar=c(1,1,2,1))
plot(Arg(hofr12.res), stack=TRUE, shrink=1., main="ho12 ~ fr12")
plot(Arg(hofr06.res), stack=TRUE, shrink=1., main="ho06 ~ fr06")
plot(Arg(ho1206.res), stack=TRUE, shrink=1., main="ho12 ~ ho06")
plot(Arg(fr1206.res), stack=TRUE, shrink=1., main="fr12 ~ fr06")
###############################################################################

## tests
# full datasets
pval.tab = matrix(nrow=12, ncol=5)
par = c(1,3,0.5,0.5)
B = 10000
pval.tab[1,] = test.result(frwe12, howe12, B, par)
pval.tab[2,] = test.result(frwe06, howe06, B, par)
pval.tab[3,] = test.result(howe06, howe12, B, c(5,5,0.5,0.5))
pval.tab[4,] = test.result(frwe06, frwe12, B, par)
# 
n = length(frwe12); n1 = n-199
pval.tab[5,] = test.result(frwe12[n1:n], howe12[n1:n], B, par)
pval.tab[6,] = test.result(frwe06[n1:n], howe06[n1:n], B, par)
pval.tab[7,] = test.result(howe06[n1:n], howe12[n1:n], B, c(5,5,0.5,0.5))
pval.tab[8,] = test.result(frwe06[n1:n], frwe12[n1:n], B, par)

n = length(frwe12); n1 = n-99
pval.tab[9,] = test.result(frwe12[n1:n], howe12[n1:n], B, par)
pval.tab[10,] = test.result(frwe06[n1:n], howe06[n1:n], B, par)
pval.tab[11,] = test.result(howe06[n1:n], howe12[n1:n], B, c(5,5,0.5,0.5))
pval.tab[12,] = test.result(frwe06[n1:n], frwe12[n1:n], B, par)
#
xtable(pval.tab)
detach(data.bf)
###############################################################################
