library(foreach)
library(doRNG)
library(doFuture)
library(circular)
library(cylcop)

test_sdj = function(theta, phi, K, weight, ...) {
  Rn = In = rep(0, K+1)
  for (i in 1:(K+1)) {
    Rn[i] = (mean(cos(as.vector(theta)*(i-1))) - phi^(i-1))^2
    In[i] = (mean(sin(as.vector(theta)*(i-1))))^2
  }
  ts = length(theta) * sum((Rn + In) * weight(0:K, ...))
  return(ts)
}

p_val = function(theta_x, theta_y, par0 = c(1, 2.5, 0.5, 0.5),
                 lower = c(0, 0, 0.001, 0.001), upper = c(2*pi, 2*pi, 1, Inf),
                 lambda0 = 0.3, K = 10, N = 10000) {
  x = exp(1i*theta_x)
  y = exp(1i*theta_y)
  n = length(x)
  
  logL = function(par) {
    theta_0_ = par[1]; theta_1_ = par[2]; delta_ = par[3]; r_ = par[4]
    -n*log(1-delta_^2) +
      sum(log(1 - 2*delta_*cos(theta_y - theta_0_ - theta_x + 
                                 2*Arg(1 + r_*exp(1i*(theta_x - theta_1_)))) + delta_^2))
  }
  
  estim = try(nlminb(par0, logL, lower = lower, upper = upper))
  theta_0 = estim$par[1]
  theta_1 = estim$par[2]
  delta = estim$par[3]
  r = estim$par[4]
  
  beta_0 = exp(1i*theta_0)
  beta_1 = r*exp(1i*theta_1)
  y_fit = beta_0*(exp(1i*theta_x) + beta_1)/(1 + Conj(beta_1)*exp(1i*theta_x))
  resid = y / y_fit
  
  ts1 = test_sdj((Arg(resid)+2*pi) %% (2*pi), delta, K = K, dpois, lambda = lambda0)
  
  registerDoFuture()
  plan(multisession, workers = 5)
  lista = foreach(j = 1:N) %dorng% {
    ok = FALSE
    while(!ok) {
      eps.0 = circular::rwrappedcauchy(n, mu = circular(0), rho = delta, control.circular = list())
      y.0 = beta_0*(x + beta_1)/(1 + Conj(beta_1)*x) * exp(1i*eps.0)
      theta_y.0 = Arg(y.0)
      
      logL.0 = function(par) {
        theta_0_ = par[1]; theta_1_ = par[2]; delta_ = par[3]; r_ = par[4]
        -n*log(1-delta_^2) +
          sum(log(1 - 2*delta_*cos(theta_y.0 - theta_0_ - theta_x + 
                                     2*Arg(1 + r_*exp(1i*(theta_x - theta_1_)))) + delta_^2))
      }
      
      estim.0 = try(nlminb(par0, logL.0, lower = lower, upper = upper))
      if (('try-error' %in% class(estim.0)) | any(is.na(estim.0$par))) {
        next
      } else {
        ok = TRUE
      }
      
      theta_0.0 = estim.0$par[1]
      theta_1.0 = estim.0$par[2]
      delta.0 = estim.0$par[3]
      r.0 = estim.0$par[4]
      
      beta_0.0 = exp(1i*theta_0.0)
      beta_1.0 = r.0 * exp(1i*theta_1.0)
      y_fit.0 = beta_0.0 * (exp(1i*theta_x) + beta_1.0)/(1 + Conj(beta_1.0)*exp(1i*theta_x))
      resid.0 = y.0 / y_fit.0
      
      ts0 = test_sdj((Arg(resid.0)+2*pi) %% (2*pi), delta.0, K = K, dpois, lambda = lambda0)
    }
    return(ts0)
  }
  
  ts0 = unlist(lista)
  
  return(list(
    p_value = 1 - ecdf(ts0)(ts1),
    test_stat = ts1,
    coefficients = c(
      theta_0 = theta_0,
      theta_1 = theta_1,
      delta = delta,
      r = r
    )))
}
