# residualFunctions.R
# Calculates residuals

#######################################
## Calc transformed residual (eta)

calc_eta = function(Qobs,Qh,param,heteroModel){
  eta = calc_tranz(Q=Qobs,Qh=Qh,param=param,heteroModel=heteroModel) - calc_tranz(Q=Qh,Qh=Qh,param=param,heteroModel=heteroModel)

  return(eta)
}

#######################################
## Calc standardised residual

calc_std_resids = function(data,param,heteroModel,opt,strat=NULL){
  Qobs = data[[opt$obs]]
  Qh = data[[opt$pred]]
  Qh_T = calc_tranz(Q=Qh,heteroModel=heteroModel,param=param) # The transformed simulated streamflow
  eta = calc_eta(Qobs,Qh,param,heteroModel)

  nT = length(eta)

  if (is.null(strat)){strat=set_strat_all(nT)}
  
  #n = length(eta)-sum(is.na(eta))

  mu0 = param$mean_eta_0 # normalising
  mu1 = param$mean_eta_1
  sigma = param$sigma
  rho = param$rho

  mu0_vec = mu1_vec = vector(length = nT)
  for (k in 1:length(strat$index[[strat$type$mean]])){
    keep = strat$index[[strat$type$mean]][[k]]
    mu0_vec[keep] = mu0[k]; mu1_vec[keep] = mu1[k]
  }

  sigma_vec = vector(length = nT)
  for (k in 1:length(strat$index[[strat$type$sigma]])){
    keep = strat$index[[strat$type$sigma]][[k]]
    sigma_vec[keep] = sigma[k]
  }

  rho_vec = vector(length = nT)
  for (k in 1:length(strat$index[[strat$type$rho]])){
    keep = strat$index[[strat$type$rho]][[k]]
    rho_vec[keep] = rho[k]
  }

  mu = mu0_vec+(Qh_T*mu1_vec)
  eta.star = eta-mu
  sigmaEta = sigma_vec/(sqrt(1-rho_vec^2))
  nu = (eta.star)/sigmaEta

  return(nu)
}

#######################################
## Calc innovations

calc_innovations = function(data,param,heteroModel,opt){
  Qobs = data[[opt$obs]]
  Qh = data[[opt$pred]]
  n=length(Qobs)
  y = vector(length=n)
  eta = vector(length=n)
  y[1] = eta[1]
  eta = calc_eta(Qobs,Qh,param,heteroModel)
  mu0 = param$mean_eta_0
  mu1 = param$mean_eta_1
  mu = mu0+Qh*mu1
  phi = param$rho

  for (i in 2:n){
    y[i] = (eta[i]-mu[i]) - (phi*(eta[i-1]-mu[i]))
  }

  return(y)
}
