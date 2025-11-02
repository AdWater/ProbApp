# predictFunctions.R
# Predictive and back-transformation functions

#######################################
## back-transformation from Box-Cox transformed predictions

calc_BC_pred_from_eta = function(Qh,A,lambda,eta,Qmin=0.,Qmax=999.,truncType='spike'){



  if (truncType == 'spike'){
    if (lambda>=0){
      Qundef = Qmin
    } else {
      Qundef = Qmax
    }
  } else if (truncType == 'resample') {
    Qundef = NA
  }
  if (lambda==0){
    pred = exp(eta)*(Qh+A)-A
  } else {

    Y = lambda*eta+(Qh+A)^lambda
    pred = Y^(1/lambda)-A

  }
  if (truncType == 'spike'){
    pred[pred<Qmin]=Qmin
    pred[pred>Qmax]=Qmax
  } else if (truncType == 'resample') {
    pred[pred<Qmin]=NA
    pred[pred>Qmax]=NA
  }

  return(pred)
}

#######################################
## back-transformation from WLS predictions

calc_WLS_pred_from_eta = function(Qh,A,eta,Qmin=0.,Qmax=999.,truncType='spike'){
  pred = Qh + eta*(Qh+A)
  if (truncType == 'spike'){
    pred[pred<Qmin]=Qmin
    pred[pred>Qmax]=Qmax
  } else if (truncType == 'resample') {
    pred[pred<Qmin]=NA
    pred[pred>Qmax]=NA
  }
  return(pred)
}

#######################################
## back-transformation from LogSinh predictions

calc_LogSinh_pred_from_eta = function(Qh,A,B,eta,Qmin=0.,Qmax=999.,truncType='spike'){
  Y = eta + calc_LogSinh_tranz(Q=Qh,A=A,B=B)
  pred = calc_inv_LogSinh_tranz(Y=Y,A=A,B=B)
  if (truncType == 'spike'){
    pred[pred<Qmin]=Qmin
    pred[pred>Qmax]=Qmax
  } else if (truncType == 'resample') {
    pred[pred<Qmin]=NA
    pred[pred>Qmax]=NA
  }
  return(pred)
}

#######################################
## Select back-transformation function

calc_pred_from_eta = function(Qh,heteroModel,param,eta,Qmin=0.,Qmax=999.,truncType='spike'){



  if (heteroModel == 'BC'){
    if (is.list(param)){
      pred = calc_BC_pred_from_eta(Qh=Qh,A=param$A,lambda=param$lambda,eta=eta,Qmin=Qmin,Qmax=Qmax,truncType=truncType)
    } else {
      pred = calc_BC_pred_from_eta(Qh=Qh,A=param['A'],lambda=param['lambda'],eta=eta,Qmin=Qmin,Qmax=Qmax,truncType=truncType)
    }
  } else if (heteroModel == 'LogSinh'){
    if (is.list(param)){
      pred = calc_LogSinh_pred_from_eta(Qh=Qh,A=param$A,B=param$B,eta=eta,Qmin=Qmin,Qmax=Qmax,truncType=truncType)
    } else {
      pred = calc_LogSinh_pred_from_eta(Qh=Qh,A=param['A'],B=param['B'],eta=eta,Qmin=Qmin,Qmax=Qmax,truncType=truncType)
    }
  } else if (heteroModel == 'WLS'){
    if (is.list(param)){
      pred = calc_WLS_pred_from_eta(Qh=Qh,A=param$A,eta=eta,Qmin=Qmin,Qmax=Qmax,truncType=truncType)
    } else {
      pred = calc_WLS_pred_from_eta(Qh=Qh,A=param['A'],eta=eta,Qmin=Qmin,Qmax=Qmax,truncType=truncType)
    }
  }



  return(pred)
}

#######################################
## calculate innovations

sim_AR1 = function(nT,mu,sigma,rho){
  if (is.vector(mu)){
    mu_vec = mu
  } else {
    mu_vec = rep(mu,nT)
  }
  if (is.vector(sigma)){
    sigma_vec = sigma
  } else {
    sigma_vec = rep(sigma,nT)
  }
  if (is.vector(rho)){
    rho_vec = rho
  } else {
    rho_vec = rep(rho,nT)
  }

  eta = vector(length = nT)
  eta[1] = mu_vec[1] + rnorm(n=1,mean=0.,sd=sigma_vec[1]/sqrt(1-rho_vec[1]^2))
  flag = 0
  for (t in 2:nT){
    if (is.na(mu_vec[t-1]) || is.na(mu_vec[t])) {
      eta[t] = NA
      flag = 1
    } else {
      if (flag == 1) {
        eta[t] = mu_vec[t] + rnorm(n=1,mean=0.,sd=sigma_vec[t]/sqrt(1-rho_vec[t]^2))
        flag = 0
        next
      }
      eta[t] = rho_vec[t]*(eta[t-1]-mu_vec[t-1]) + mu_vec[t] + rnorm(n=1,mean=0.,sd=sigma_vec[t])
    }

  }

  return(eta)
}

#######################################
## generate replicates

calc_pred_reps = function(Qh,heteroModel,param,nReps=1e2,Qmin=0.,Qmax=999.,truncType='spike',validate=F,strat=NULL){
  
  nT = length(Qh)
  if (is.null(param$mean_eta_0)){
    mean_eta_0 = 0.
  } else {
    mean_eta_0 = param$mean_eta_0
  }
  if (is.null(param$mean_eta_1)){
    mean_eta_1 = 0.
  } else {
    mean_eta_1 = param$mean_eta_1
  }

  rho_eta = param$rho
  sigma_eta = param$sigma_y

  Qh_T = calc_tranz(Q=Qh,heteroModel=heteroModel,param=param) # The transformed simulated streamflow

  mean_eta_0_vec = mean_eta_1_vec = vector(length = nT)
  if (!is.null(strat$index[[strat$type$mean]])){
    for (k in 1:length(strat$index[[strat$type$mean]])){
      keep = strat$index[[strat$type$mean]][[k]]
      mean_eta_0_vec[keep] = mean_eta_0[k]; mean_eta_1_vec[keep] = mean_eta_1[k]
    }    
  } else {
    mean_eta_0_vec = rep(mean_eta_0,nT)
    mean_eta_1_vec = rep(mean_eta_1,nT)
  }


  sigma_eta_vec = vector(length = nT)
  if (!is.null(strat$index[[strat$type$sigma]])){
    for (k in 1:length(strat$index[[strat$type$sigma]])){
      keep = strat$index[[strat$type$sigma]][[k]]
      sigma_eta_vec[keep] = sigma_eta[k]
    }
  } else {
    sigma_eta_vec = rep(sigma_eta,nT)
  }

  rho_eta_vec = vector(length = nT)
  if (!is.null(strat$index[[strat$type$rho]])){
    for (k in 1:length(strat$index[[strat$type$rho]])){
      keep = strat$index[[strat$type$rho]][[k]]
      rho_eta_vec[keep] = rho_eta[k]
    }
  } else {
    rho_eta_vec = rep(rho_eta,nT)
  }

  mean_eta = mean_eta_0_vec+mean_eta_1_vec*Qh_T

  predReps = matrix(nrow=nT,ncol=nReps)

    for (r in 1:nReps){
      eta = sim_AR1(nT,mu=mean_eta,sigma=sigma_eta_vec,rho=rho_eta_vec)
      predReps[,r] = calc_pred_from_eta(Qh=Qh,heteroModel=heteroModel,param=param,
                                        eta=eta,Qmin=Qmin,Qmax=Qmax,truncType=truncType)
    }

    if (nReps==1){predReps=predReps[,1]}
  colnames(predReps)=paste("rep",seq(1:nReps),sep="")

    return(predReps)
}
