# calibrateFunctions.R
# Calibrates the model

#######################################

cal_mean = function(eta,Qh_T,meantype){

  if (meantype=="linear"){
    m = lm(eta~Qh_T,na.action=na.omit)
    mu0 = m$coefficients[1]
    mu1 = m$coefficients[2]
  } else if (meantype=="constant"){
    mu0 = mean(eta,na.rm=T)
    mu1 = 0.
  } else if (meantype=="zero"){
    mu0 = 0.
    mu1 = 0.
  } else{
    mu0 = 0.
    mu1 = 0.
    print("WARNING: unrecognised mean parameter type provided - zero mean used.")
  }

  return(list(mu0=mu0,mu1=mu1))

}

#######################################
## calibrate model parameters (call)

calibrate_hetero = function(data,param,heteroModel,calc_rho=F,meantype,opt,strat=NULL){

  # need mean_type

  Qobs=data[[opt$obs]]
  Qh=data[[opt$pred]]
  #  Qh_T = vector(length=length(Qh))

  eta = calc_eta(Qobs=Qobs,Qh=Qh,param=param,heteroModel=heteroModel) # obs - simulated
  Qh_T = calc_tranz(Q=Qh,heteroModel=heteroModel,param=param) # The transformed simulated streamflow

  N = length(Qobs)

  date = as.Date(data$date,format='%d/%m/%Y')

  if (is.null(strat)){strat=set_strat_all(N)}
  
  # if (method=='MoM'){
  #   p = AR1_MoM(eta=eta,Qh=Qh_T,calc_rho=calc_rho,meantype)
  #   param$mean_eta_0 = p$mu0
  #   param$mean_eta_1 = p$mu1
  #   param$rho = p$rho
  #   param$sigma_y = p$sigma
  # } else {
  #   print("Invalid method selected - use MoM only")
  #   browser()
  # }

  # keep_list = list()
  # if (strat_mu_type == 'month'){
  #   for (m in 1:12){
  #     keep_list[[m]] = which(format(date,'%m')==m)
  #   }
  # } else {
  #   keep_list[[1]] = seq(1:length(date))
  # }

  # for (k in 1:length(keep_list)){
  #
  #   keep = keep_list[[k]]



  # }

  # strat_mu_type='month'
  #
  # mu0_vec = mu1_vec = vector(length = N)
  # mu0 = mu1 = c()
  # if (is.null(strat_mu_type)){
  #   k = 1
  #   keep = 1:N
  #   o = cal_mean(eta=eta[keep],Qh_T=Qh_T[keep],meantype=meantype); mu0[k] = o$mu0; mu1[k]=o$mu1
  #   mu0_vec[keep] = mu0; mu1_vec[keep] = mu1
  # } else if (strat_mu_type=='month'){
  #   for (k in 1:12){
  #     keep = which(as.integer(format(date,'%m'))==k)
  #     o = cal_mean(eta=eta[keep],Qh_T=Qh_T[keep],meantype=meantype); mu0[k] = o$mu0; mu1[k]=o$mu1
  #     mu0_vec = rep(mu0,N); mu0_vec = rep(mu1,N)
  #     mu0_vec[keep] = mu0; mu1_vec[keep] = mu1
  #   }
  # }

  mu0 = mu1 = c()
  mu0_vec = mu1_vec = vector(length = N)
  for (k in 1:length(strat$index[[strat$type$mean]])){
    keep = strat$index[[strat$type$mean]][[k]]
    omit = which(!(1:N)%in%keep)
    eta_tmp = eta; Qh_T_tmp = Qh_T
    eta_tmp[omit]= NA; Qh_T_tmp[omit] = NA
    o = cal_mean(eta=eta_tmp,Qh_T=Qh_T_tmp,meantype=meantype)
#    o = cal_mean(eta=eta[keep],Qh_T=Qh_T[keep],meantype=meantype)
    mu0[k] = o$mu0; mu1[k]=o$mu1
    mu0_vec[keep] = mu0[k]; mu1_vec[keep] = mu1[k]
  }

#    n = length(eta)-sum(is.na(eta))
  mu = mu0_vec+mu1_vec*Qh_T
  eta.star = eta-mu
  #s = sqrt((sum((eta.star)^2,na.rm=T))/n) # sigmaEta

  s = c()
  s_vec = vector(length = N)
  for (k in 1:length(strat$index[[strat$type$sigma]])){
    keep = strat$index[[strat$type$sigma]][[k]]
    omit = which(!(1:N)%in%keep)
    eta.star_tmp = eta.star
    eta.star_tmp[omit]= NA
    s[k] = sd(eta.star_tmp,na.rm=T)
    s_vec[keep] = s[k]
  }

#  s = sd(eta.star,na.rm=T)

  if (calc_rho){

    # ErrorlagForward <- eta.star[2:nlen]
    # ErrorlagBackward <- eta.star[1:nlen-1]
    # sb = sqrt((sum((ErrorlagBackward)^2,na.rm=T))/n)
    # sf = sqrt((sum((ErrorlagForward)^2,na.rm=T))/n)

    #    rho = (sum((ErrorlagForward)*(ErrorlagBackward),na.rm=T))/((n-1)*sb*sf) # autocorrelation
    #rho = cor(eta.star[2:N],eta.star[1:(N-1)],use = 'pairwise.complete.obs')

    rho = c()
    rho_vec = vector(length = N)
    for (k in 1:length(strat$index[[strat$type$rho]])){
      keep = strat$index[[strat$type$rho]][[k]]
      omit = which(!(1:N)%in%keep)
      eta.star_tmp = eta.star
      eta.star_tmp[omit]= NA
      rho[k] = cor(eta.star_tmp[2:N],eta.star_tmp[1:(N-1)],use = 'pairwise.complete.obs')
    }

  } else {
    rho = 0.
  }

  sigma = sqrt((s^2)*(1-(rho^2))) # sigmaY

  #    return(list(mu0=mu0,mu1=mu1,rho=rho,sigma=sigma,mu=mu))
  #  }

  param$mean_eta_0 = mu0
  param$mean_eta_1 = mu1
  param$rho = rho
  param$sigma_y = sigma

  return(param)
}

#######################################
## calibrate AR1 parameters

# AR1_MoM = function(eta,Qh=NULL,calc_rho=F,meantype){
#
#   if (meantype=="linear"){
#     m = lm(eta~Qh,na.action=na.omit)
#     mu0 = m$coefficients[1]
#     mu1 = m$coefficients[2]
#
#   } else if (meantype=="constant"){
#     mu0 = mean(eta,na.rm=T)
#     mu1 = 0.
#   } else if (meantype=="zero"){
#     mu0 = 0.
#     mu1 = 0.
#   } else{
#     mu0 = 0.
#     mu1 = 0.
#     print("WARNING: unrecognised mean parameter type provided - zero mean used.")
#   }
#
#   n = length(eta)-sum(is.na(eta))
#   nlen = length(eta)
#   mu = mu0+mu1*Qh
#   eta.star = eta-mu
#   #s = sqrt((sum((eta.star)^2,na.rm=T))/n) # sigmaEta
#   s = sd(eta.star,na.rm=T)
#
#   if (calc_rho){
#
#     # ErrorlagForward <- eta.star[2:nlen]
#     # ErrorlagBackward <- eta.star[1:nlen-1]
#     # sb = sqrt((sum((ErrorlagBackward)^2,na.rm=T))/n)
#     # sf = sqrt((sum((ErrorlagForward)^2,na.rm=T))/n)
#
# #    rho = (sum((ErrorlagForward)*(ErrorlagBackward),na.rm=T))/((n-1)*sb*sf) # autocorrelation
#     rho = cor(eta.star[2:nlen],eta.star[1:nlen-1],use = 'pairwise.complete.obs')
#
#   } else {
#     rho = 0.
#   }
#   sigma = sqrt((s^2)*(1-(rho^2))) # sigmaY
#
#   return(list(mu0=mu0,mu1=mu1,rho=rho,sigma=sigma,mu=mu))
# }
