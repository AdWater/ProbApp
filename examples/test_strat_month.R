rm(list=ls())

devtools::load_all('C:/Users/a1065639/Work/ProbApp/')

################################################

setup_par_mat_list = function(o,opt){
  par_mat_list = list()
  par_mat_list$mean_eta_0 = matrix(nrow=length(o$param$mean_eta_0),ncol=opt$reps)
  par_mat_list$mean_eta_1 = matrix(nrow=length(o$param$mean_eta_1),ncol=opt$reps)
  par_mat_list$rho = matrix(nrow=length(o$param$rho),ncol=opt$reps)
  par_mat_list$sigma_y = matrix(nrow=length(o$param$sigma_y),ncol=opt$reps)
  return(par_mat_list)
}

add_par_mat_list = function(par_mat_list,o,r){
  par_mat_list$mean_eta_0[,r] = o$param$mean_eta_0 
  par_mat_list$mean_eta_1[,r] = o$param$mean_eta_1 
  par_mat_list$rho[,r] = o$param$rho 
  par_mat_list$sigma_y[,r] = o$param$sigma_y 
  return(par_mat_list)
}

plot_par_mat_list = function(par_mat_list,o){
  par(mfrow=c(2,2),mar=c(1,3,3,1))
  boxplot.ext(t(par_mat_list$mean_eta_0),colouring = 'white',main='mean_eta_0',xaxt='n')
  points(x=1:nrow(par_mat_list$mean_eta_0),y=o$param$mean_eta_0,col='red',pch=16)
  boxplot.ext(t(par_mat_list$mean_eta_1),colouring = 'white',main='mean_eta_1',xaxt='n')
  points(x=1:nrow(par_mat_list$mean_eta_1),y=o$param$mean_eta_1,col='red',pch=16)
  boxplot.ext(t(par_mat_list$rho),colouring = 'white',main='rho',xaxt='n')
  points(x=1:nrow(par_mat_list$rho),y=o$param$rho,col='red',pch=16)
  boxplot.ext(t(par_mat_list$sigma_y),colouring = 'white',main='sigma_y',xaxt='n')
  points(x=1:nrow(par_mat_list$sigma_y),y=o$param$sigma_y,col='red',pch=16)
}

synth_test = function(o,opt,param){
  
  par_mat_list = setup_par_mat_list(o,opt)
  for (r in 1:opt$reps){
    print(r)
    data.synth.o = data; data.synth.o$obs = o$pred.reps[,r]
    o.synth = probPred(data=data.synth.o,opt=opt,param=param)
    par_mat_list = add_par_mat_list(par_mat_list,o.synth,r)
  }
  plot_par_mat_list(par_mat_list,o)
  
}

################################################

# Read example data file
fname = system.file("extdata", "402204_SLS.csv", package = "ProbPred")
data = read.csv(fname) # Example of setting up the data file.

#data1 = rbind(data,data,data,data,data,data,data,data,data,data)
#data1$date = as.Date(data1$date[1],format ='%d/%m/%Y') + 0:(length(data1$date)-1)
#data = data1

# -------------------------------------------------------------------------
# Example 1 - calibration and prediction using all data
# Setup options list. (All other options use default values)
opt = list(title='myProbPredictions.all',   # title of output
           dirName='./',                    # output in cwd
           pdfOutput=FALSE,                     # generate PDF output
           returnOutput=TRUE,               #don't output data list
           strat_mean=NULL,reps=20)

# Setup (transformation) parameter values
param = list(A=0,lambda=0.2)

# Run probPred
# This produces PDF output file in current working directory
o = probPred(data=data,opt=opt,param=param)

std_resid = o$std.resids
date = as.Date(data$date,format='%d/%m/%Y')
std_resid_zoo = zoo(std_resid,order.by=date)
std_resid_month_zoo= aggregate(std_resid_zoo, format(time(std_resid_zoo), "%m"), mean)
plot(std_resid_month_zoo,type='o')
abline(h=0,lty=2)

synth_test(o,opt,param)

# mean_eta_0_mat = matrix(nrow=length(o$param$mean_eta_0),ncol=nReps)
# mean_eta_1_mat = matrix(nrow=length(o$param$mean_eta_1),ncol=nReps)
# rho_mat = matrix(nrow=length(o$param$rho),ncol=nReps)
# sigma_y_mat = matrix(nrow=length(o$param$sigma_y),ncol=nReps)

# par_mat_list = setup_par_mat_list(o)
# for (r in 1:nReps){
#   print(r)
#   data.synth.o = data; data.synth.o$obs = o$pred.reps[,r]
#   o.synth = probPred(data=data.synth.o,opt=opt,param=param)
#   par_mat_list = add_par_mat_list(par_mat_list,o.synth,r)
#   # mean_eta_0_mat[,r] = o.synth$param$mean_eta_0 
#   # mean_eta_1_mat[,r] = o.synth$param$mean_eta_1 
#   # rho_mat[,r] = o.synth$param$rho 
#   # sigma_y_mat[,r] = o.synth$param$sigma_y 
# }
# plot_par_mat_list(par_mat_list,o)
 
  
# par(mfrow=c(2,2),mar=c(1,3,3,1))
# 
# boxplot.ext(t(mean_eta_0_mat),colouring = 'white',main='mean_eta_0')
# # abline(h=o$param$mean_eta_0)
# points(x=1,y=o$param$mean_eta_0,col='red',pch=16)
# boxplot.ext(t(mean_eta_1_mat),colouring = 'white',main='mean_eta_1')
# #abline(h=o$param$mean_eta_1)
# points(x=1,y=o$param$mean_eta_1,col='red',pch=16)
# boxplot.ext(t(rho_mat),colouring = 'white',main='rho')
# #abline(h=o$param$rho)
# points(x=1,y=o$param$rho,col='red',pch=16)
# boxplot.ext(t(sigma_y_mat),colouring = 'white',main='sigma_y')
# #abline(h=o$param$sigma_y)
# points(x=1,y=o$param$sigma_y,col='red',pch=16)

opt1 = opt; opt1$strat_mean = 'month'
o1 = probPred(data=data,opt=opt1,param=param)
synth_test(o1,opt1,param)

# data.synth.o1 = data; data.synth.o1$obs = o1$pred.reps[,1]
# o1.synth = probPred(data=data.synth.o1,opt=opt1,param=param)


opt2 = opt1; opt2$strat_sigma = 'month'
o2 = probPred(data=data,opt=opt2,param=param)
synth_test(o2,opt2,param)

# data.synth.o2 = data; data.synth.o2$obs = o2$pred.reps[,1]
# o2.synth = probPred(data=data.synth.o2,opt=opt2,param=param)

opt3 = opt2; opt3$strat_rho = 'month'
o3 = probPred(data=data,opt=opt3,param=param)

std_resid = o3$std.resids
date = as.Date(data$date,format='%d/%m/%Y')
std_resid_zoo = zoo(std_resid,order.by=date)
std_resid_month_zoo= aggregate(std_resid_zoo, format(time(std_resid_zoo), "%m"), mean)
plot(std_resid_month_zoo,type='o')
abline(h=0,lty=2)


synth_test(o3,opt3,param)


# data.synth.o3 = data
# data.synth.o3$obs = o3$pred.reps[,1]
# 
# o3.synth = probPred(data=data.synth.o3,opt=opt3,param=param)

#####################
# check rho

dates = as.Date(data$date,format='%d/%m/%Y')
mons = as.integer(format(dates,'%m'))

make.na = which(!(mons%in%c(6,7,8)))
data.winter = data
data.winter$obs[make.na] = NA
o.winter = probPred(data=data.winter,opt=opt,param=param)

make.na = which(!(mons%in%c(12,1,2)))
data.summer = data
data.summer$obs[make.na] = NA
o.summer = probPred(data=data.summer,opt=opt,param=param)
