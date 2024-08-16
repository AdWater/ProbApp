rm(list=ls())

devtools::load_all()


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
           strat_mean=NULL,reps=100)

# Setup (transformation) parameter values
param = list(A=0,lambda=0.2)

# Run probPred
# This produces PDF output file in current working directory
o = probPred(data=data,opt=opt,param=param)

data.synth.o = data; data.synth.o$obs = o$pred.reps[,1]
o.synth = probPred(data=data.synth.o,opt=opt,param=param)

opt1 = opt; opt1$strat_mean = 'month'
o1 = probPred(data=data,opt=opt1,param=param)

data.synth.o1 = data; data.synth.o1$obs = o1$pred.reps[,1]
o1.synth = probPred(data=data.synth.o1,opt=opt1,param=param)


opt2 = opt1; opt2$strat_sigma = 'month'
o2 = probPred(data=data,opt=opt2,param=param)

data.synth.o2 = data; data.synth.o2$obs = o2$pred.reps[,1]
o2.synth = probPred(data=data.synth.o2,opt=opt2,param=param)

opt3 = opt2; opt3$strat_rho = 'month'
o3 = probPred(data=data,opt=opt3,param=param)

data.synth.o3 = data
data.synth.o3$obs = o3$pred.reps[,1]

o3.synth = probPred(data=data.synth.o3,opt=opt3,param=param)
