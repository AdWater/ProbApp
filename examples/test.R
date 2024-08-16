rm(list=ls())
devtools::load_all()

data = read.csv('data/402204_SLS.csv',as.is=T) ## Example of setting up the data file.

opt = list(obs='obs',
           pred='pred',
           date='date',
           meantype='linear',
           dirname='.',
           pdfOutput=F)                        ## Example of setting up the user options.

param = list(A=0.,
             lambda=0.2)                     ## Example of setting up the input parameters.

opt$title = 'myProbPredictions.all'
opt$pdfOutput = T
#out=probPred(data=data,opt=opt,param=param)
out=probPred(data=data)

nT = nrow(data)
data.cal = data[1:(nT/2),]                 # calibration data
data.val = data[(nT/2+1):nT,]            # validation data

opt.cal = opt
opt.cal$title = 'myProbPredictions.cal'
opt.cal$return.output = T # return output data from probPred (including parameters)
out.cal=probPred(data=data.cal,opt=opt.cal,param=param)

opt.val = opt
opt.val$title = 'myProbPredictions.val'
out.val=probPred(data=data.val,opt=opt.val,param=out.cal$param) # run probPred, using parameters from calibration period
