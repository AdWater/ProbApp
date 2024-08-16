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

out=probPred(data=data,opt=opt,param=param)
#out=probPred(data=data)

data.synth = data
rep = 1
data.synth$obs = out$pred.reps[,rep]
out.synth=probPred(data=data.synth,opt=opt,param=param)
