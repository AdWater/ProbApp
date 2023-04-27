#' probPred function
#'
#' Calibrates error model and generates probabilistic predictions from a given observed and simulated series.
#' @param data  Input data, given as a dataframe. Must contain unbroken, aligned columns of a) dates, b) observed streamflow, and c) simulated streamflow from a deterministic model. Missing values for observed data can be represent by \code{NA} or negative numbers.)
#' @param opt List of options - see details.
#' @param param List of error model parameters - see details.
#' @returns Output includes (depending on options specified in \code{opt})
##' \itemize{
##' \item{list containing calibrated error model parameters (\code{param}),
##'       predictive replicates (\code{pred.reps}), standardized residuals (\code{std.resids})
##'       and metrics (\code{metrics}) (when \code{opt$returnOutput=T})}
##' \item{.pdf summary file containing plots and metrics (when \code{opt$pdfOutput=T})}
##' \item{.csv file containing the probabilistic predictions (when \code{opt$repPrint=T})}
##' \item{.csv file containing probability limits of predictions (when \code{opt$plPrint=T})}
##' }
##'
##' @details The \code{opt} argument is a list with the following options
##' \itemize{
##'  \item{\code{reps}}{: number of probabilistic replicates to be generated. Recommended (and default) number is 100. Higher numbers of replicates require more computing time and lower numbers risk innacurate predictions.}
##'  \item{\code{title}}{: Default 'replicate'. The title (minus the extension) of the printed .pdf that is output from this function..}
##'  \item{\code{dirName}}{: directory into which the output files are to be saved.}
##'  \item{\code{obsName}}{: Column name of the observed streamflow in the input data. Default is 'obs'. This input must match the column header of the observed streamflow.}
##'  \item{\code{predName}}{: Column name of the simulated, deterministic streamflow in the input data. Default is 'pred'. This input must match the column header of the simulated streamflow.}
##'  \item{\code{dateName}}{: Column name of the dates in the input data. Default is 'date'. This input must match the column header of the dates.}
##'  \item{\code{meanType}}{: Options are "zero" or "linear", with "linear" being default. Defines the structure of the mean parameter in the probabilistic model. "zero" uses traditional modelling assumptions, and "linear" represents an innovative new approach that is more generally-applicable to a wider range of objective functions. Please refer to Hunter et al. 2021 for a comprehensive demonstration of the difference in predictive quality between the two.}
##'  \item{\code{unit}}{: Units of the input streamflow (used in plot labels). Default is 'mmd'. Alternatives are 'mm/d', 'm3s', 'm3/s','ML/d', 'MLd'.}
##'  \item{\code{repPrint}}{: TRUE / FALSE to print a .csv containing the probabilistic replicates. Default is FALSE.}
##'  \item{\code{plPrint}}{: TRUE / FALSE to print a .csv containing the probability limits. Default is FALSE.}
##'  \item{\code{plot}}{: TRUE / FALSE to produce figures. Default is FALSE.}
##'  \item{\code{pdfOutput}}{: TRUE / FALSE to print figures to pdf. Default is FALSE.}
##'  \item{\code{returnOutput}}{: TRUE / FALSE to return parameters, predictive replicates, standardized residuals and metric values. Default is F.}
##' }
##' @details The \code{param} argument is a list of transformation parameters used in the error model, including
##' \itemize{
##'  \item{\code{A}}{: Box-Cox shift parameter. Default is 0.}
##'  \item{\code{lambda}}{: Box-Cox power parameter. Default is 0.2.}
##' }
##' Note \code{param} can also include other calibrated error model parameters obtained from another call to probPed() - see example 2 below.
#' @export
#' @examples
#' # Read example data file
#' fname = system.file("extdata", "402204_SLS.csv", package = "ProbPred")
#' data = read.csv(fname) # Example of setting up the data file.
#'
#' # -------------------------------------------------------------------------
#' # Example 1 - calibration and prediction using all data
#' # Setup options list. (All other options use default values)
#' opt = list(title='myProbPredictions.all',   # title of output
#'            dirName='./',                    # output in cwd
#'            pdfOutput=TRUE,                     # generate PDF output
#'            returnOutput=FALSE)                    # don't output data list
#'
#' # Setup (transformation) parameter values
#' param = list(A=0,lambda=0.2)
#'
#' # Run probPred
#' # This produces PDF output file in current working directory
#' probPred(data=data,opt=opt,param=param)
#'
#' # -------------------------------------------------------------------------
#' # Example 2 - separate calibration and validation (prediction) periods
#'
#' # Use 1st half of data for calibration, and 2nd half for validation
#' nT = nrow(data)
#' data.cal = data[1:(nT/2),]               # calibration data
#' data.val = data[(nT/2+1):nT,]            # validation data
#'
#' # ---- Calibration ----
#' # Setup calibration run options. Return output data from probPred
#' # including parameters (but don't produce PDF or save replicates to file)
#' opt.cal = list(returnOutput = TRUE)
#'
#' # Setup (transformation) parameter values
#' param.cal = list(A=0,lambda=0.2)

#' # Run probPred for calibration period (producing parameter values)
#' out.cal=probPred(data=data.cal,opt=opt.cal,param=param.cal)
#'
#' # ---- Validation ----
#' # Setup validation run options.
#' opt.val = list(title = 'myProbPredictions.val', # output name for val period
#'                dirName='./',                    # output in cwd
#'                pdfOutput=TRUE,                      # create PDF output this time
#'                returnOutput = FALSE)               # don't return data for val
#'
#' # Use calibrated parameter values from calibration period
#' param.val = out.cal$param
#'
#' # Run probPred for validation period (producing PDF output)
#' out.val=probPred(data=data.val,opt=opt.val,param=param.val)

probPred = function(data,opt=NULL,param=NULL) {

#######################################
## Inputs & parameters

  if(is.null(opt$reps)){opt$reps=100}
  if(is.null(opt$title)){opt$title='replicate'}
  if(is.null(opt$dirName)){opt$dirName='./'}
  if(is.null(opt$obsName)){opt$obsName='obs'}
  if(is.null(opt$predName)){opt$predName='pred'}
  if(is.null(opt$dateName)){opt$dateName='date'}
  if(is.null(opt$meanType)){opt$meanType='linear'}
  if(is.null(opt$unit)){opt$unit='mmd'}
  if(is.null(opt$repPrint)){opt$repPrint=F}
  if(is.null(opt$plPrint)){opt$plPrint=F}
  if(is.null(opt$plot)){opt$plot=F}
  if(is.null(opt$pdfOutput)){opt$pdfOutput=F}
  if(opt$pdfOutput){opt$plot=T}
  if(is.null(opt$returnOutput)){opt$returnOutput=T}

  if(is.null(param$A)){param$A=0}
  if(is.null(param$lambda)){param$lambda=0.2}

  setwd(opt$dirName)
  data_dirname = system.file("shiny",package="ProbPred")

#######################################
## Determine if we calculate parameters

  heteroModel = 'BC'

  if (!(is.null(param$mean_eta_0)|is.null(param$mean_eta_1)|is.null(param$rho)|is.null(param$sigma_y))){
    calc.params = F
  } else {
    calc.params = T
    paramFix = list(A=param$A,lambda=param$lambda)
    meantype = opt$meantype
  }

#######################################
## Error checks on the input data

  data.headers = colnames(data)

  if(!is.element(opt$obsName,data.headers) |
     !is.element(opt$predName,data.headers) |
     !is.element(opt$dateName,data.headers)) {
    xerr(flag=1,shiny=F)# headers not in the data file
    return()
  } else if(all(is.na(data[opt$obsName][[1]])) | all(is.na(data[opt$predName][[1]])) | all(is.na(data[opt$dateName][[1]]))) {
    xerr(flag=2,shiny=F) # checks for empty data vectors
    return()

  } else if(is.character(data[opt$obsName][[1]]) | is.character(data[opt$predName][[1]]) | is.numeric(data[opt$dateName][[1]])) {
    xerr(flag=3,shiny=F) # check for characters (e.g. dates) in the obs or pred
    return()
  } else if(sum(data[opt$obsName][[1]]-data[opt$predName][[1]],na.rm=T)==0) {
    xerr(flag=4,shiny=F) # check for obs and pred being the same vector
    return()
  } else if(length(data[opt$obsName][[1]])!=length(data[opt$predName][[1]]) |
            length(data[opt$obsName][[1]])!=length(data[opt$dateName][[1]])) {
    xerr(flag=5,shiny=F) # checks that the input vectors are all the same length
    return()
  }

  # set missing data to 'NA' - important for plotting
  data[[opt$obsName]][data[[opt$obsName]]<0 | is.infinite(data[[opt$obsName]]) | is.nan(data[[opt$obsName]])] = NA
  data[[opt$predName]][data[[opt$predName]]<0 | is.infinite(data[[opt$predName]]) | is.nan(data[[opt$predName]])] = NA

######################################
## Calculations

  if (calc.params){
    # calc parameters
    print('Calibrating parameters')
    param = calibrate_hetero(data=data,param=paramFix,heteroModel=heteroModel,calc_rho=T,meantype=opt$meanType,opt=opt)
  } else {
    print('Using provided parameter values (i.e. not calibrating error model)')
  }

  # calc eta_star
  std.resids = calc_std_resids(data=data,param=param,heteroModel=heteroModel,opt=opt)

  print("Starting calculation of probabilistic replicates...")
  # calc predictive replicates
  pred.reps = calc_pred_reps(Qh=data[[opt$pred]],heteroModel=heteroModel,param=param,nReps=opt$reps,Qmin=0.,Qmax=999.,truncType='spike')

  # calc probability limits
  pred.pl = calc.problim(pred.reps,percentiles=c(0.05,0.25,0.5,0.75,0.95))

  print("Starting calculation of metrics...")
  # generating metrics (reliability, precision, bias)
  metrics = calc_metrics(data=data,pred.reps=pred.reps,opt=opt)

  # opening pdf
  if(opt$pdfOutput){
    print("Printing to pdf...")
    pdf(paste(opt$title,".summary.pdf",sep=""))
  }

  # decide whether to plot
  if(opt$plot){

######################################
## producing plots

    # Front page
    output.main(param=param,metrics=metrics,data=data,is.data=T,opt=opt,dir.loc=data_dirname)

    # Boxplots
    boxplotter(data_dirname=data_dirname,catchmentMetric=metrics$reliability,metric="reliability",boxColour="pink")
    boxplotter(data_dirname=data_dirname,catchmentMetric=metrics$sharpness,metric="sharpness",boxColour="white")
    boxplotter(data_dirname=data_dirname,catchmentMetric=metrics$bias,metric="bias",boxColour="lightblue")

    # PQQ plot
    plot.performance(data=data,pred.reps=pred.reps,type='PQQ',opt=opt)

    # Residual plot #1
    plot.residuals(data=data,std.resids=std.resids,type='pred',opt=opt)

    # Residual plot #2
    plot.residuals(data=data,std.resids=std.resids,type='prob(pred)',opt=opt)

    # Density plot
    plot.residuals(data=data,std.resids=std.resids,type='density',opt=opt)

    # standardised residual plot
    tranzplotter(data=data,param=param,heteroModel=heteroModel,add.legend=T,add.title=T,opt=opt)

    # auto & partial correlation plots - temporarily unable to handle missing data
    #if (!is.na(min(data[[opt$obs]])) && !is.na(min(data[[opt$pred]]))) {
      #acfplotter(data=data,acfType='acf',param=param,heteroModel=heteroModel,opt=opt)
      #acfplotter(data=data,acfType='pacf',param=param,heteroModel=heteroModel,opt=opt)
    #}

    # Timeseries
    timeseries(data=data,pred.reps=pred.reps,opt=opt)

  }

  # terminate PDF
  if(opt$pdfOutput){dev.off()}

######################################
## Printing .csv

  if(opt$repPrint) { # Print out a .csv with the replicates in it
    write.csv(x=pred.reps,file=paste(title,"_replicates.csv",sep=""))
  }
  if(opt$plPrint) { # Print out a .csv with the probability limits in it
    write.csv(x=pred.pl,file=paste(title,"_probLimits.csv",sep=""))
  }
  print("Run complete!")

  if(opt$returnOutput){
    return(list(param=param,pred.reps=pred.reps,std.resids=std.resids,metrics=metrics))
  }

}


