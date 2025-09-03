#dyn.load('D:\\Data\\Work\\BOM_project\\UoA_Data70k_20140514\\c_crps.dll')


#' Cumulative Rank Probability Score (CRPS) 
#'
#' Calculate the CRPS as defined by Hersbach, H. Decomposition of the continuous ranked probability score for ensemble prediction systems Weather and Forecasting, 2000, 15, 559-570
#'
#' @param obs Vector of observed data (nval x 1)
#' @param pred Matrix of simulated ensemble data (nval x ne)
#' @return a vector of 14 values containing the performance criteria
#' @export
#' @useDynLib hydrodiy
#' @examples
#' n <- 5000
#' p <- 100
#' o <- arima.sim(n=n,list(ar=0.9),innov=rexp(n)) # Observation
#' e <- matrix(rnorm(n,0,sd(o)/3),n,p) # error
#' s <- matrix(o,n,p) + e # Ensemble 
#' cr <- crps(o,s)
#' # Reliability diagram
#' plot(cr$crpsmat[,5],cr$crpsmat[,1],main="Reliability",xlab="Obs. Freq",ylab="Sim. Freq")
crpsFunc <- function(obs,pred){#,perturb=F){

#browser()
  keep = !is.na(obs)
  obs = obs[keep]
  pred = pred[keep,]
  
#  if (perturb){
#    nObs = dim(pred)[1]
#    nReps = dim(pred)[2]
#    pred = pred + matrix(runif(length(pred),0,1e-8),ncol=nReps)
#    obs = obs + runif(nObs,0,1e-8)
#  }
  
	# Input check
	nval<-as.integer(length(obs))
	if(is.null(dim(pred))) pred <- matrix(pred,length(pred),nrow=1)
	if(nrow(pred)!=nval) stop("[0001] nrow(pred)!=nval")
	ncol<-as.integer(ncol(pred))
	
	if(ncol>0 & nval>0){
		# inputs
		obs<-as.double(obs)
		pred<-as.double(pred)
		w<-rep(1,nval)
	
		#outputs
		reliability_table<-as.double(rep(0,(ncol+1)*7))
		crps_decompos<-as.double(rep(0.0,5))
    out <- .C("c_crps", nval=nval,ncol=ncol,use_weights=1,is_sorted=0,obs=obs,sim=pred,weights_vector=w,
			reliability_table=reliability_table,crps_decompos=crps_decompos)
    out$reliability_table = matrix(out$reliability_table,ncol+1,7)
		colnames(out$reliability_table) = c("f","a","b","g","o","r","c")
		out$crps_decompos = matrix(out$crps_decompos,1,5)
		colnames(out$crps_decompos) = c("CRPS","Reliability","Resolution","Uncertainty","CRPSpot")
	}
	return(out);
}
