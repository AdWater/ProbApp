
set_strat_all = function(N){
  strat = list()
  strat$type = strat$index = list()
  ######
  strat$index$all = list()
  strat$index$all[[1]] = 1:N
  strat$type$mean=strat$type$sigma=strat$type$rho='all'
  return(strat)
}


