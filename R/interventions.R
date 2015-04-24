getCauses <- function(g, v){  
  if(is.null(V(g)$intervention)) stop("Interventions have not been set")
  if(V(g)[v]$intervention) return(V(g)[0])
  V(g)[nei(v, mode="in")]
}

setValuesFunc <- function(val){
  function(g, v){
    if(V(g)[v]$intervention){
      V(g)[v]$output.signal <- list(rep(val, g$n))
    }else{
      g <- calculateVals(g, v)
    }
    g
  }
}

fixNodes <- function(g, v.set, val){
  V(g)$intervention <- FALSE
  V(g)[v.set]$intervention <- TRUE
  input.intervention <- as.logical((V(g)$type == "input") * V(g)$intervention) 
  V(g)[input.intervention]$output.signal <- list(rep(val, g$n))
  setValCallback <- setValuesFunc(val)
  g <- updateVertices(g, getDeterminers = getCauses, callback = setValCallback)
  g
}