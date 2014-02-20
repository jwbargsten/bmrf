transformNet <- function(bmrf, method=c("degree.max","degree.prod", "degree.add", "degree.meandev", "degree.component"), factor=1 ) {
  if(is.function(method)) {
    weight.fun <- method
  } else {
    method <- match.arg(method)

    weight.fun <- switch(method,
      degree.max = .weight.degree.max,
      degree.prod = .weight.degree.prod,
      degree.add = .weight.degree.add,
      degree.meandev = .weight.degree.meandev,
      degree.component = .weight.degree.component
    )
  }

  deg <- rowSums(bmrf@net)
  bmrf@net <- bmrf@net * matrix(mapply(
      weight.fun,
      row(bmrf@net),
      col(bmrf@net),
      MoreArgs=list(net=bmrf@net, deg=deg, factor=factor, deg.mean=mean(deg), deg.median=median(deg))),
    nrow=nrow(bmrf@net)
    )

  bmrf
}

clean <- function(bmrf) {
  edgeless.idcs <- which(rowSums(bmrf@net) == 0)

  if(length(edgeless.idcs ) > 0) {
    bmrf@net <- bmrf@net[-edgeless.idcs, -edgeless.idcs]
    bmrf@fd <- bmrf@fd[-edgeless.idcs,]
    bmrf@go <- bmrf@go[-edgeless.idcs,]
    bmrf@go <- bmrf@go[,colSums(bmrf@go) > 0]
    bmrf@unknown.idcs <- which(rowSums(bmrf@go) == 0)
  }
  bmrf
}


.weight.degree.max <- function(i, j, net, deg) { 1/max(deg[i], deg[j]) }
.weight.degree.prod <- function(i, j, net, deg) { 1/(deg[i] * deg[j]) }
.weight.degree.add <- function(i, j, net, deg) { 1/(deg[i] + deg[j]) }
.weight.degree.meandev <- function(i, j, net, deg, factor=3, deg.mean, deg.median) { 
  if(max(deg[i], deg[j]) > factor*deg.mean) {
    return(0)
  }
  return(1)
}
.weight.degree.mediandev <- function(i, j, net, deg, factor=3, deg.mean, deg.median) { 
  if(max(deg[i], deg[j]) > deg.median*factor) {
    return(0)
  }
  return(1)
}

.weight.degree.component <- function(i, j, net, deg, factor=3, deg.mean, deg.median) { 
  ## two high degree nodes
  if(deg[i] > factor*deg.median && deg[j] > factor*deg.median)  {
    return(1)
  ## two low degree nodes
  } else if(deg[i] <= factor*deg.median && deg[j] <= factor*deg.median)  {
    return(1)
  ## two different degree nodes
  } else {
    return(0)
  }
}
