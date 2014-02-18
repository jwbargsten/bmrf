as.weighted <- function(bmrf, method=c("degree.max","degree.prod", "degree.add", "degree.meandev") ) {

  if(is.function(method)) {
    weight.fun <- method
  } else {
    method <- match.arg(method)

    weight.fun <- switch(method,
      degree.max = .weight.degree.max,
      degree.prod = .weight.degree.prod,
      degree.add = .weight.degree.add,
      degree.meandev = .weight.degree.meandev
    )
  }

  deg <- rowSums(bmrf@net)
  bmrf@net <- bmrf@net * matrix(mapply(
      weight.fun,
      row(bmrf@net),
      col(bmrf@net),
      MoreArgs=list(bmrf@net, deg)),
    nrow=nrow(bmrf@net)
    )
  singelton.idcs <- which(rowSums(bmrf@net) == 0)
  bmrf@net <- bmrf@net[-singelton.idcs, -singelton.idcs]
  bmrf@fd <- bmrf@fd[-singelton.idcs,]
  bmrf@go <- bmrf@go[-singelton.idcs,]
  bmrf@unknown.idcs <- which(rowSums(bmrf@go) == 0)
  bmrf
}

.weight.degree.max <- function(i, j, net, deg) { 1/max(deg[i], deg[j]) }
.weight.degree.prod <- function(i, j, net, deg) { 1/(deg[i] * deg[j]) }
.weight.degree.add <- function(i, j, net, deg) { 1/(deg[i] + deg[j]) }
.weight.degree.meandev <- function(i, j, net, deg) { 
  m <- mean(deg)
  if(max(deg[i], deg[j]) > 3*m) {
    return(0)
  }
  return(1)
}
