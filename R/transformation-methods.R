as.weighted <- function(bmrf, method=c("degree.max","degree.prod") ) {
  method <- match.arg(method)

  weight.fun <- switch(method,
    degree.max = .weight.degree.max,
    degree.prod = .weight.degree.prod
  )

  cs <- colSums(bmrf@net)
  bmrf@net <- bmrf@net * matrix(mapply(
      weight.fun,
      row(bmrf@net),
      col(bmrf@net),
      MoreArgs=list(cs)),
    nrow=nrow(bmrf@net)
    )
  bmrf
}

.weight.degree.max <- function(i, j, cs) { 1/max(cs[i], cs[j]) }
.weight.degree.prod <- function(i, j, cs) { 1/(cs[i] * cs[j]) }
