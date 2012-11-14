## predict using a bmrf object
predict.bmrf <- function(b) {

  for(i in 1:ncol(b@go)) {
    go.proteins <- b@go[,i]
    go.proteins[b@unknown.idcs] = -1

    glm.enet.pred <- predict.glmnet(bmrf=b, go.idx=i, 
    dfmax = (ncol(b@fd)-1)
    )

    if(is.null(glm.enet.pred))
      next
  }
}

predict.bmrf.glmnet <- function(bmrf, go.idx, dfmax) {
  yf <- bmrf@go[-bmrf@unknown.idcs,go.idx]
  xf <- bmrf@fd[-bmrf@unknown.idcs,]

  #Fit for the common set of proteins (since for them there is Y and X for the regression fitting step)
  f = try(
    cv.glmnet(
      y = yf,
      x = xf,
      family = "binomial",
      dfmax = dfmax,
      alpha=0.5,
      maxit=1000,
      type.measure="auc"
    )
  )

  if(class(f) == 'try-error')
    return(NULL)

  yhat = try(predict(f, bmrf@fd, s = "lambda.min"))
  names(Yhat) = rownames(bmrf@fd);

  if(class(yhat) == 'try-error')
    return(NULL)

  return(yhat);
}
