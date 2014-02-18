## predict using a bmrf object
predict.bmrf <- function(b, burnin=20, niter=20, file=NULL, format=c("3col", "long", "h5"), verbose=FALSE ) {

  format <- match.arg(format)
  o <- options()

  result <- NULL
  con <- NULL

  ## write header to output file if desired
  if(!is.null(file) && format != "h5") {
    con <- file(file, "w")
    if(format == "long")
      cat("go", paste0("\t", rownames(b@go)), "\n", sep="", file=con)
    else
      cat("", sep="", file=con)
    flush(con)
  } else {
    result <- Matrix(0, ncol=ncol(b@go), nrow=nrow(b@go), dimnames=dimnames(b@go))
  }

  if(verbose)
    message(paste("Number of GO-terms to predict: ", ncol(b@go), sep=""))

  for(i in 1:ncol(b@go)) {
    go.proteins <- b@go[,i]
    go.proteins[b@unknown.idcs] = -1
    term_name <- colnames(b@go)[i]

		## make the elastic net fitting and predictions for the functional domains (fd).
    ## set it to NULL if error, NA is used for situations where the fd stuff is excluded on purpose
    if(all(is.na(b@fd))) {
      warning(term_name, ": functional domains are not used for prediction", call.=FALSE)
      glm.enet.pred <- rep(0, length(go.proteins))
      names(glm.enet.pred) <- names(go.proteins)
    } else {
      glm.enet.pred <- predict.bmrf_glmnet(bmrf=b, go.idx=i, dfmax = (ncol(b@fd)-1))
      if(is.null(glm.enet.pred)) {
        warning(term_name, ": could not fit elastic net, skipping", call.=FALSE)
        next
      }
    }


    ## run BMRFz
    p = try(predict.bmrf_go(bmrf=b, go.proteins=go.proteins, fd.predicted=glm.enet.pred, burnin=burnin, niter=niter, term_name=term_name))
    if(class(p) == 'try-error')
      next

    p = .calibrate_logitR(p)

    if(is.null(file) || is.null(con)) {
      result[,i] <- p
    } else {
      switch(format,
        "3col" = cat(
            paste(names(p), rep(term_name, length(p)),  p, sep="\t"),
            sep="\n",
            file=con
        ),
        "long" = cat(
            colnames(b@go)[i],
            paste0("\t", p),
            "\n",
            sep="",
            file=con
        )
      )
      flush(con)
    }
    if(verbose) cat(".", sep="", file=stderr())
  }
  options(o)

  if(verbose) message(paste(" finished\n", sep=""))

  if(!is.null(con))
    close(con)

  if(format=="h5" && !is.null(file))
    .write_h5(file, result)

  result
}

predict.bmrf_go <- function(bmrf, go.proteins, fd.predicted, burnin, niter, term_name, initZlength=30) {
#
#	#---------------------------------------------------------
#	# Construct the schedules for reporting things
#	# genrep.schd: report the iteration
#	# scaleupd.schd: update the scale parameter for the proposal distribution
#	# simstor.schd: store the values
#	#--------------------------------------------------------
#
  A <- bmrf@net
  L <- go.proteins
	Dalpha = fd.predicted
	titer = burnin + niter;
	genrep.schd    = seq(from = 1, to = titer, by = 100);
	simstor.schd   = seq(from = 1, to = titer, by = 5);
	Z.schd 		   = seq(from = 1, to = titer, by = 1);

	# Degree of each protein
	dA = rowSums(A);

	names(dA) = rownames(A);
	unknowns = which(L == -1);
	knowns	= which(L >= 0);


	#Domains
  #D <- fd.predicted
#	Dalpha = glmnetDalpha(L,D, 10)
	# Initialize MRF parameters

  ## the number of neightbours with label
	M1 = as.vector(A[knowns,knowns] %*% L[knowns]);
  ## NS = neighbours
	NS = as.vector(dA[knowns]);
	Lk = as.vector(L[knowns])
	Dalphak = as.vector(Dalpha[knowns]);

	regtable = as.data.frame(cbind(Lk, M1,NS, Dalphak));
	colnames(regtable) = c("L", "M1", "NS", "Dalpha");


	regtable.fit = brglm(L ~ M1 + NS + Dalpha,
						family=binomial(link = "logit"),
						method = "brglm.fit", data=regtable);

  if(any(is.na(regtable.fit$coefficients)))
    warning(term_name, ": brglm has NA coefficients - ",
      paste0(names(regtable.fit$coefficients)[is.na(regtable.fit$coefficients)], collapse=", "),
      call.=FALSE
    )

  if (!regtable.fit$converged) {
    warning(term_name, ": brglm fitting did not converge", call.=FALSE)
  }

  rt.coeff <- regtable.fit$coefficients
  rt.vcov <- vcov(regtable.fit)

  Z <- NULL

  ## if Dalpha is NA, it contains no information and gets excluded from the model (aka set to 0)
  if(is.na(rt.coeff["Dalpha"])) {
    warning(term_name, ": rmvnorm dimensions of mean and sigma are not consistent, Dalpha is NA. Excluding fd.", call.=FALSE)
    Z = rmvnorm(mean=rt.coeff[- which(names(rt.coeff) == "Dalpha")], sigma=rt.vcov, n=initZlength);
    Z <- cbind(Z, Dalpha=0)

    rt.coeff["Dalpha"] <- 0
  } else {
    Z = rmvnorm(mean=rt.coeff, sigma=rt.vcov, n=initZlength);
  }

	MRFparams = as.vector(rt.coeff);


	#Initialization of Labelling
	Auk = A[unknowns,knowns];

	M1 = as.vector(Auk %*% L[knowns]);
	NSu  =  dA[unknowns];

	y = MRFparams[1] + M1*MRFparams[2] + NSu*MRFparams[3] + Dalpha[unknowns]*MRFparams[4];

	d = (1/(1 + exp(-y)))  - runif(length(unknowns));
	L[unknowns[d >= 0]] = 1;
	L[unknowns[d  < 0]] = 0;

	e.sigma = matrix(ncol= ncol(Z), nrow= ncol(Z));
	e.sigma[] = 0;
	diag(e.sigma) = 0.0001;

	Au = A[unknowns,];

	probs = vector(mode = "numeric", length = length(L));
	names(probs) = names(L);
	probs[]  = 0;
	counter  = 0;

	for(t in 1:titer)
	{

		#t1 = proc.time();

		M1 = as.vector(Au %*% L);

#XXX

		y =  MRFparams[1] + M1*MRFparams[2] + NSu*MRFparams[3] + Dalpha[unknowns]*MRFparams[4];

		d = (1/(1 + exp(-y)))  - runif(length(unknowns));
		L[unknowns[d >= 0]] = 1;
		L[unknowns[d  < 0]] = 0;

		# Update MRFparams. Propose a candidate

		s = sample(1:nrow(Z),2, replace=F);
		e = rmvnorm(mean = c(rep(0,ncol(Z))), sigma = e.sigma,n=1);
		gamma = runif(min = 0.2, max = 1,n=1);
		MRFparamsP = as.vector(MRFparams + (gamma*(Z[s[1],] - Z[s[2],])) + e);

		M1 = as.vector (A %*% L);
		yc = as.vector(MRFparams[1]  + M1*MRFparams[2] +  dA*MRFparams[3] + Dalpha*MRFparams[4]);
		yp = as.vector(MRFparamsP[1] + M1*MRFparamsP[2] + dA*MRFparamsP[3] + Dalpha*MRFparamsP[4]);

		expc = exp(-yc);
		cpc = 1/(1 + expc);
		cpc.1 = cpc;
		cpc[L == 0] = 1 - cpc[L == 0];
		cpc[cpc < 1e-12] = 1e-12;
		logc = log(cpc);

		expp = exp(-yp);
		cpp = 1/(1 + expp);
		cpp.1 = cpp;
		cpp[L == 0] = 1 - cpp[L == 0];
		cpp[cpp < 1e-12] = 1e-12;
		logp = log(cpp);
		log_r = sum(logp)  - sum(logc);
		if(log_r >= log(runif(1))){MRFparams = MRFparamsP;}

		# Check if it is time to report:
		if(is.element(t,Z.schd)){Z = rbind(Z, MRFparams);}
		if(is.element(t,simstor.schd))
		{
			if(t > burnin)
			{
				probs = probs + cpc.1;
				counter = counter + 1;
			}
		}
#		print(proc.time() - t1);
	}

	probs = probs/counter;
	return(probs);
}

predict.bmrf_glmnet <- function(bmrf, go.idx, dfmax) {
  yf <- bmrf@go[-bmrf@unknown.idcs,go.idx]
  xf <- bmrf@fd[-bmrf@unknown.idcs,]

  if(diff(range(yf)) == 0) {
    warning(
            "cv.glmnet: ",
            ifelse(yf[1] == 1, "ALL", "NO"),
            " proteins with existing labels (known proteins) are associated with ",
            colnames(bmrf@go)[go.idx],
            ", skipping term.",
            call.=FALSE
    )
    return(NULL)
  }
  #cat(sum(yf), sum(xf), "\n", sep="  ")

  #Fit for the common set of proteins (since for them there is Y and X for the regression fitting step)
  f = try(
    suppressWarnings(
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
  )

  if(class(f) == 'try-error') {
    warning("cv.glmnet fitting not successful", call.=FALSE)
    return(NULL)
  }


  yhat = try(suppressWarnings(predict(f, bmrf@fd, s = "lambda.min")))
  names(yhat) = rownames(bmrf@fd);

  if(class(yhat) == 'try-error') {
    warning("cv.glmnet prediction not successful", call.=FALSE)
    return(NULL)
  }

  return(yhat);
}

.calibrate_logitR = function(P)
{
	#Calibrate the input probabilities
    P[P >= 1 - (1e-10)] = 1- (1e-10);
    P[P <= (1e-10)] =  (1e-10);
	Ppriors = mean(P);
	logitP = log(P/(1-P))
	logitPpriors = log(Ppriors/(1-Ppriors));
	a = 2;
	P2 = logitPpriors + a*(logitP - logitPpriors);
	P2 = exp(P2)/(1+exp(P2));
	P = P2;
	return(P);

}


.write_h5 <- function(file, m) {
  if(!require(rhdf5))
      stop("could not load package rhdf5")

  h5createFile(file)
  h5write(as.matrix(m),file,"predictions")
  h5write(dimnames(m)[[1]], file,"pid")
  h5write(dimnames(m)[[2]], file,"go")
}
