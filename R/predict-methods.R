## predict using a bmrf object
predict.bmrf <- function(b, burnin=20, niter=20, file=NULL, format=c("3col", "long"), verbose=FALSE ) {

  format <- match.arg(format)
  o <- options()
  result <- list()

  ## write header to output file if desired
  if(!is.null(file))
    if(format == "long")
      cat("go", paste0("\t", rownames(b@go)), "\n", sep="", file=file)
    else
      cat("", sep="", file=file)

  for(i in 1:ncol(b@go)) {
    go.proteins <- b@go[,i]
    go.proteins[b@unknown.idcs] = -1

		## make the elastic net fitting and predictions for the functional domains (fd).
    glm.enet.pred <- predict.bmrf_glmnet(bmrf=b, go.idx=i, dfmax = (ncol(b@fd)-1))

    if(is.null(glm.enet.pred))
      next

    ## run BMRFz
    p = predict.bmrf_go(bmrf=b, go.proteins=go.proteins, fd.model=glm.enet.pred, burnin=burnin, niter=niter)
    p = .calibrate(p)

    if(is.null(file)) {
      result[[i]] <- list(go=colnames(b@go)[i], p=p)
    } else {
      switch(format, 
        "3col" = cat(
            paste(names(p), rep(colnames(b@go)[i], length(p)),  p, sep="\t"),
            sep="\n",
            file=file,
            append=TRUE
        ),
        "long" = cat(
            colnames(b@go)[i], 
            paste0("\t", p), 
            "\n", 
            sep="", 
            file=file, 
            append=TRUE
        )
      )
    }
    if(i == 5) break
    if(verbose) cat(".", sep="", file=stderr())
  }
  options(o)

  if(verbose) cat(" finished\n", sep="", file=stderr())

  result
}

predict.bmrf_go <- function(bmrf, go.proteins, fd.model, burnin, niter, initZlength=30) {
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
  D <- fd.model
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
#	Dalpha = glmnetDalpha(L,D, 10)
	Dalpha = D;	
	# Initialize MRF parameters

	M1 = as.vector(A[knowns,knowns] %*% L[knowns]);
	NS = as.vector(dA[knowns]);
	Lk = as.vector(L[knowns])
	Dalphak = as.vector(Dalpha[knowns]);

	regtable = as.data.frame(cbind(Lk, M1,NS, Dalphak));
	colnames(regtable) = c("L", "M1", "NS", "Dalpha");
	
	regtable.fit = brglm(regtable$L ~ regtable$M1 + regtable$NS + regtable$Dalpha,
						family=binomial(link = "logit"), 
						method = "brglm.fit");

	
	Z = rmvnorm(mean = regtable.fit$coefficients, 
							sigma=vcov(regtable.fit), n=initZlength);
							
	MRFparams = as.vector(regtable.fit$coefficients);


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

		t1 = proc.time();
		
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

gibbssample.labels <- function(params, data, num.unknowns) {
		y =  params[1] + data[,1]*params[,2] + data[,2]*params[,3] + data[,3]*params[,4];

		d = (1/(1 + exp(-y)))  - runif(num.unknowns);
		L[unknowns[d >= 0]] = 1;
		L[unknowns[d  < 0]] = 0;	

}

predict.bmrf_glmnet <- function(bmrf, go.idx, dfmax) {
  yf <- bmrf@go[-bmrf@unknown.idcs,go.idx]
  xf <- bmrf@fd[-bmrf@unknown.idcs,]

  #cat(sum(yf), sum(xf), "\n", sep="  ")

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
  names(yhat) = rownames(bmrf@fd);

  if(class(yhat) == 'try-error')
    return(NULL)

  return(yhat);
}

.calibrate = function(P)
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
