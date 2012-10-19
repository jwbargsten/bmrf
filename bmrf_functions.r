# Functions used by bmrf.r script  



BMRFz = function(A, L, D, burnin, niter)
{
    ###########################################################
    # Input arguments:
    # A: The LOWER TRIANGULAR adjacency matrix of the network
    # D: The N by D binary domains matrix
    # L: the labellings vector (0,1,-1)
    # burnin, niter: MCMC iterations 
    ###########################################################

	#---------------------------------------------------------	
	# Load dependencies
	#---------------------------------------------------------
	
	require(brglm);
	require(mvtnorm);	
	require(Matrix);

	A = A + t(A);	# In case it is not symmetric, symmetrize it
	A@x[] = 1; # In case it is not binary, binarize it.


	# Set some internal parameters	
	initZlength = 30; 
#	
#	#---------------------------------------------------------
#	# Construct the schedules for reporting things
#	# genrep.schd: report the iteration
#	# scaleupd.schd: update the scale parameter for the proposal distribution
#	# simstor.schd: store the values 
#	#---------------------------------------------------------
#	
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

loadSingleNet = function(FILE)
{

    # FIle format: 
    # protA protB conf netId
	data = read.table(FILE, sep = "\t");
    proteins  = c(as.character(data[,1]), as.character(data[,2]));
    proteins = unique(proteins);
    proteins = sort(proteins);
    ids = 1:length(proteins)
    names(ids) = proteins
    protA = as.character(data[,1]);
    protB = as.character(data[,2]);
    l = length(proteins);
	ids = 1:l;
	names(ids) = proteins;
    #print(proteins)
    #print(ids)
    iv = ids[protA];
    jv = ids[protB];
    #print(iv)
    #print(jv)
	#Interactions
	Net = sparseMatrix(i=iv,j=jv,x=1, dims=c(l,l));
	Net = Net + t(Net);
	rownames(Net) = proteins;
	colnames(Net) = proteins;
	print("loadSingleNet... ok");
	return(Net);
	
}

loadAnn = function(FILE)
{

	data = read.table(FILE, sep = "\t");
    proteins = data[,1];
    proteins = unique(proteins);
    proteins = sort(proteins);
    
    ids = 1:length(proteins)
    names(ids) = proteins;
    
    ann = as.character(data[,2]);
    ann = unique(ann);
    ann = sort(ann);

    ida = 1:length(ann);
    names(ida) = ann;    

    p = as.character(data[,1]);
    a = as.character(data[,2]);
    
    iv = ids[p];
    jv = ida[a];
    
	L = sparseMatrix(i = iv, j = jv , x = 1, dims = c(length(proteins), length(ann)));	
	L@x[] = 1;
	
	rownames(L) = proteins;
	colnames(L) = ann;
    print("loadAnn... ok");	
	return(L);

}

loadAnnFixedProteins = function(FILE, A)
{

	data = read.table(FILE, sep = "\t");
    proteins = rownames(A)
    #print("proteins")
    #print(proteins)
    
    ids = 1:length(proteins)
    names(ids) = proteins;
    
    ann = as.character(data[,2]);
    ann = unique(ann);
    ann = sort(ann);

    ida = 1:length(ann);
    names(ida) = ann;    

    p = as.character(data[,1]);
    a = as.character(data[,2]);
    #print("p")
    #print(p)
    #print("a")
    #print(a)
    
    iv = ids[p];
    jv = ida[a];
    #print("iv")
    #print(iv)
    #print("jv")
    #print(jv)
    
	L = sparseMatrix(i = iv, j = jv , x = 1, dims = c(length(proteins), length(ann)));	
	L@x[] = 1;
	
	rownames(L) = proteins;
	colnames(L) = ann;
    print("loadAnn... ok");	
	return(L);

}



glmnetDalpha = function(Y,X, MAXVAR)
{
        require(glmnet);
        U = which(Y == -1);
        
        Yf = Y[-U];
        Xf = X[-U,];
		#Fit for the common set of proteins (since for them there is Y and X for the regression fitting step)
        f = cv.glmnet(y = Yf, x = Xf, family = "binomial" , dfmax = MAXVAR, alpha=0.5, maxit=1000, type.measure="auc");
        Yhat = predict(f, X, s = "lambda.min")
        names(Yhat) = rownames(X);
		return(Yhat);
}


calibrate = function(P)
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



