glmnetDalpha = function(L,ProgMap, MAXVAR)
{
        require(glmnet);

		#First find the proteins that both have annotations + also progmap profiles
		common = intersect(names(L), rownames(ProgMap));
		Y = L[common];
		X = ProgMap[common,];
	
		#Fit for the common set of proteins (since for them there is Y and X for the regression fitting step)
        f = cv.glmnet(y = Y, x = X, family = "binomial" , dfmax = MAXVAR, alpha=0.5, maxit=1000, type.measure="auc");
        D = predict(f, ProgMap, s = "lambda.min")
        names(D) = rownames(ProgMap);
		return(D);
}

#BMRF functions for reading files and preparing inputs
loadNet = function(FILE)
{
    # Loads a network from a file in triplet format, i.e:
    # A \t B \t Value
    
	data = as.data.frame(read.table(FILE, sep = "\t"));
	proteins = as.vector(unique(c(as.character(data[,1]), as.character(data[,2]))));	# Set of proteins
	l = length(proteins); # #of proteins
	ids = vector(mode = "numeric", length = length(proteins));
	ids = 1:length(proteins);
	names(ids) = proteins;

	iv = ids[as.vector(data[,1])];
	jv = ids[as.vector(data[,2])];
	xv = as.vector(data[,3]);


	Net = sparseMatrix(i=iv,j=jv,x=xv, dims=c(l,l));
	Net = Net + t(Net);
	rownames(Net) = proteins;
	colnames(Net) = proteins;
	
	return(Net);
	
}

loadAnn = function(FILE)
{

	# Takes an annotation file in triplet format i.e:
	# protein \t annotation 
	# and outputs a binary sparse matrix with the information

	data = as.data.frame(read.table(FILE, sep = "\t"));

	proteins = unique(data[,1]);
	ids = vector(mode = "numeric", length = length(proteins));
	ids = 1:length(proteins);
	names(ids) = proteins;
	
	annotations = unique(data[,2]);
	idsc = vector(mode = "numeric", length = length(annotations));
	idsc = 1:length(annotations);
	names(idsc) = annotations;	
	
	igo = ids[as.vector(data[,1])];
	jgo = idsc[as.vector(data[,2])];
	L = sparseMatrix(i = igo, j = jgo , x = 1, dims = c(length(proteins), length(annotations)));	
	
	rownames(L) = proteins;
	colnames(L) = annotations;
	
	return(L);

}


