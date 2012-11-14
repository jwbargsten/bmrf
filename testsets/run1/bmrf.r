# Functions used by bmrf.r script  


###########################################################
# BMRFz for one network and a domain matrix
# BMRF with network uncertainty
# Author: Yiannis Kourmpetis (yiannis.kourmpetis@wur.nl)
###########################################################
# Input arguments: 
# 1. File with the network topology, with format:
# protA protB confidence network_ID
#
# 2. GO Annotation file:
# protA GOterm evidence
#
# 3. Domain file

library(Matrix);
#library(ROCR);
library(brglm);
library(mvtnorm);	
library(glmnet);
#library(GO.db);
#library(igraph);	


source("bmrf_functions.r");


# Load the data

bmrf = function(FILEnet, FILEann, FILEclust, outprobsfile) 
{
    burnin = 20;
    niter = 20;
    

    # Load, network, annotations and clusters

    A = loadSingleNet(FILEnet); # Loads the data in a single network (so netid is not used at the moment)    
    L = loadAnnFixedProteins(FILEann, A);
    D = loadAnnFixedProteins(FILEclust, A);
    
    # Take the goterms that are neither very general nor very sparse
    minGOsize = 20;
    maxGOsize = (0.9*nrow(A));

    Ls = colSums(L);
    sel = which(Ls >= minGOsize & Ls <= maxGOsize);
    L = L[,sel];
    
    # Do the same for clusters
    minClustsize = 20; 
    maxClustsize = (0.9*nrow(A));

    Ds = colSums(D);
    sel = which(Ds >= minClustsize & Ds <= maxClustsize);
    D = D[,sel];
    
    u = rowSums(L);
    U = which(u == 0);
        
    
    for(i in 1:ncol(L))
	{
		cat(paste(i, " "));			
		#Make the elastic net fitting and predictions:
		Lsingle = L[,i];
		Lsingle[U] = -1;
		
		glmnetpred = try(glmnetDalpha(Y = Lsingle, X = D, MAXVAR = (ncol(D)-1)), silent=FALSE);
		if(class(glmnetpred) == 'try-error')
		{
			next;
		}
        
		# Prepare the Labellings vector
        

		# Run BMRFz 
		posteriors = try(BMRFz(A=A, L = Lsingle, D = glmnetpred, burnin = burnin, niter = niter), silent=FALSE);
		if(class(posteriors) == 'try-error')
		{
			# If BMRFz does not run for any reason, report the EN probs only
		 	next;
		}
		
		posteriors = calibrate(posteriors);
		
		names(posteriors) = rownames(A);
		out = cbind(rownames(A), round(posteriors,5));
		go = colnames(L)[i];

		go = rep(x=colnames(L)[i], times=length(posteriors));
		probs = cbind(names(posteriors), go, posteriors);
		write.table(probs, file= outprobsfile, append=TRUE, quote=F, row.names=F, col.names=F, sep = "\t");			
	}
	
}

    
    
         

