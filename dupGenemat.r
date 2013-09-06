#this function handles allele vectors. It duplicates an allele in the vector, the input matrix does not need to be changed,
#no input expressions are affected. Should be a breeze :-)
dupGenemat<-function(alvec,parentAllele=NULL,printIt=TRUE){

	#pick random node if no target node is specified, choose a node that is not an environment node	
	if (is.null(parentAllele)){
		targetAllele<-sample(alvec[-grep("e",alvec)],1)
		
		if (printIt){cat("duplicating allele",targetAllele,"(index: ",which(alvec==targetAllele),")\n",sep=" ")}
		
		alvec<-append(alvec,targetAllele)

	}
	else {
		cat("you have specified a target allele, if you write me a function, I can duplicate it!\n")
	}
		
	return(alvec)

}
