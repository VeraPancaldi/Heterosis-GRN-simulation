#this function handles allele vectors, all it is doing is removing an allele from the vector, so it should be beautifully simple
killGenemat<-function(alvec,targetAllele=NULL,printIt=TRUE){

	#pick random node if no target node is specified, choose a node that is not an environment node
	
	if (is.null(targetAllele)){
		targetAllele<-sample(alvec[-grep("e",alvec)],1)
		
		if (printIt){cat("removing allele ",targetAllele," (index:",which(alvec==targetAllele),")\n",sep="")}
		
		alvec<-alvec[-which(alvec==targetAllele)]

	}
	else {
		cat("you have specified a target allele, if you write me a function, I can remove it!\n")
	}
		
	return(alvec)

}

get.genes<-function(alvec){
	g<-strsplit(alvec,split="_")
	return(as.character(g[1]))}
