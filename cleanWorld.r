#This function removes alleles that have been completely lost from the population from the allele space and the input Matrix
#this is not necessary, but could become important if evolution runs on for a while and lots of new alleles are generated.
#handling a matrix with thousands of columns and rows (which we might get in a decent sized experiment) might be slow.

cleanWorld<-function(population,alspace,inMat){
	
	#make up a vector of all the alleles in the population, get rid of duplicates
	#the new vector will thus not include all the alleles of the allele space that have been lost from the population
	keep<-unique(unlist(population))
	cat("cleanup will remove",length(alspace$name)-length(keep),"alleles from the allele space of length", length(alspace$name),"\n")
	#find out which indices in the original allele space these alleles correspond to
	keepIndices<-which(alspace$name %in% keep)

	#generate a new allele space from the selected alleles, cut down the input matrix accordingly
	newspace<-data.frame(name=vector('character',length(keep)),exp=vector('character',length(keep)))
	newspace$name<-alspace$name[keepIndices]
	newspace$exp<-alspace$exp[keepIndices]
	newspace$group<-alspace$group[keepIndices]
	
	newMat<-inMat[keepIndices,keepIndices]
	
	return(list(alspace=newspace,inMat=newMat))
}