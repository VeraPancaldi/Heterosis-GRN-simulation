#this function performs homologous recombination of allele vectors. It requires linkage groups to work
Recmat<-function(Falvec,Malvec,alspace,printIt=TRUE){
	#find the number of linkage groups
	linkNum<-max(as.numeric(alspace$group))
	Ralvec<-vector('character',0)
	
	#go through linkage groups, for every group pick whether it will be inherited from the male or the female
	for (i in 1:linkNum){
		from<-sample(c("M","F"),1)
		if (from=="M"){
			if (printIt==TRUE){
				cat("\n")	
				cat("inherit group",i,"from M, alleles:  ")
			}
			#if linkage group i is inherited from the male, go through the alleles in this vector and copy all
			#alleles of this linkage group over into the new (recombined) allele vector
			for (j in Malvec){
				if (alspace$group[which(alspace$name==j)]==i){
					if (printIt==TRUE){
						cat(j," ")
					}
					Ralvec<-append(Ralvec,j)
				}
			}
	
		}

		if (from=="F"){
			if (printIt==TRUE){
				cat("\n")
				cat("inherit group",i,"from F, alleles:  ")
			}
			for (j in Falvec){
				if (alspace$group[which(alspace$name==j)]==i){
					if (printIt==TRUE){
						cat(j," ")
					}
					Ralvec<-append(Ralvec,j)
				}
			}

		}
	}
	return(Ralvec)
}