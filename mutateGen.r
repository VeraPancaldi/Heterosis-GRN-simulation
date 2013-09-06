#this function handles the input matrix, it deals with a population of allele vectors and performs one action on each of them
mutateGen<-function(pop,alspace,inMat,type=NULL,prob=c(1,1,1,1,1,1),plotIt=FALSE,printIt=F){
	#load mutation functions, might need to change working directory
	probArray<-list(mut=c("loseIn","loseOut","dupGene","killGene","addEdge","none"),prob=prob)
	
	for (i in 1:length(pop)){	
		#if no mutation type specified, pick one at random
		if (is.null(type)){
			r<-sample(1:sum(probArray$prob),1)}
		else {
			r<-type}
		
		#cat("alvec ",i,": ",sep="")
		#perform the selected mutation
		if	(r>sum(probArray$prob[1:5])){						#no mutation
		#	if (printIt){cat("no mutation, the allele vector is passed on unchanged.\n")} disable for now
		}
			
		else if	(r>sum(probArray$prob[1:4])){					#addEdge
			if (printIt){cat("addEdge:    ")}
			newStuff<-addEdgemat(pop[[i]],alspace,inMat,printIt=printIt)
			pop[[i]]<-newStuff$alvec
			alspace<-newStuff$alspace
			inMat<-newStuff$inMat
		}
			
		else if	(r>sum(probArray$prob[1:3])){					#killGene	
			if (printIt){cat("killGene:   ")}
			pop[[i]]<-killGenemat(pop[[i]],printIt=printIt)
		}
			
		else if	(r>sum(probArray$prob[1:2])){					#dupGene	
			if (printIt){cat("dupGene:    ")}
			pop[[i]]<-dupGenemat(pop[[i]],printIt=printIt)
		}
			
		else if(r>sum(probArray$prob[1])){						#loseOut
			if (printIt){cat("loseOut:    ")}
			newStuff<-loseOutmat(pop[[i]],alspace,inMat,printIt=printIt)
			pop[[i]]<-newStuff$alvec
			alspace<-newStuff$alspace
			inMat<-newStuff$inMat
		}
			
		else {													#loseIn
			if (printIt){cat("loseIn:     ")}
			newStuff<-loseInmat(pop[[i]],alspace,inMat,printIt=printIt)
			pop[[i]]<-newStuff$alvec
			alspace<-newStuff$alspace
			inMat<-newStuff$inMat
		}
		
		
		
		if (plotIt){
			network<-build.ntwk(pop[[i]],alspace,inMat)

			netPlotmat(network,pop[[i]])
		}
	}

return(list(pop=pop,alspace=alspace,inMat=inMat))
}