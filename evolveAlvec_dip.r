evolve.diploid<-function(pop, alspace, inMat, envmts, emph=NULL, startStates, output.genes,reps=50, steps=5, daughters=6, max.popn = 4, cull.frac = 0,prob=c(1,1,1,1,0,1),printMuts=FALSE,beta=500){
	diploid.list<-list()
	
	#ignoring the cull fraction
	num<-max.popn
	
	for (j in 1:steps){
		cat("starting generation number:",j,"\n")
		
		if (length(pop)<num){
			cat("population in generation",j,"is very small, skip mutations\n")
			#just so we don't mutate in the very first generation when the population is very small,
			#this should not be essential, but otherwise we would end up with having fixed mutations throughout the population.
			mutpop<-pop
		}
		else{
			#this mutates the allele vectors in the population. Depending on the values in the vector prob, 
			#some allele vectors may be handed down without mutating.
			newWorld<-mutateGen(pop,alspace,inMat,prob,plotIt=F,printIt=TRUE)
			mutpop<-newWorld$pop
			alspace<-newWorld$alspace
			inMat<-newWorld$inMat
		}
		
		#generate vectors for which allele vectors will be used for each hybridisation.
		#These vectors just contain the indices of the female and male parent allele vector, respectively
		Findex<-sample(length(pop),daughters*num,replace=TRUE)
		Mindex<-sample(length(pop),daughters*num,replace=TRUE)
		
		for (i in 1:(daughters*num)){
			#build diploid networks from each pair of allele vectors. The build.ntwk function makes a network out of any
			#allele vector, regardless of ploidy, so all we need to do is concatenate the two parents.
			diploid.list[[i]]<-build.ntwk(c(mutpop[[Findex[i]]],mutpop[[Mindex[i]]]),alspace,inMat,plotIt=F)
		}
		
		#could change natural selection to not return networks anymore, that isn't actually necessary anymore
		results.of.selection<-natural.diploid.selection(diploid.list,reps,max.popn,cull.frac,output.genes,startStates,envmts=envmts,emph=emph)
		keepers<-results.of.selection$keep
		
		#generate a new population by homologously recombining the parents of the selected networks  
		count<-0
		for (k in keepers){			
			for (l in 1:daughters){
				pop[[(count)*daughters+l]]<-Recmat(mutpop[[Findex[k]]],mutpop[[Mindex[k]]],alspace)
			}
		count<-count+1
		}
	}
	
	
	return(list(pop=pop,alspace=alspace,inMat=inMat))
}

#in this function we want to make sure that homologous recombination and natural selection only happens within separate populations and not between,
#This function now goes through 1.Mutate, 2.Combine & reproduce & hand down 3.Select (repeat) where as the function above does 1.Mutate 2.Select 3.Combine (repeat)

#Because we're using this order of things, this is going to have to be done in pairs. Therefore every list in pops will have size 2*max.popn after the first generation
evolve.sep.pops<-function(pops, alspace, inMat, envmts, 
												emph=NULL, 
												startStates, 
												output.genes,
												reps=50, 
												steps=5, 
												daughters=6, 
												max.popn = 4, 
												cull.frac = 0,
												prob=c(1,1,1,1,0,1),
												printMuts=FALSE,
												printInherit=TRUE,
                        clean=NULL,
            						evolve="MH",
            						old.pheno=NULL,
                        testEpi=FALSE,
                        beta=500){											
	num.pops<-length(pops)
	selectedPheno<-matrix(0,steps,num.pops)	
	wholePheno<-matrix(0,steps,num.pops)
  epis <- matrix(0,2,3)
	#combine all the pops to make the world (called pop)
	pop<-list()
	curr<-1
	popn.size<-length(pops[[1]])#we assume they're all the same size
	for (i in 1:num.pops){
		pop[curr:(curr+popn.size-1)]<-pops[[i]]
		curr<-curr+popn.size
	}
	diploid.list<-list()#will have one for each population
	popAlvecs<-list()
	#ignoring the cull fraction
	num<-max.popn
	
	for (j in 1:steps){
		cat("starting generation number:",j,"\n")
		if (!(is.null(clean))){
      if (j/clean==ceiling(j/clean)){
        cleaned<-cleanWorld(pop,alspace,inMat)
        alspace<-cleaned$alspace
        inMat<-cleaned$inMat
      }
		}
		#1. Mutate
		if (length(pop)<num*num.pops*2){
			cat("population in generation",j,"is very small, skip mutations\n")
			#just so we don't mutate in the very first generation when the population is very small,
			#this should not be essential, but otherwise we would end up with having fixed mutations throughout the population.
			mutpop<-pop
		}
		else{
			#this mutates the allele vectors in the population. Depending on the values in the vector prob, 
			#some allele vectors may be handed down without mutating.
			newWorld<-mutateGen(pop,alspace,inMat,type=NULL,prob,plotIt=F,printIt=TRUE)#need to make sure the mutpop corresponds to the pop so we don't mix populations
			mutpop<-newWorld$pop
			alspace<-newWorld$alspace
			inMat<-newWorld$inMat
		}
		
		#2. Recombine
		#Directly handing down plants might not work because alspace has changed so in order to hand down change params to mutGen so that some things aren't mutated - done
		#num.hand.down<-sample(1:(popn.size/2),1)
		#for (i in 1:num.pops){
		#	to.hand.down<-sample(1:(popn.size/2),num.hand.down)
		#	for (k in 1:num.hand.down){
		#		diploid.list[[i]][[k]]<-build.ntwk(c())
		#	}
		#}
		
		#all the randomness of which plants reproduce is introduced here - we'll just make diploids out of the pairs of things next to each other (arbitrary)
		homo.rec<-list()
		Comb.index<-list()
		for (i in 1:num.pops){
			homo.rec[[i]]<-list()
			Comb.index[[i]]<-sample(1:(popn.size/2),2*daughters*num,replace=TRUE)
			for (k in 1:(2*daughters*num)){
				homo.rec[[i]][[k]]<-Recmat(mutpop[[(i-1)*popn.size+2*Comb.index[[i]][[k]]-1]],mutpop[[(i-1)*popn.size+2*Comb.index[[i]][[k]]]],alspace,printIt=printInherit)
			}
		}
		
		#3. Selection
		#make diploids out of pairs of things that occur next to each other in homo.rec
		#this assumes that all the populations are the same size
		
		for (i in 1:num.pops){
			diploid.list[[i]]<-list()
      popAlvecs[[i]] <- list()
			for (k in 1:(daughters*num)){
				#build diploid networks from each pair of allele vectors. The build.ntwk function makes a network out of any
				#allele vector, regardless of ploidy, so all we need to do is concatenate the two parents.
				#we combine the pair of things that are next to each other in the list
				diploid.list[[i]][[k]]<-build.ntwk(c(homo.rec[[i]][[2*k-1]],homo.rec[[i]][[2*k]]),alspace,inMat,plotIt=F)
			  popAlvecs[[i]][[k]] <- union(homo.rec[[i]][[2*k-1]],homo.rec[[i]][[2*k]])
      }
		
			#could change natural selection to not return networks anymore, that isn't actually necessary anymore
			if (evolve != "MH"){
				results.of.selection<-natural.diploid.selection(diploid.list[[i]],reps,max.popn,cull.frac,output.genes,startStates,envmts=envmts,emph=emph)
			  if (testEpi){
          if (j==steps){
            epis[i,] <- findEpis(popAlvecs[[i]],popAlvecs[[i]][results.of.selection$keep],alspace)
  			  }
			  }
      }
			else if (j==1 & is.null(old.pheno)){
				results.of.selection<-natural.diploid.selection(diploid.list[[i]],reps,max.popn,cull.frac,output.genes,startStates,envmts=envmts,emph=emph)
  		  if (testEpi){
          if (j==steps){
            epis[i,] <- findEpis(popAlvecs[[i]],popAlvecs[[i]][results.of.selection$keep],alspace)
  			  }
			  }			}
			else if (j==1){
				results.of.selection<-natural.selection.MH(diploid.list[[i]],reps,max.popn,cull.frac,output.genes,startStates,envmts=envmts,emph=emph,old.pheno=old.pheno[i])
  		  if (testEpi){
          if (j==steps){
            epis[i,] <- findEpis(popAlvecs[[i]],popAlvecs[[i]][results.of.selection$keep],alspace)
  			  }
			  }			}
			else{
				results.of.selection<-natural.selection.MH(diploid.list[[i]],reps,max.popn,cull.frac,output.genes,startStates,envmts=envmts,emph=emph,old.pheno=selectedPheno[j-1,i])
  		  if (testEpi){
          if (j==steps){
            epis[i,] <- findEpis(popAlvecs[[i]],popAlvecs[[i]][results.of.selection$keep],alspace)
  			  }
			  }			}
			keepers<-results.of.selection$keep
			
			#get the allele vectors (gametes) that make up the selected population
			count<-1
			for (k in keepers){	#length keepers should be max.popn, if not there will probably be a bug somewhere!
				pop[[(i-1)*2*length(keepers)+2*count-1]]<-homo.rec[[i]][[2*k-1]]
				pop[[(i-1)*2*length(keepers)+2*count]]<-homo.rec[[i]][[2*k]]
				count<-count+1
			}
			popn.size<-2*length(keepers) #hopefully this is the same for all i
			pops[[i]]<-pop[((i-1)*popn.size+1):(i*popn.size)]
			selectedPheno[j,i]<-results.of.selection$selectedAvg #each column is a population
			wholePheno[j,i]<-results.of.selection$wholeAvg      
		}
	}
	return(list(pop=pops,mutpop=mutpop,intraCombs=Comb.index,alspace=alspace,inMat=inMat,pheno=selectedPheno,wholePheno=wholePheno,epis1=epis[1,],epis2=epis[2,]))
}

