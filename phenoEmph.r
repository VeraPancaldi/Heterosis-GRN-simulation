#have only changed line 25, not sure rest of it needs changing for the allele matrix

calculate.diploid.phenotype<-function(hybrid,reps,startStates=NULL,gene.list,func="avg",envmts,emph){
	#this is a phenotype function in which several sets of output nodes (as specified in gene.list) 
	#are tested under a list of environmental conditions (envmts).
	#if startStates are given, the phenotypes for all the environments will be calculated with the same startStates, 
	#only the environment nodes will be set according to the environment
	#now editing this function so that we can change the environment emphasis, as opposed to having to completely change the environments
	#emph should be a vector of length length(envmts)
	N<-length(hybrid$genes)
	emph<-emph/sum(emph) #it's now normalised
	envNumber<-length(envmts[[1]])
	
	#we now have a diploid, so we need to work out where the output nodes actually are in the male and the female
	new.gene.list<-gene.list
	#cat("length of gene.list",length(gene.list))
	for (i in 1:length(gene.list)){
		module<-gene.list[[i]]
		new.module<-vector('numeric')
		index<-1
		for (j in 1:N){
			letter<-strsplit(hybrid$genes[j],split="")[[1]][1]#this should be "e" or "g", but beware of caplitals!
			gene.homologue<-as.integer(strsplit(strsplit(hybrid$genes[j],split="_")[[1]][1],split=letter)[[1]][2])
			if (any(gene.homologue==module)){
				new.module[index]<-j
				index<-index+1
			}
		}
		new.gene.list[[i]]<-new.module
		#print(hybrid$genes[new.gene.list[[i]]])
	}
	
	#envPheno contains phenotype values for all the environments
	envPheno<-vector("numeric",length(envmts))
	for (envIndex in 1:length(envmts)){
		env<-envmts[[envIndex]]
		attr<-calculate.diploid.attractors(hybrid,reps,startStates,env)
		#attrPheno contains phenotype values for all the attractors of the network in the current env 	
		attrPheno<-vector('numeric',length(attr$attr$attractors))
		#calculate sum of basin sizes of all the found attractors (for weighting)
		allBasins<-0
		for (a in 1:length(attr$attr$attractors)){
			allBasins<-allBasins+attr$attr$attractors[[a]]$basinSize}

		for (i in 1:length(attr$attr$attractors)){
			attractor<-attr$attr$attractors[[i]]
			#statePheno contains phenotype values for all the states in the current attractor  
			statePheno<-vector('numeric',length(attractor$involvedStates))
			
			for (j in 1:dim(attractor$involvedStates)[2]){
				#get status of all nodes in this state
				#not as simple as it seems - we have a matrix and some of the integers are negative
				state<-vector('numeric')
				for (k in 1:length(attractor$involvedStates[,j])){
					if (attractor$involvedStates[k,j]>=0){
						state.section<-digitsBase(attractor$involvedStates[k,j],base=2,N)
					}
					else{
					#see documentation for why we get negative states
						state.section<-digitsBase(attractor$involvedStates[k,j] + 2^32,base=2,N)
					}
						#unfortunately this vector is back to front, so we should reverse it
					state<-c(state,rev(state.section))
				}
				if (func=="avg"){
					#this calculates the phenotype value of the current state and writes it into the statePheno vector
					otherEnv<-0
					thisEnv<-0
					for (e in 1:length(envmts)){
						if(length(new.gene.list[[e]])!=0){
							if (all(envmts[[e]]==env)){
								thisEnv<-mean(state[new.gene.list[[e]]])
							}
							else {
								otherEnv<-otherEnv+mean(state[new.gene.list[[e]]])
							}
						}
					}
					statePheno[j]<-thisEnv-otherEnv/(length(envmts)-1)
				}
			}
			attrPheno[i]<-mean(statePheno)
			#cat("attrPheno[",i,"]: ",attrPheno[i],"\n")
		}
		#cat("env",envIndex, " attrphenos ", attrPheno, "\n")
		envPheno[envIndex]<-mean(attrPheno*attractor$basinSize/allBasins)
	}
	netPheno<-mean(envPheno*emph)
	#cat("envPheno",envPheno,"\n")
	return(netPheno)
}
