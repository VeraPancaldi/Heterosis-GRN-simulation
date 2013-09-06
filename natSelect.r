

#have edited lines 16 and 23  for the allele matrix

natural.diploid.selection<-function(ntwk.list,reps,max.popn,cull.frac,gene.list,startStates,func="avg",envmts,emph=NULL){
	#ntwk.list is expected to be a list of boolean networks
	n<-length(ntwk.list)
	phenotype<-vector('numeric',n)
	if (is.null(emph)){
		emph<-rep(1,length(envmts))
	}

	#current startStates that are input here should have the same length as the original gametes
	for (i in 1:n){
		env.positions<-grep("e",ntwk.list[[i]]$genes)
		#calculate the startStates for this diploid
		new.startStates<-matrix(0,length(ntwk.list[[i]]$genes),reps)
		for (j in 1:length(ntwk.list[[i]]$genes)){
			if (sum(j==env.positions)>0){#set env startStates
				#which envnode is this?
				envnode<-as.integer(strsplit(strsplit(ntwk.list[[i]]$genes[j],split="_")[[1]][1],split="e")[[1]][2])
				new.startStates[j,]<-startStates[envnode,]
			}
			else{
				allele<-as.integer(strsplit(strsplit(ntwk.list[[i]]$genes[j],split="_")[[1]][1],split="g")[[1]][2])#this is what might also be called the gene 'group'
				new.startStates[j,]<-startStates[allele,]
			}
		}
		#convert startStates to list form
		new.startStates<-as.list(as.data.frame(new.startStates))
		phenotype[i]<-calculate.diploid.phenotype(ntwk.list[[i]],reps,new.startStates,gene.list=gene.list,func=func,envmts,emph)
	}
	num<-min(max.popn,ceiling(n*(1-cull.frac)))
	ordered.phtype<-sort(phenotype,decreasing=TRUE)
	threshold<-ordered.phtype[num]
	keep<-which(phenotype>=threshold)[1:num]
	new.ntwk.list<-vector('list',num)
	for (i in 1:num){
		new.ntwk.list[[i]]<-ntwk.list[[keep[i]]]
	}
	#cat("average phenotype of this generation:   ",mean(phenotype),"\n")
	#cat("average phenotype of selected networks: ",mean(phenotype[keep]),"\n")
	return(list(phenotypes=phenotype,ntwks=new.ntwk.list,keep=keep,selectedAvg=mean(phenotype[keep]),wholeAvg=mean(phenotype)))
}

natural.selection.MH<-function(ntwk.list,reps,max.popn,cull.frac,gene.list,startStates,func="avg",envmts,emph=NULL,old.pheno){#not suitable for diploids atm
	#ntwk.list is expected to be a list of boolean networks
	n<-length(ntwk.list)
	phenotype<-vector('numeric',n)
	if (is.null(emph)){
		emph<-rep(1,length(envmts))
	}
	#current startStates that are input here should have the same length as the original gametes
	for (i in 1:n){
		env.positions<-grep("e",ntwk.list[[i]]$genes)
		#calculate the startStates for this diploid
		new.startStates<-matrix(0,length(ntwk.list[[i]]$genes),reps)
		for (j in 1:length(ntwk.list[[i]]$genes)){
			if (sum(j==env.positions)>0){#set env startStates
				#which envnode is this?
				envnode<-as.integer(strsplit(strsplit(ntwk.list[[i]]$genes[j],split="_")[[1]][1],split="e")[[1]][2])
				new.startStates[j,]<-startStates[envnode,]
			}
			else{
				allele<-as.integer(strsplit(strsplit(ntwk.list[[i]]$genes[j],split="_")[[1]][1],split="g")[[1]][2])#this is what might also be called the gene 'group'
				new.startStates[j,]<-startStates[allele,]
			}
		}
		#convert startStates to list form
		new.startStates<-as.list(as.data.frame(new.startStates))
		phenotype[i]<-calculate.diploid.phenotype(ntwk.list[[i]],reps,new.startStates,gene.list=gene.list,func=func,envmts,emph)
	}

	num<-min(max.popn,ceiling(n*(1-cull.frac)))
	ordered.phtype<-sort(phenotype,decreasing=TRUE)
	#we find out the average selected phenotype of the previous generation. 
	#then pick a network, accept it with exponential probability depending on the difference between it's pheno type and the old average - definitely select if pheno is better.
	#keep doing this until we've selected max.popn networks to keep.
	#Trying different things out atm
	kept<-0
	keep<-vector('numeric',num)
	new.ntwk.list<-vector('list',num)
	while (kept<num){
		rand<-runif(1)
		new.ntwk.index<-sample(1:length(phenotype),1)
		new.pheno<-phenotype[new.ntwk.index]
		if (exp(1000*(new.pheno-old.pheno))>rand){
			keep[kept+1]<-new.ntwk.index
			new.ntwk.list[[kept+1]]<-ntwk.list[[new.ntwk.index]]
			kept<-kept+1
		}
		else{
			keep[kept+1]<-which(phenotype==ordered.phtype[kept+1])[1]
			new.ntwk.list[[kept+1]]<-ntwk.list[[keep[kept+1]]]	
			kept<-kept+1
		}
	}
	#cat("average phenotype of this generation:   ",mean(phenotype),"\n")
	#cat("average phenotype of selected networks: ",mean(phenotype[keep]),"\n")
	return(list(phenotypes=phenotype,ntwks=new.ntwk.list,keep=keep,selectedPhenos=phenotype[keep],selectedAvg=mean(phenotype[keep]),wholeAvg=mean(phenotype)))
}