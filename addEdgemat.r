#This function handles allele vectors. It generates a copies an allele, and adds another input to the new version's input expression.
#this can be done by 4 modes (gates): AND, AND NOT, OR and OR NOT.
addEdgemat<-function(alvec,alspace,inMat,targetAllele=NULL,inputGene=NULL,gate=NULL,printIt=TRUE){
	alnames<-as.character(alspace$name)
	#pick random node if no target node is specified, choose a node that is not an environment node
	if (is.null(targetAllele)){
		targetAllele<-sample(alvec[-grep("e",alvec)],1)
	}

   	#choose a gene to make a new imput from
	if (is.null(inputGene)){
    	alsplit<-strsplit(alvec,split="_")
	  	genes<-unique(sapply(alsplit, function(vec) vec[1]))
		inputGene<-sample(genes,1)
	}
  	#pick a logical gate (AND, AND NOT, OR, OR NOT)
	if (is.null(gate)){	
		gate<-sample(c("&","& !","|","| !"),1)
	}
  	#load input expression
	inExp<-alspace$exp[alnames==targetAllele]
	#get factors (individual input genes) and bracket groups (e.g. "( e1 & ! g4 )" )
  	bgroups<-union(unlist(get.bracket.groups(inExp)),inExp)
	expFacs<-get.factors(inExp)
  	#depending on whether a factor or a bracket group is chosen as the parent, 
  	#slightly different algorithms need to be used. If just gsub is used you sometimes get
  	#unwanted behaviour, like the "g4" in "g47" being replace with a new expression

  	if (sample(length(expFacs)+length(bgroups),1)<=length(expFacs)){
    	#pick a random factor
    	parent<-sample(expFacs,1)
		inExpVec<-strsplit(inExp,split=" ")[[1]]
    	#construct the new term and replace the corresponding element in the vectorised form of the inExp
		newTerm<-paste("(",parent,gate,inputGene,")")
    	inExpVec[inExpVec==parent]<-newTerm
    	#make up a new input expression
    	inExp<-paste(inExpVec,collapse=" ")
  	}
  	else {
    	#alternatively, pick a random bracket group and replace it in the expression
    	parent<-sample(bgroups,1)
	  	inExpVec <- strsplit(inExp,split=" ")[[1]]
      newTerm<-paste(" (",parent,gate,inputGene,") ")
    	inExpVec[inExpVec==parent]<-newTerm
      inExp<-paste(inExpVec,collapse=" ")
  	}
	#make up the new allele name
	splitNames<-strsplit(alnames,split="_")		
	lastIndex<-max(sapply(splitNames,function(vec) as.integer(vec[2])))
	
	#put the new allele with its changed name and input expression into the allele space
	newName<-paste(strsplit(targetAllele,split="_")[[1]][1],lastIndex+1,sep="_")
	group<-alspace$group[alnames==targetAllele]
	alspace[length(alspace$name)+1,]<-c(newName,inExp,group)
	
	#copy column and row in the input Matrix
	inMat<-rbind(inMat,inMat[targetAllele,])
	inMat<-cbind(inMat,inMat[,targetAllele])
	colnames(inMat)<-as.character(alspace$name)
	rownames(inMat)<-as.character(alspace$name)
	
	#finally turn on the inputs from from all alleles of the new input in this allele vector
	#this is quite arbitrary, other options would be to set inputs from all alleles (use alspace$name instead of alvec)
	#using a random sample (like all, but set values to 0 or 1) or just one allele (pick one)
	alsPresent<-grep(inputGene,alvec,value=T)
	inMat[newName,alsPresent]<-1
	
	#now put the new allele into the allele vector
	alvec[which(alvec==targetAllele)[1]]<-newName
	
	if (printIt){cat("adding edge from ",inputGene," to ",targetAllele," gate(",gate,") the new allele is named ",newName,"\n",sep="")}

	return(list(inMat=inMat,alspace=alspace,alvec=alvec))
}
