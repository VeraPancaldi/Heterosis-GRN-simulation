#This function handles allele vectors. It generates a new allele by setting a factor in the input expression to ON or OFF
loseInmat<-function(alvec,alspace,inMat,targetAllele=NULL,inputGene=NULL,type=NULL,printIt=TRUE){
	alnames<-as.character(alspace$name)
	#pick random node if no target node is specified, choose a node that is not an environment node
	#it is possible, though highly unlikely, that no gene will have an input, 
	#the function stops after 20 tries to prevent an infinite loop
	envnumber<-length(grep("e",alvec))
	foundOne=T
	count<-1
	if (is.null(targetAllele)){
		repeat {
			targetAllele<-sample(alvec[(envnumber+1):length(alvec)],1)
			count<-count+1
			if (nchar(alspace$exp[targetAllele])>1){
				break
			}
			else if (count>20){
				cat("no appropriate target for input removal found \n")
				foundOne=F
				break
			}
		}
	}
	if (foundOne){
		#load the input expression of the target allele and get the factors that are involved in the input expression
		inExp<-alspace$exp[alnames==targetAllele]
    
    factors<-get.factors(inExp)		
		
		#pick an input gene to be removed and a mode of removal
		if (is.null(inputGene)){
			inputGene<-sample(factors,1)
		}	
		if (is.null(type)){
			type<-sample(c("ON","OFF"),1)
		}
		
		#simply replace the target input with ON or OFF. when the network is built, these will be added as nodes, 
		#but the transition functions are then processed to exclude these inputs. 
		#use stringsplit to make sure the right elements are replaced, not things like g4 in g47
    inExpVec <- strsplit(inExp,split=" ")[[1]]
    inExpVec[which(inExpVec==inputGene)] <- type
    inExp <- paste(inExpVec,collapse=" ")
		
		#finding the last index in the allele space, the ID of the newly formed allele will be that +1
		splitNames<-strsplit(alnames,split="_")		
		lastIndex<-max(sapply(splitNames, function(vec) as.integer(vec[2])))
		
		#the new allele is in the same linkage group as the original, but name and input esxpression are different
		newName<-paste(strsplit(targetAllele,split="_")[[1]][-2],lastIndex+1,sep="_")
		group<-alspace$group[alnames==targetAllele]
		
		#put the new allele into the allele space
		alspace[length(alnames)+1,]<-c(newName,inExp,group)
		
		#copy column and row in the input Matrix
		inMat<-rbind(inMat,inMat[targetAllele,])
		inMat<-cbind(inMat,inMat[,targetAllele])
		colnames(inMat)<-as.character(alspace$name)
		rownames(inMat)<-as.character(alspace$name)
		
		#finally turn off the inputs from all alleles of the input gene that was just lost, 
		#this is mainly for clearing up, it should not even be necessary
		alsWorld<-grep(inputGene,as.character(alspace$name),value=T)
		inMat[alsWorld,newName]<-0
		
		#now put the new allele into the allele vector
		alvec[which(alvec==targetAllele)[1]]<-newName
		
	if (printIt){cat(targetAllele," loses its input from ",inputGene," the new allele is named ",newName,"\n",sep="")}
	}
	return(list(inMat=inMat,alspace=alspace,alvec=alvec))
}
