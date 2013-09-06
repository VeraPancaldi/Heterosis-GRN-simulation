#This function handles allele vectors. It generates a new allele and sets all the inputs of the alleles
#of the lost output gene from the new allele to 0.
loseOutmat<-function(alvec,alspace,inMat,targetAllele=NULL,outputGene=NULL,printIt=TRUE){
	alnames<-as.character(alspace$name)
	#pick random node if no target node is specified, choose a node that is not an environment node
	envnumber<-length(grep("e",alvec))

	if (is.null(targetAllele)){
    skip<-F
    tryCatch({notEnv<-alvec[-grep("e",alvec)]
            },
            error=function(err){
              if (printIt){cat("there are no alleles in this alvec with outputs to each other. Skip mutation\nThe error message is:",err)}
              skip<-T
            },
            finally={
              outNum<-apply(inMat[notEnv,],1,sum)
              targetAllele<-sample(notEnv[which(outNum>0)],1)
            })

	}
	if (!skip){
	#pick a gene that this allele is currently outputting (putting out? - probably not) to.
	outputals<-which(inMat[targetAllele,]==1)	
	outsplit<-strsplit(alnames[outputals],split="_")
  genes<-unique(sapply(outsplit, function(vec) vec[1]))
	
	if (is.null(outputGene)){
		outputGene<-sample(genes,1)
	}
	
	#finding the last index in the allele space, the ID of the newly formed allele will be that +1
	splitNames<-strsplit(alspace$name,split="_")		
	lastIndex<-max(sapply(splitNames,function(vec) as.integer(vec[2])))
	
	#the new allele has the same in input expression and linkage group to the original, only the name is different
	inExp<-alspace$exp[alnames==targetAllele]
	newName<-paste(strsplit(targetAllele,split="_")[[1]][-2],lastIndex+1,sep="_")
	group<-alspace$group[alnames==targetAllele]
	alspace[length(alnames)+1,]<-c(newName,inExp,group)
	
	#copy column and row in the input Matrix
	inMat<-rbind(inMat,inMat[targetAllele,])
	inMat<-cbind(inMat,inMat[,targetAllele])
	#and re-write the column and row names, they will now include the new allele name 
	colnames(inMat)<-as.character(alspace$name)
	rownames(inMat)<-as.character(alspace$name)
	
	#This is the important bit: find all the alleles of the picked output gene and set their inputs from the new allele to 0
	alsWorld<-grep(outputGene,alspace$name,value=T)
	inMat[newName,alsWorld]<-0
	
	#now put the new allele into the allele vector
	alvec[which(alvec==targetAllele)[1]]<-newName
	
	if (printIt){cat(targetAllele," loses its output to ",outputGene," the new allele is named ",newName,"\n",sep="")}
	}
	return(list(inMat=inMat,alspace=alspace,alvec=alvec))
}
