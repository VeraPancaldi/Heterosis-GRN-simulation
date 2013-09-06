#now this function is quite important. It builds a new network given an allele vector, the allele space and the input matrix

build.ntwk<-function(alvec,alspace,inMat,plotIt=TRUE){
	#this function will complain if the allele vector contains several expressions for one allele.
	#in this case it goes on to assume it is dealing with a probabalistic network
	#an easy solution is to just ingnore all exact duplicates, can't think of a reason why not.
	alvec<-unique(alvec)

	#we need two more genes in the network, which will just be constitutively on and off. 
	#These can serve as placeholders if an input expression ends up with missing elements after processing or if an 
	#input expression has lost inputs
	netList<-list()
  alsplit<-strsplit(alspace$name,split="_")
  
  #go through the alleles in this allele vector
	for (i in alvec){
    #load the input Expression that corresponds to the current allele
		alIndex<-which(alspace$name==i)
		inExp<-alspace$exp[alIndex]

    #extract the factors (i.e. the genes that input into this allele) from the input Expression
		factors<-unique(get.factors(inExp))
    inExpVec<-strsplit(inExp,split=" ")[[1]]
    
    #go through the factors of this inExp
    for (j in factors){
      #find all the alleles in the world that are of this gene
			alsWorld<-alspace$name[which(sapply(alsplit,match,j,nomatch=0)[1,]>0)]
			#which of these input into the current allele?
      hasInputFrom<-inMat[alsWorld,alIndex]
      isPresent<-alsWorld %in% alvec
			ins<-vector('numeric',0)
			#which alleles of this gene are present in the alvec AND input into the current allele?
      ins<-alsWorld[hasInputFrom & isPresent]
      if (length(ins)==0){ #if none, turn input off
				newfac<-"OFF"
			}
			else if (length(ins)==1){ #if one, take that as input
				newfac<-ins[1]			
			}
			else { #if several, use all of them as inputs, connected by OR gates
				newfac<-paste("(",paste(ins,collapse=" | "),")")
			}     
      #replace the factors in the input expression with the appropriate alleles 
			inExpVec[which(inExpVec==j)]<-newfac
		}
    #and write the new input expression into a nice list
		inExp<-paste(inExpVec,collapse=" ")
		netList[[i]] <- c(i,inExp)
	}
	
  #add placeholder alleles ON and OFF
	netList[[length(alvec)+1]]<-c("OFF", "0")
	netList[[length(alvec)+2]]<-c("ON", "1")
	
  #construct a Boolean network from the list of alleles and correctly formed input expressions
  network <- toNet(netList)
	
	
	#Now that the network is constructed, we need to remove the placeholder nodes ON and OFF
  #otherwise the network will look very strange as ON and OFF will appear as nodes with lots of connections
	N<-length(network$genes)-2
	InputFunction<-vector('numeric',0)
	
	for (g in 1:N){
		changed<-FALSE
    #remove OFFs
		if ((N+1) %in% network$interactions[[g]]$input){
			oldInputIndex<-which(network$interactions[[g]]$input==(N+1))
      		for (ind in oldInputIndex){	
        		InputNumber<-length(network$interactions[[g]]$input)
			  	  InputFunction<-remPH(network$interactions[[g]]$func,ind,InputNumber,0)
			    	network$interactions[[g]]$input<-network$interactions[[g]]$input[-ind]
      		}
			changed<-TRUE
		}
		
		if ((N+2) %in% network$interactions[[g]]$input){
			oldInputIndex<-which(network$interactions[[g]]$input==(N+2))
      		for (ind in oldInputIndex){  
        		InputNumber<-length(network$interactions[[g]]$input)
			    	InputFunction<-remPH(network$interactions[[g]]$func,ind,InputNumber,1)
			    	network$interactions[[g]]$input<-network$interactions[[g]]$input[-ind]
      		}
			changed<-TRUE
		}

		if (length(network$interactions[[g]]$input)==0){
			network$interactions[[g]]$input<-0
		}
		if (changed==TRUE){
			network$interactions[[g]]$func<-InputFunction
		}
	}
	
	#for some weird reason, BoolNet often puts the first environment node at the last index in the network. This is annoying because it means the off and on nodes are not removed properly. <-- This problem should hopefully be solved
	network$genes<-network$genes[1:N]
	network$interactions<-network$interactions[1:N]
	network$fixed<-network$fixed[1:N]

	if(plotIt){
		netPlotmat(network,alvec)
	}

	return(network)
}


get.factors<-function(inExp){
	op <- c("!", "&", "\\|", "\\(", "\\)")
	for (sym in op){
		inExp<-gsub(sym,"",inExp)
	}
	splitted<-strsplit(inExp,split=" ")[[1]]
	factors<-splitted[splitted>0]
	return(factors)
}
