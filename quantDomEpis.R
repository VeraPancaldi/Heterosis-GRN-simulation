#This function tries to quantify the effects of different modes of heterosis.
#The analysis is done for dominance, overdominance, underdominance, positive and negative epistasis.
#This is opposed to the older functions we used to asses these mechanisms, 
#as these only counted the number of occurences of each mechanism.
testDomEpi = function(nodes,pops,alspace,inMat,envmts=list(c(0,1,0),c(1,0,0),c(0,0,1)),newemph=c(1,1,2),rep,population,max.popn,gene.list,startStates,reps=100){
  #initialise
  hybrid_results<-data.frame()
  difmatlist <- list()
  combs<-population[[rep+2]]$interCombs
  
  #go through each hybrid
  for (i in 1:length(combs[,1])){
    #these two counting variables track the number of heterozygous and homozygous linkage groups in the hybrid
    hemicount<-0
    homocount<-0
    
    #this should now be irrelevant  
    if (combs[i,1]==combs[i,3]){
      mainlabel <- paste("hybrid",i,"taken from same population")
    }
    else {
      mainlabel <- paste("hybrid",i,"taken from different populations")
    }
    
    #load the allele vectors of the grandparental gametes
    MgrandDad <- population[[rep+2]]$plants[[(combs[i,1]-1)*max.popn+2*combs[i,2]]]
    MgrandMa <- population[[rep+2]]$plants[[(combs[i,1]-1)*max.popn+2*combs[i,2]-1]] 
    FgrandDad <- population[[rep+2]]$plants[[(combs[i,3]-1)*max.popn+2*combs[i,4]]]
    FgrandMa <- population[[rep+2]]$plants[[(combs[i,3]-1)*max.popn+2*combs[i,4]-1]]
    #and recombine them into the parental ones. Note: because of the randomness of this process,
    #we won't get the same hybrid networks as were used in the last step of the simulation.
    #the set of hybrids here should however be a non-biased sample, so we will go with it.
    Falvec <- Recmat(FgrandMa,FgrandDad,alspace,F)
    Malvec <- Recmat(MgrandMa,MgrandDad,alspace,F)
    fullHybrid <- union(Falvec,Malvec)
    
    #start off by calculating the phenotype of the hybrid. For this we have to construct the network
    #and generate appropriate start states for the calculation of attractors. 
    #This bit of code is copied from natural.diploid.selection
    network<-build.ntwk(fullHybrid,alspace,inMat,plotIt=F)
    env.positions<-grep("e",network$genes)
    #calculate the startStates for this diploid
    new.startStates<-matrix(0,length(network$genes),reps)
    for (j in 1:length(network$genes)){
      if (j %in% env.positions){#set env startStates
        #which envnode is this?
        envnode<-as.integer(strsplit(strsplit(network$genes[j],split="_")[[1]][1],split="e")[[1]][2])
        new.startStates[j,]<-startStates[envnode,]
      }
      else{
        allele<-as.integer(strsplit(strsplit(network$genes[j],split="_")[[1]][1],split="g")[[1]][2])#this is what might also be called the gene 'group'
        new.startStates[j,]<-startStates[allele,]
      }
    }
    #convert startStates to list form
    new.startStates<-as.list(as.data.frame(new.startStates))
    fullPheno<-round(calculate.diploid.phenotype(network,reps,new.startStates,gene.list,envmts=envmts,emph=newemph),10)
    
    #initialise matrix
    hemiPhenos<-matrix(nrow=nodes,ncol=2)
    doublePhenos <- matrix(nrow=nodes,ncol=nodes)
    
    #now go through linkage groups. For every linkage group, generate a hybrid allele vector that inherits this linkage group
    #from only one of its parents
    for (lingroup in 1:nodes) {
      
      knockouts<-alspace$name[which(alspace$group==lingroup)]   #get all the alleles in this linkage group
      lostFromFalvec<-union(Malvec,Falvec[-which(Falvec %in% knockouts)])
      lostFromMalvec<-union(Falvec,Malvec[-which(Malvec %in% knockouts)])     
      depleted<-union(Falvec[-which(Falvec %in% knockouts)],Malvec[-which(Malvec %in% knockouts)])  #this is a hybrid allele vector that doesn't have the linkage group at all. It is used to detect hemizygotes
      
      if (all(lostFromFalvec %in% depleted) | all(lostFromMalvec %in% depleted)){ #if the hybrid is hemizygotic for this linkage group
        hemiPhenos[lingroup,]<-c(NA,NA)
        hemicount=hemicount+1
      }
      else if (all(lostFromFalvec %in% lostFromMalvec) & length(lostFromFalvec)==length(lostFromFalvec)){ #if the hybrid is homozygotic for this linkage group
        hemiPhenos[lingroup,]<-c(NA,NA)
        homocount<-homocount+1      
      }
      else { #if the hybrid is heterozygotic for this linkage group. Most of this code appears twice because it is executed for the networks that have lost the paternal or the maternal linkage group
        if (all(fullHybrid %in% lostFromFalvec)){ #if the hybrid is not actually affected by this loss, e.g. if all alleles that are excluded are inherited from the other parent, along with others (still heterozygotic)
          hemiPhenos[lingroup,1]<-fullPheno
        }
        else {
          network<-build.ntwk(lostFromFalvec,alspace,inMat,plotIt=F)
          env.positions<-grep("e",network$genes)
          #calculate the startStates for this diploid
          new.startStates<-matrix(0,length(network$genes),reps)
          for (j in 1:length(network$genes)){
            if (j %in% env.positions){#set env startStates
              #which envnode is this?
              envnode<-as.integer(strsplit(strsplit(network$genes[j],split="_")[[1]][1],split="e")[[1]][2])
              new.startStates[j,]<-startStates[envnode,]
            }
            else{
              allele<-as.integer(strsplit(strsplit(network$genes[j],split="_")[[1]][1],split="g")[[1]][2])#this is what might also be called the gene 'group'
              new.startStates[j,]<-startStates[allele,]
            }
          }
          #convert startStates to list form
          new.startStates<-as.list(as.data.frame(new.startStates))
          hemiPhenos[lingroup,1]<-round(calculate.diploid.phenotype(network,reps,new.startStates,gene.list,envmts=envmts,emph=newemph),10)
          
          
          #now this bit is for the calculation of epistasis. We calculate the fitness of each double mutant in the network.
          #i.e. we run through each linkage group in order (the loop we are already in). For each lingroup we then also remove each
          #group from the other parent in turn. The network for each double mutant is generated and the fitness is calculated.
          #That matrix is then compared to the matrix of predicted fitnesses (the fitnesses if the effects of both mutations were multiplicative)
          for (othergroup in 1:nodes){
            
            otherknockouts<-alspace$name[which(alspace$group==othergroup)]   #get all the alleles in this linkage group
            if (length(Malvec[-which(Malvec %in% otherknockouts)])>0 & length(Falvec[-which(Falvec %in% knockouts)])>0){
              lostpair<-union(Malvec[-which(Malvec %in% otherknockouts)],Falvec[-which(Falvec %in% knockouts)])
              network<-build.ntwk(lostpair,alspace,inMat,plotIt=F)
              env.positions<-grep("e",network$genes)
              #calculate the startStates for this diploid
              new.startStates<-matrix(0,length(network$genes),reps)
              for (j in 1:length(network$genes)){
                if (j %in% env.positions){#set env startStates
                  #which envnode is this?
                  envnode<-as.integer(strsplit(strsplit(network$genes[j],split="_")[[1]][1],split="e")[[1]][2])
                  new.startStates[j,]<-startStates[envnode,]
                }
                else{
                  allele<-as.integer(strsplit(strsplit(network$genes[j],split="_")[[1]][1],split="g")[[1]][2])#this is what might also be called the gene 'group'
                  new.startStates[j,]<-startStates[allele,]
                }
              }
              #convert startStates to list form
              new.startStates<-as.list(as.data.frame(new.startStates))
              doublePhenos[lingroup,othergroup]<-calculate.diploid.phenotype(network,reps,new.startStates,gene.list,envmts=envmts,emph=newemph)
            }
          }
        }
        
        
        if (all(fullHybrid %in% lostFromMalvec)){ #same again
          hemiPhenos[lingroup,2]<-fullPheno
        }
        else {
          network<-build.ntwk(lostFromMalvec,alspace,inMat,plotIt=F)
          env.positions<-grep("e",network$genes)
          #calculate the startStates for this diploid
          new.startStates<-matrix(0,length(network$genes),reps)
          for (j in 1:length(network$genes)){
            if (j %in% env.positions){#set env startStates
              #which envnode is this?
              envnode<-as.integer(strsplit(strsplit(network$genes[j],split="_")[[1]][1],split="e")[[1]][2])
              new.startStates[j,]<-startStates[envnode,]
            }
            else{
              allele<-as.integer(strsplit(strsplit(network$genes[j],split="_")[[1]][1],split="g")[[1]][2])#this is what might also be called the gene 'group'
              new.startStates[j,]<-startStates[allele,]
            }
          }
          #convert startStates to list form
          new.startStates<-as.list(as.data.frame(new.startStates))
          hemiPhenos[lingroup,2]<-round(calculate.diploid.phenotype(network,reps,new.startStates,gene.list,envmts=envmts,emph=newemph),10)
        }
      }
    }
    #now produce a matrix of predicted phenotypes if there were no epistasis whatsoever. 
    #In this case (theoretically) the fitness of the double mutant should be the product of the fitnesses of the two single mutants,
    #provided the fitness of the non-mutated network has a fitness of 1.
    #thus we calculate by:  mut1pheno/fullPheno * mut2pheno/fullPheno * fullPheno
    #this simplifies to:
    noEpiPhenos <- matrix(nrow=nodes,ncol=nodes)
    for (m in 1:nodes){
      noEpiPhenos[,m] <- hemiPhenos[m,1]*hemiPhenos[,2]/fullPheno
      noEpiPhenos[m,m] <- NA
    }
    #The comparison between the two matrices is a bit hairy. The prediction of fitnesses by multiplying the effects of both mutations is quite idealised.
    #Since there is bound to be some variation, a simple smaller/greater between the two fitnesses might give spurious results. 
    #You could even argue that in a network like ours, there will always be some degree of interaction between two alleles, but when do we start counting it as epistasis?
    #Might do by defining a threshold, say 0.5% of full network fitness. Then we could use the same threshold for the other mechanisms.
    
    #these three lines compare the matrices of calculated and predicted double mutant fitnesses
    posepistable <- doublePhenos>(noEpiPhenos+0.005*fullPheno)
    posepistable[which(is.na(posepistable))] <- 0
    posepis <- sum(posepistable)
    
    #same again for negative epistasis, i.e. double mutant performing better than predicted
    #negepistable <- doublePhenos>(noEpiPhenos+0.005*fullPheno)
    #negepistable[which(is.na(negepistable))] <- 0
    #negepis <- sum(negepistable)
    
    #split negative epistasis up into two factors: epistatic complementation and genetic incompatibilty. 
    #The first should be irrelevant for heterosis, the other is directly opposed and probably causes the fitness collapse
    underdommatrix <- matrix(0,ncol=length(hemiPhenos[,1]),nrow=length(hemiPhenos[,2]))
    underdommatrix[which(hemiPhenos[,1]>fullPheno),] <- 1
    underdommatrix[,which(hemiPhenos[,2]>fullPheno)] <- 1
    
    
    completable<-(doublePhenos<(noEpiPhenos-0.005*fullPheno))*(!(underdommatrix))
    completable[which(is.na(completable))] <- 0
    epicomps <- sum(completable)
    
    incomptable<-doublePhenos<(noEpiPhenos-0.005*fullPheno)*(underdommatrix)
    incomptable[which(is.na(incomptable))] <- 0
    epincomp <- sum(incomptable)
#if(any(underdommatrix==1)){    
 #   print("fullPheno")
#    print(fullPheno)
#    print("hemiphenos")
#    print(hemiPhenos)
#    print("doublephenos")
#    print(doublePhenos)
#    print("underdommatrix")
#    print(underdommatrix)}

    
    #------------------#
    #quantitative calculation of effects of heterosis mechanisms: for dominance and overdominance, 
    #take difference between whole network fitness and average of single mutant fitness.
    #But DOES THAT MAKE SENSE??? Is there a bias in favor of overdominance, since you only count half the difference
    #between the two mutants in the case of dominance? I think not, it's probably a fair comparison
    #underdominance should work just the same.
    #But how can we do the same for epistasis???
    #maybe we can look at each occurrence of epistasis, but for each allele involved in an epistatic interaction 
    #only count the interaction where the effect of epistasis is the strongest. This is to exclude cases where one allele
    #is epistatic to an entire 'pathway', so the epistasis would end up being counted for each pair, giving spurious replication
    #but this would be fairly complicated.
    #------------------#
    
    
    minPhenos<-round(na.omit(mapply(min,hemiPhenos[,1],hemiPhenos[,2])),10)
    maxPhenos<-round(na.omit(mapply(max,hemiPhenos[,1],hemiPhenos[,2])),10)
    #calculate the number of occurences of different effects that removing either allele can have:
    #dominance: removing one of the alleles doesn't cause a change, but removing the other leads to a reduction in the phenotype
    #overdominance: removing either allele reduces the phenotype
    #underdominance: removing one allele increases the phenotype
    #irrelevance: removing either allele does not affect the phenotype
    doms<-sum((maxPhenos>=fullPheno)*(minPhenos<(fullPheno)))
    overdoms<-sum((minPhenos<(fullPheno))*(maxPhenos<(fullPheno)))
    underdoms<-sum(minPhenos>(fullPheno))
    irrels<-sum((minPhenos==fullPheno)*(maxPhenos==fullPheno))
    negligs <- nodes-homocount-hemicount-doms-overdoms-underdoms-irrels
    
    #boxplot(list(na.omit(minPhenos),na.omit(maxPhenos),fullPheno),main=mainlabel)
    #cat(mainlabel,": \n")
    #cat(hemicount,"genes are hemizygous\n")
    #cat(homocount,"genes are homozygous\n")
    
    #quantitative calculation of effects. We don't need the thresholds anymore, 
    #since loci with tiny effects won't make much of a difference to the result.
    #by this point, all the fitness values will have gone through a round(...,10), 
    #so numerical error should not be an issue anymore.
    domeffect <- sum(  ((maxPhenos>=fullPheno)*(minPhenos<fullPheno))*  (fullPheno-(maxPhenos+minPhenos)/2)/fullPheno)
    overdomeffect <- sum(  ((minPhenos<fullPheno)*(maxPhenos<fullPheno))*(fullPheno-(maxPhenos+minPhenos)/2)/fullPheno)
    underdomeffect<- sum(  (minPhenos>fullPheno)*(fullPheno-(maxPhenos+minPhenos)/2)/fullPheno)

    #epistasis: take actual double mutant fitness minus predicted.
    #list differences in order of size. Accept the greatest value.
    #look at next greatest value, accept only if neither of the alleles is in an epistatic pair that has already been accepted.
    #repeat for all other cases of epistasis

    #or: make a triangular matrix of differences, select biggest one in each column --> that sounds like it would cause bias
    #or: start off with a symmetric matrix of differences and select the biggest one from each column
    #then blank (e.g. set to 0) the value in the reverse field, so the same epistatic interaction does not get selected twice.
    #^^that may be the best option and wouldn't be a nightmare to implement
    posepifect <- 0
    negepifect <- 0
    incompepifect <- 0

    #fullPheno: fitness of entire hybrid network, with no genes taken out
    #noEpiPhenos: n x n matrix of all the predicted fitnesses of the double mutants based on the fitnesses of the single mutants.
    #   calculated as mut1fit*mut2fit/fullPheno
    #incompatibility: define as cases in which at least one of the single mutants has a higher fitness than the doubles AND doublePhenos < noEpiPhenos
    
    difmat <- doublePhenos - noEpiPhenos
    difmatsave <- difmat
    
    for (k in 1:ncol(doublePhenos)){
      if (!(all(is.na(difmat[,k])))){
      #cat("col ",k,":   ")
      if (max(na.omit(difmat[,k]))>0){
        posepifect <- posepifect + max(na.omit(difmat[,k]))/fullPheno
        j <- which(difmat[,k]==max(na.omit(difmat[,k])))[1]
#        cat("posepifect: ",posepifect,"  blanking: ",difmat[k,j]," at [",k,",",j,"] and ",difmat[j,k]," at [",j,",",k,"]\n")
        difmat[k,j] <- 0
        difmat[j,k] <- 0
      }
      if (min(na.omit(difmat[,k]))<0){
        j <- which(difmat[,k]==min(na.omit(difmat[,k])))[1]
        if (hemiPhenos[k,1]>fullPheno | hemiPhenos[k,2]>fullPheno){
          incompepifect <- incompepifect + min(na.omit(difmat[,k]))/fullPheno
          }
        else {
          negepifect <- negepifect + min(na.omit(difmat[,k]))/fullPheno
          }
#       cat("negepifect: ",negepifect,"  blanking: ",difmat[k,j]," at [",k,",",j,"] and ",difmat[j,k]," at [",j,",",k,"]\n")
        difmat[k,j] <- 0
        difmat[j,k] <- 0
      }
      }
    }
  hybrid_results<-rbind(hybrid_results,c(fullPheno,hemicount,homocount,doms,overdoms,underdoms,irrels,negligs,posepis,epicomps,epincomp,domeffect,overdomeffect,underdomeffect,posepifect,negepifect,incompepifect))
#  cat("saving into difmatlist[[",i,"]]\n")
    difmatlist[[i]] <- difmatsave
}
  colnames(hybrid_results)<-c("fullPheno","hemis","homos","doms","overdoms","underdoms","irrels","negligs","posepis","epicomps","epincomp","domeffect","overdomeffect","underdomeffect","posepifect","negepifect","incompepifect")
  return(list(hybrid_results,difmatlist,gene.list))
}
