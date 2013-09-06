#we want to write a function that is able to assess whether our networks show hybrid vigour.
#So, starting from one plant we let a population evolve and adjust to one or a couple of different environments, then we change the weighting of the environments (emph)
#take x individuals from this point and evolve populations from all of them. After every y generations hybridise between two different populations, record differences in phenotype between parents and hybrid (potentially record the hybrids as well)
#hybridising consists of taking two diploids (or rather the corresponding gametes) doing homologous recombination in the diploid and then picking resulting gametes to cross.
#x=pops
#y=gap
#time for adjusting to initial environment = adjust
#I think we should have adjust>>gaps
#eventually, would be good to change startEnvs and newEnvs to different weightings
#if we keep a record of the hybrids/parent groups at each stage then we can look at how quickly the hybrids adapt to new environments compared to the parents
#might be good to write some sort of script that measure the rate of phenotype evolution in order to do the above - this would be a v complicated way of doing things though,
#to start with perhaps just change emph and keep a record of what the average average phtype is after a set number of generations - 1 or 2 should suffice
library("BoolNet")
library("sfsmisc")

source("sf.sw.directed_matrix.r")
source("cleanWorld.r")
source("evolveAlvec_dip.r")
 source("setEnv.r")
 source("natSelect.r")
  source("phenoEmph.r")
   source("calcAttr.r")
source("Recmat.r")
source("mutateGen.r")
 source("killGenemat.r")
 source("dupGenemat.r")
 source("addEdgemat.r")
 source("loseOutmat.r")
 source("loseInmat.r")
source("build.ntwk.r")
 source("netPlotmat.r")
 source("remPH.r")
 source("toNet.r")
#source("newEpis2.R")
source("quantDomEpis.R")

#to do - edit so that it displays pheno from start, and correct timescale
hybrid.vigour<-function(nodes,adjust,pops,gaps,rep,max.popn,prob=c(1,1,1,1,1,1),envmts=NULL,emph=NULL,fname="test",printInherit=FALSE, evolve="MH",testDom=FALSE,savePops=TRUE){
	#rep is the number of hybrid-making stages the program goes through
	if (is.null(emph)){
		startEmph<-c(1,1,1)
		newEmph<-c(1,1,1)
		newnewEmph<-c(1,1,1)
	}
	else{
		startEmph<-emph[[1]]
		newEmph<-emph[[2]]
		newnewEmph<-emph[[3]]
	}
	if (is.null(envmts)){
		envmts<-list(c(0,1,0),c(1,0,0),c(0,0,1))
	}
	envNumber<-length(envmts[[1]])
	reps<-100
	daughters<-4
	
	#fix starting states of length N=nodes
	startStates<-matrix(0,nodes,reps)
	for (i in 1:reps){
		startStates[,i]<-sample(0:1,nodes,replace=T)
	}
	
	#create initial gametes, set output genes etc
	startWorld<-create.allele.ntwk(nodes,2,2,0.5,0.5,envNumber)
	inMat<-startWorld$inMat
	origamete<-startWorld$alvec
	alspace<-startWorld$alspace
	ntwk<-build.ntwk(origamete,alspace,inMat,plotIt=TRUE)
	output.genes<-select.modules(ntwk,envNumber)
	
	#record phenotypes in a nice vector (instead of yucky list)
	intra.selected.pheno<-matrix(0,nrow=adjust+rep*gaps,ncol=pops)
	intra.whole.pheno<-matrix(0,nrow=adjust+rep*gaps,ncol=pops)
	inter.selected.pheno<-vector('numeric',rep)
	inter.whole.pheno<-vector('numeric',rep)
	
	#record how well hybrids cope with change in environment relative to the parents - not implemented yet
	inter_parent<-0
	intra_parent<-0
	
	#a population adjusts to the environment.
	firstWorld<-evolve.sep.pops(pops=list(list(origamete,origamete)), alspace, inMat, envmts, emph=startEmph, startStates, output.genes,
								reps, steps=adjust, daughters=daughters, max.popn = max.popn, cull.frac = 0,prob=prob,printMuts=FALSE,printInherit=printInherit,clean=25,evolve=evolve)
	#take best performing diploid as the start
	firstPlant<-firstWorld$pop[[1]]
	alspace<-firstWorld$alspace
	inMat<-firstWorld$inMat
	parents<-list()
	for (i in 1:pops){
		parents[[i]]<-firstPlant
		intra.selected.pheno[1:adjust,i]<-firstWorld$pheno
		intra.whole.pheno[1:adjust,i]<-firstWorld$wholePheno
	}
		
	#At each stage we need to save the plants, the allele space and inMat and the combinations of hybrids that are made (actual hybrids probably unnecessary). The phenotypes should be both saved and written to file
	#Am saving and outputting these things just so that extra unforseen analysis is possible after the program is run
	#So, populations will have one entry for each stage (i.e. one for startWorld, firstWorld, and then one for each rep)
	#from the third entry onwards (the first rep) we will have $plants (mutpop) $alspace $inMat $intraCombs $interCombs $intraPheno $interPheno in that order.
	#note $plants is a list of length pops*max.popn*2 (for i>2)
	#parents will not be saved but it is a list of length pops of lists of length 2*max.popn which contains the selected plants (have decided not to save this) 
	#to start off with, the parents in firstWorld are pops identical copies of the same plant
	populations<-list()
	#unnecessary output here for testing purposes
	populations[[1]]<-list(plants=list(origamete,origamete),alspace=startWorld$alspace,inMat=startWorld$inMat,sS=startStates,env=envmts,emph=newEmph,mod=output.genes)
	populations[[2]]<-list(plants=parents,alspace=alspace,inMat=inMat)
	
  hybAlvecs<-list()
  hybEpis <- matrix(0,rep,3)
  pop1Epis <- matrix(0,rep,3)
  pop2Epis <- matrix(0,rep,3)
  if (rep >=1){
    hybrid_results<-list()
		for (i in 1:rep){
			cat("starting repetition number",i,"\n")
			#evolve these separate populations for y=gaps generations
			newWorld<-evolve.sep.pops(parents,alspace=alspace,inMat=inMat,envmts=envmts,emph=newEmph,startStates,output.genes,
										reps=reps,steps=gaps,daughters=daughters,max.popn=max.popn,cull.frac=0,prob=prob,printMuts=FALSE,printInherit=printInherit,
										evolve=evolve,old.pheno=intra.selected.pheno[adjust+(i-1)*gaps,],testEpi=F)
			parents<-newWorld$pop #these are the same as the intra hybrids and are also used to carry on evolution
			plants<-newWorld$mutpop #we will make inter hybrids from these
			clean<-cleanWorld(plants,newWorld$alspace,newWorld$inMat)
      alspace<-clean$alspace
			inMat<-clean$inMat
			pop1Epis[i,] <- newWorld$epis1
      pop2Epis[i,] <- newWorld$epis2
			#New way of making hybrids - at each stage 1.Mutate 2.Combine (as in evolution) but instead of selection just take the hybrid phenotype
			#one way of doing this (and still recording what's going on) would be to take the mutpop vector from evolveAlvec, and use this for the hybrids just as we used it for the parents
			#then in order to make direct comparison we would want to select on hybrids. This would mean that parentPheno wouldn't exist - parents would be the intra hybrids (which perhaps makes more sense)
			#for the interCombs we will select random pairs where each element is from ({1,2,..,max.popn},{1,2,...,pops})
			#in think it is good that some of the interCombs will in fact be intra
			#same number of pairs in interCombs as in intraCombs - makes easy printing and means that we select the same proportion for each, though in different ways
			#not sure that selection is really a good idea. Perhaps we should look at population phenotype as well as selected phenotype
			intraCombs<-newWorld$intraCombs 
			#convert so that it can be written to file in way that makes sense
			intra<-matrix(0,4,pops*length(intraCombs[[1]])/2)
			for (j in 1:pops){
				for (k in 1:(length(intraCombs[[1]])/2)){
					intra[,(j-1)*pops+k]<-c(j,intraCombs[[j]][[2*k-1]],j,intraCombs[[j]][[2*k]])
				}
			}
			
			#Make inter hybrids
			#do we really want this many on the first run?
      if (pops==2){
        inter<-rbind(rep(1,daughters*pops*max.popn),
  					sample(1:(length(plants)/(2*pops)),daughters*pops*max.popn,replace=TRUE),
						rep(2,daughters*pops*max.popn),
						sample(1:(length(plants)/(2*pops)),daughters*pops*max.popn,replace=TRUE))
      }
      else {
        firstpop<-sample(1:pops,daughters*pops*max.popn,replace=TRUE)
			  inter<-rbind(firstpop,
						sample(1:(length(plants)/(2*pops)),daughters*pops*max.popn,replace=TRUE),
						sapply(firstpop,function(x) sample((1:pops)[-x],size=1)),
						sample(1:(length(plants)/(2*pops)),daughters*pops*max.popn,replace=TRUE))
      }
      #the interCombs is a little unnecessary here, I think it just makes things a bit clearer. Note length of plants should be 2*max.popn.pops
			interCombs<-vector('numeric',0)
			homo.rec<-list()
			hyb.list<-list()
      hybAlvecs <- list()
			for (j in 1:(daughters*pops*max.popn)){
				interCombs[2*j-1]<-(inter[1,j]-1)*length(plants)/(2*pops)+inter[2,j]
				interCombs[2*j]<-(inter[3,j]-1)*length(plants)/(2*pops)+inter[4,j]
				#to clarify, both the 2j-1 th and 2j th plants will undergo homologous recombination and then the resulting allele vectors from each will be combined to make a diploid:
				homo.rec[[2*j-1]]<-Recmat(plants[[2*interCombs[2*j-1]-1]],plants[[2*interCombs[2*j-1]]],alspace,printIt=printInherit)
				homo.rec[[2*j]]<-Recmat(plants[[2*interCombs[2*j]-1]],plants[[2*interCombs[2*j]]],alspace,printIt=printInherit)
				hyb.list[[j]]<-build.ntwk(c(homo.rec[[2*j-1]],homo.rec[[2*j]]),alspace,inMat,plotIt=F)
        hybAlvecs[[j]] <- union(homo.rec[[2*j-1]],homo.rec[[2*j]])			
      }
			
			#select hybrids and get phenotypes
			intraPheno<-newWorld$pheno #a gaps by pops matrix (selected and population pheno)
			intra.selected.pheno[(adjust+(i-1)*gaps+1):(adjust+i*gaps),]<-intraPheno
			intra.whole.pheno[(adjust+(i-1)*gaps+1):(adjust+i*gaps),]<-newWorld$wholePheno
			if (evolve=="MH"){
				old.pheno<-mean(intra.selected.pheno[adjust+i*gaps-1,])
				#I think this is fair - use same old.pheno as was used for most recent generation of normal plants and average over all populations. Might be better way
				results.of.selection<-natural.selection.MH(hyb.list,reps,pops*max.popn,cull.frac=0,output.genes,startStates,envmts=envmts,emph=newEmph,old.pheno=old.pheno)
        #hybEpis[i,] <- findEpis(hybAlvecs,hybAlvecs[results.of.selection$keep],alspace)			
      }	
			
      else{
				results.of.selection<-natural.diploid.selection(hyb.list,reps,pops*max.popn,cull.frac=0,output.genes,startStates,envmts=envmts,emph=newEmph)
			}	
			interPheno<-c(results.of.selection$wholeAvg,results.of.selection$selectedAvg)
			inter.selected.pheno[i]<-interPheno[2]
			inter.whole.pheno[i]<-interPheno[1]

			#write the crossing information to file	- may well actually move this to later in order to add the phenotypes in each row
			intra<-t(intra)
			inter<-t(inter)	
			crosses<-paste("(",intra[,1],",",intra[,2],")x(",intra[,3],",",intra[,4],")","  (", inter[,1],",",inter[,2],")x(", inter[,3],",",inter[,4],")",sep="")
			filename<-paste(fname, "/hybridresults.txt")
			#write.table(crosses,filename,row.names=FALSE,append=TRUE)
			#write.table(intraPheno,filename,append=TRUE)
			#write.table(interPheno,filename,append=TRUE)
			if (savePops){
        populations[[i+2]]<-list(plants=plants,alspace=alspace,inMat=inMat,intraCombs=intra,interCombs=inter,intraPheno=intraPheno,interPheno=interPheno, interHybs=hyb.list)
			}
      if (testDom){
        hybrid_results[[i]]<-testDomEpi(nodes,pops,alspace,inMat,rep=i,population=populations,max.popn=max.popn,gene.list=output.genes,startStates=startStates)
      }
    }
	}
	#present phenotypes graphically
	#selected.pheno<-inter.selected.pheno/intra.selected.pheno
	#whole.pheno<-inter.whole.pheno/intra.whole.pheno
  	subtitle<-paste(nodes,"nodes,",gaps,"gaps,",rep,"rep, max.popn=",max.popn,"prob=(",paste(prob,collapse=","),")")
  	plot(0,0,col="white",ylim=c(min(intra.whole.pheno),max(intra.whole.pheno,inter.whole.pheno)),xlim=c(1,adjust+rep*gaps),ylab="average phenotype",xlab="generation",main="phenotype evolution",sub=subtitle)
	lines(1:(adjust+rep*gaps),intra.whole.pheno[,1],col=1,type="l")
  	lines(c(adjust,adjust),c(min(intra.whole.pheno),intra.whole.pheno[adjust,1]),lty=2)
  
	for(i in 2:pops){
		lines(1:(adjust+rep*gaps),intra.whole.pheno[,i],col=i,type="l")
	}
	points(adjust+gaps*(1:rep),inter.whole.pheno,col="darkgreen")

  return(list(populations,intra.whole.pheno,inter.whole.pheno,hybrid_results,pop1Epis,pop2Epis,hybEpis))
  

}
