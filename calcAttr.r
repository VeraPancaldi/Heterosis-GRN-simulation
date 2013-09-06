calculate.diploid.attractors<-function(hybrid,reps,startStates=NULL,envFunc){
	N<-length(hybrid$genes)
	if (is.null(startStates)){
		startStates<-vector('list',reps)
		for (i in 1:reps){
			startStates[[i]]<-sample(0:1,N,replace=T)
		}
	}
	hybrid<-setEnv(hybrid,envFunc)
	for (i in 1:reps){		
		#set environment nodes to the current environment values
		startStates[[i]][1:length(envFunc)]<-envFunc}	
	
	fixed.on<-which(hybrid$fixed==1)
	fixed.off<-which(hybrid$fixed==0)
	for (i in 1:reps){
		if (length(fixed.on)!=0){
			startStates[[i]][fixed.on]<-1
		}
		if (length(fixed.off)!=0){
			startStates[[i]][fixed.off]<-0
		}
	}
	attr<-getAttractors(hybrid,method="chosen",startStates=startStates)
	return(list(attr=attr))
}
