setEnv<-function(network,envFunc){

oldGenes<-network$genes

for (i in 1:length(envFunc)){
	network$genes[i]<-paste("env ",i," ",sep="")
	network$interactions[[i]]$input<-0
	network$interactions[[i]]$func<-envFunc[i]
	network$interactions[[i]]$expression<-"environment"
	network$fixed[i]<-envFunc[i]
	}
return(network)
}