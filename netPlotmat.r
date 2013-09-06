#this function handles networks generated from an allele vector
netPlotmat<-function(network,alvec=NULL,geneNumber=NULL){
	N<-length(network$genes)
	if (is.null(alvec)){
		alvec<-network$genes
	}
	if (is.null(geneNumber)){
		#giving the number of alleles as the gene number, that isn't quite right, but only matters for the colours
		geneNumber<-length(alvec)
#		geneNumber<-max(as.integer(sub("g","",sapply(splitName,get.first))))
	}

	geneColours<-rainbow(geneNumber,0.8,0.8)
	
	nodeSize<-rep(round(50/sqrt(length(network$genes))),N)
	nodeColor<-rep("white",N)
	nodeShape<-rep("circle",N)
	frameColor<-rep("black",N)
	textSize<-round(4/sqrt(N))
	
	splitName<-strsplit(network$genes,split="_")
	nodeLabel<-vector('character',length(network$genes))
	
	for (i in 1:N){
		if (length(grep("e",network$genes[i]))>0){
			
			copyNumber<-length(which(tolower(alvec)==network$genes[i]))
			nodeColor[i]<-colors()[145]  #environment nodes --> gold
			nodeSize[i]<-round(50/sqrt(length(network$genes))*sqrt(copyNumber))
			nodeShape[i]<-"square"
			network$genes[i]<-toupper(splitName[[i]][1])
		}
		else if (network$genes[i] %in% c("on","off")){
			nodeLabel[i]<-network$genes[i]
		}
		else {
			copyNumber<-length(which(tolower(alvec)==network$genes[i]))
			nodeSize[i]<-round(nodeSize[i]*sqrt(copyNumber))
			nodeColor[i]<-geneColours[as.integer(sub("g","",splitName[[i]][1]))]
			network$genes[i]<-toupper(splitName[[i]][1])
			}
		
		if (length(network$interactions[[i]]$func)==1 & sum(network$interactions[[i]]$func)==1){
			frameColor[i]<-colors()[630]} #nodes that are always on --> pale red
		else if (length(network$interactions[[i]]$func)==1 & sum(network$interactions[[i]]$func)==0){
			frameColor[i]<-colors()[240]} #nodes that are always off --> light gray
		}

	plotNetworkWiring(network, vertex.color=nodeColor, vertex.size=nodeSize, vertex.frame.color=frameColor, vertex.shape=nodeShape, vertex.label.cex=textSize)
}
get.first<-function(vec){
	return(as.character(vec[1]))
}