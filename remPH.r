#adapted from the kill edge function, this removes the placeholders on and off
remPH<-function(oldInputFunction,oldInputIndex,oldInputNumber,type){
	
	runLength<-length(oldInputFunction)/2^(oldInputIndex)
	runNumber<-2^(oldInputIndex-1)
	digitOrder<-0

	for (i in 1:runNumber){
		digitOrder[((i-1)*runLength+1):(i*runLength)]<-c(((i-1)*runLength*2+1):((i-1)*runLength*2+runLength))
	}		
		
	newInputFunction<-0	
		
	for (j in 1:length(digitOrder)){
		newInputFunction[j]<-oldInputFunction[digitOrder[j]+runLength*type]
	}
	
	
	return(newInputFunction)
}
