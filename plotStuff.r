plotStuff <- function(stuff,adjust,gaps,rep,fromGen=1,toGen=NULL,ymin=0,main="Fitness of parent populations hybrids"){
if (is.null(toGen)){
  toGen <- adjust+gaps*rep
}
plot(0,0,col="white",ylim=c(ymin,max(stuff[[2]][fromGen:toGen,],stuff[[3]])),xlim=c(fromGen,toGen),ylab="Population average fitness",xlab="Generation",main=main)#,cex.axis=2)
#min(stuff[[2]][fromGen:toGen,],stuff[[3]])
lines(fromGen:toGen,stuff[[2]][fromGen:toGen,1],col=1,type="l",lwd=1)

lines(fromGen:toGen,stuff[[2]][fromGen:toGen,2],col=2,type="l",lwd=1)
if (fromGen<adjust){
lines(c(adjust,adjust),c(stuff[[2]][adjust,1],min(stuff[[2]][fromGen:toGen,],stuff[[3]])),lty=2,lwd=1)
}
if (toGen>adjust){
points(adjust+(1:rep)*gaps,stuff[[3]],col="darkgreen",lwd=2)

}
}

plotRestuff <- function(stuff,gaps,rep){
#subtitle<-paste(nodes,"nodes,",gaps,"gaps,",rep,"rep, max.popn=",max.popn,"prob=(",paste(prob,collapse=","),")")
pop1pheno <- vector("numeric",rep/2)
pop2pheno <- vector("numeric",rep/2)
hybpheno <- vector("numeric",rep/2)

for (i in 1:(rep/2)){
  pop1pheno[i] <- mean(stuff[[2]][((i-1)*gaps+1):((i+1)*gaps),1])
  pop2pheno[i] <- mean(stuff[[2]][((i-1)*gaps+1):((i+1)*gaps),2])
  hybpheno[i] <- mean(stuff[[3]][(i*2)-1:i*2])
}
  
plot(0,0,col="white",ylim=c(min(stuff[[2]],stuff[[3]]),max(stuff[[2]],stuff[[3]])),xlim=c(1,(rep*gaps)))


lines(2*gaps*(1:(rep/2)),pop1pheno,col=1,type="l")

lines(2*gaps*(1:(rep/2)),pop2pheno,col=2,type="l")

points(2*gaps*(1:(rep/2)),hybpheno,col="darkgreen",type="l")
}
