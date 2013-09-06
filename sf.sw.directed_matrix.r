#see generating.graphs documentation for more details on all of these variables

sf.sw<-function(N,m0,m,p){
	adj.mat<-matrix(0,N,N)
	#choose M for first time round (it must be < m0)
	r<-runif(1)
	if (m0<m){
		m0<-m
	}
	for (i in (m0+1):N){
		#i is index of current vertex
		new.edges<-0
		while (new.edges<m){
			#preferential attachment. non.adj.mat is a matrix where rows of vertices adjacent to i and with degree already > 15 have been removed
			adj.mat.crop<-adj.mat[1:i-1,]
			adj.indices<-which(adj.mat.crop[,i]==1)
			#indices2<-which(adj.mat[,i]==0)
			indices.to.remove<-which(adj.mat.crop[,i]==1 | apply(adj.mat.crop,1,sum)>=15)
			indices.to.keep<-which(adj.mat[,i]==0 & apply(adj.mat,1,sum)<15)
			if (length(indices.to.remove)==0){
				non.adj.mat<-adj.mat.crop
			}
			else{
				non.adj.mat<-adj.mat.crop[-indices.to.remove,]
			}
			if (is.null(dim(non.adj.mat))){
				relevant.degrees<-sum(non.adj.mat)
			}
			else{
				relevant.degrees<-apply(non.adj.mat,1,sum)
			}
			degrees<-apply(adj.mat[1:i-1,],1,sum)
			total.degree<-sum(relevant.degrees)
			rand<-runif(1)
			cum.deg<-relevant.degrees[1]
			j<-1
			while (cum.deg<total.degree){
				if (rand<cum.deg/total.degree){
					break
				}
				j<-j+1
				cum.deg<-cum.deg + relevant.degrees[j]
			}
			#so j should be then index in relevant.degrees of the edge we will connect to
			new.edge.vertex<-indices.to.keep[j]
			#add this new edge to the matrix
			new.edges<-new.edges+1
			adj.mat[i,new.edge.vertex]<-1
			adj.mat[new.edge.vertex,i]<-1
			#triad formation
			new.edge.vertex2<-0
			if (runif(1)<p & new.edges<m){
				#find out how many candidate vertices there are to connect to
				adj.pos.v<-adj.mat[1:i-1,new.edge.vertex]-adj.mat[1:i-1,i]
				ind.pos.v<-which(adj.pos.v==1 & apply(adj.mat[1:i-1,],1,sum)<15)
				num.pos.v<-length(ind.pos.v)
				#pick one at random and connect to i
				if (num.pos.v>0){
					rand1<-runif(1)
					ind.v<-ceiling(rand1*num.pos.v)
					new.edge.vertex2<-ind.pos.v[ind.v]
					new.edges<-new.edges+1
					adj.mat[i,new.edge.vertex2]<-1
					adj.mat[new.edge.vertex2,i]<-1
				}	
			}
			#print(c(i,new.edge.vertex,new.edge.vertex2))
		}
	}
	return(list(network=adj.mat))
}

make.directed<-function(ntwk,d){
	size<-dim(ntwk)
	N<-size[1]
	for (i in 1:N){
		nbhd<-which(ntwk[,i]==1)
		non.nbhd<-which(ntwk[,i]==0)
		for (j in nbhd){
			if (runif(1)<d){
				#remove i from the possible set of new vertices
				i.index<-which(non.nbhd==i)
				non.nbhd<-non.nbhd[-i.index]
				#randomly choose a new vertex
				total.poss<-length(non.nbhd)
				rand<-runif(1)
				ind.v<-ceiling(rand*total.poss)
				new.v<-non.nbhd[ind.v]
				#remove directed edge from i to j
				ntwk[i,j]<-0
				#add directed edge from i to new.v
				ntwk[i,new.v]<-1
			}
		}
	}
	return(list(network=ntwk))
}

create.logical<-function(ntwk,envNumber){
	size<-dim(ntwk)
	N<-size[1]
	genes<-vector('character',N)
	exp<-vector('character',N)
	colnames(ntwk)<-genes
	rownames(ntwk)<-genes
	#inMat<-ntwk
	#inMat[,1:envNumber]<-matrix(0,length(genes),envNumber)
	for (i in 1:envNumber){
		genes[i]<-paste("e",i,"_",i,sep="")
	}
	for (i in (envNumber+1):N){
		genes[i]<-paste("g",i,"_",i,sep="")
	}
	for (i in 1:envNumber){
		inputs<-0
		func<-as.character(sample(0:1,1))#this probably needs changing
		exp[i]<-func
	}
	for (i in (envNumber+1):N){
		inputs<-which(ntwk[,i]==1)
		num.in<-length(inputs)
		if (num.in ==0){
			inputs<-0
			func<-as.character(sample(0:1,1))#this probably needs changing
		}
		else if (num.in == 1){
			inputs<-genes[inputs]
			if (sample(0:1,1)==1){
				func<-paste("!",strsplit(inputs,split="_")[[1]][1])
			}
			else{
				func<-strsplit(inputs,split="_")[[1]][1]
			}
		}
		else{
			inputs<-genes[inputs]#now the inputs have their proper names
			rand.inputs<-sample(inputs, length(inputs))#randomise the order
			rand.nots<-sample(0:1,length(inputs),replace=TRUE)
			rand.andor<-sample(0:1,length(inputs)-1,replace=TRUE)#lets have 1 = & and 0 = | though for the purposes of this algorithm they're not much different
			brackets<-vector('numeric',length(inputs)*2)# we will have 0=no bracket, 1=closed bracket, -1=open bracket. odd positions should have 0/-1, even positions should have 0/1
			brackets[1]<--1
			brackets[length(inputs)*2]<-1
			#if there are more than two genes we'll need more brackets
			if (length(inputs)>2){
				for (j in 1:(length(rand.andor)-1)){
					if (abs(brackets[2])>0){
						stop("brackets error")
					}
					#if we're between an & and an |, and there are no brackets either side of the gene.
					if (rand.andor[j]!=rand.andor[j+1] & brackets[2*j+1]==0 & brackets[2*j+2]==0){
						openclose<-sample(c(-1,1),1)
						#note, the bracket is going around the (j+1)th gene
						if (openclose==-1){
							brackets[2*j+1]<-brackets[2*j+1]-1
							next.close<-match(1,brackets[(2*j+4):length(brackets)]>0)
							#total.options<-(next.close-(2*j+4))/2
							place<-sample(seq(1,next.close,2),1)
							brackets[2*j+2+place+1]<-brackets[2*j+2+place+1]+1
							#cat("have placed brackets at",place,2*j+1,"\n")
						}
						else{
							brackets[2*j+2]<-brackets[2*j+2]+1
							num.openclose<-sum(brackets[1:(2*j+1)])#this will help to make sure that all our brackets are in the pairs we intend them to be in
							#find a vector of suitable places to put brackets
							poss.options<-vector('numeric')
							index<-1
							for (k in seq((2*j+1),1,-2)){
								if (sum(brackets[1:k])>num.openclose){
									break
								}
								if (k==1){
									if (0>=num.openclose){
										poss.options[index]<-k
										index<-index+1
									}
								}
								else if (k<2*j+1){
									if (sum(brackets[1:(k-1)])>=num.openclose & (rand.andor[(k-1)/2]!=rand.andor[(k+1)/2])){
										poss.options[index]<-k
										index<-index+1
									}
								}
							}
							if (length(poss.options)==1){
								brackets[poss.options]<-brackets[poss.options]-1
							}
							else{
								place<-sample(poss.options,1)
								brackets[place]<-brackets[place]-1
							}
							#cat("have placed brackets at",place,2*j+2,"\n")
						}
					}
				}
			}
			#now write the expression as a character vector
			func<-vector('character')
			index<-1
			for (j in 1:length(inputs)){
				if (brackets[2*j]<0 | brackets[2*j-1]>0){
					stop("brackets error")
				}
				if (j==1){
					func[1:abs(brackets[1])]<-rep("(",abs(brackets[1]))
					index<-index+abs(brackets[1])
					if (rand.nots[1]==1){
						func[index]<-"!"
						index<-index+1
					}
					func[index]<-strsplit(rand.inputs[1],split="_")[[1]][1]
					index<-index+1
					if (abs(brackets[2])>0){
						cat("something has gone wrong, I shouldn't have a bracket here", "\n")
						func[index:(index-1+abs(brackets[2]))]<-rep(")",abs(brackets[2]))
						index<-index+abs(brackets[2])
					}
				}
				else{
					if (rand.andor[j-1]==0){
						func[index]<-"|"
					}
					else{
						func[index]<-"&"
					}
					index<-index+1
					if (brackets[2*j-1]<0){
						func[index:(index-1+abs(brackets[2*j-1]))]<-rep("(",abs(brackets[2*j-1]))
						index<-index+abs(brackets[2*j-1])
					}
					if (rand.nots[j]==1){
						func[index]<-"!"
						index<-index+1
					}
					func[index]<-strsplit(rand.inputs[j],split="_")[[1]][1]
					index<-index+1
					if (brackets[2*j]>0){
						func[index:(index-1+abs(brackets[2*j]))]<-rep(")",abs(brackets[2*j]))
						index<-index+abs(brackets[2*j])
					}
				}
			}
			func<-paste(func,collapse=" ")
		}	
		exp[i]<-func
	}
	alspace<-data.frame(name=genes,exp=exp,group=c(1:N),stringsAsFactors=FALSE)
	return(alspace)

}

make.boolean<-function(ntwk){
	#NB: The adjacency matrix should have inputs down columns, as detailed in documentation
	size<-dim(ntwk)
	N<-size[1]
	bool.ntwk<-list(interactions=vector('list',N),genes=as.character(1:N),fixed=rep(-1,N))
	for (i in 1:N){
		inputs<-which(ntwk[,i]==1)
		num.in<-length(inputs)
		if (num.in ==0){
			inputs<-0
		}
		func<-sample(0:1,2^(num.in),replace=T)
		bool.ntwk$interactions[[i]]<-list(input=inputs,func=func,expression="")		
	}
	names(bool.ntwk$interactions)<-bool.ntwk$genes
	class(bool.ntwk)<-"BooleanNetwork"
	return(bool.ntwk)
}


bool.to.mat<-function(bool.net){
	N<-length(bool.net$genes)
	adj.mat<-matrix(0,N,N,dimnames=list(bool.net$genes,bool.net$genes))
	func.list<-vector('list',N)
	for (i in 1:N){
		if (bool.net$fixed[i]==-1){
			ins<-bool.net$interactions[[i]]$input
			adj.mat[ins,i]<-1
			func.list[[i]]<-bool.net$interactions[[i]]$func
	    }
	}
	return(list(mat=adj.mat,func=func.list))
}

local.clustering.weighted<-function(ntwk,type="either"){
	#note that the clustering coeff is weighted by nbhd size as we divide by
	#k-1 instead of by k(k-1)
	if (class(ntwk)=="BooleanNetwork"){
		ntwk<-bool.to.mat(ntwk)$mat
	}
	N<-dim(ntwk)[1]
	coeffs<-vector('numeric',N)
	if (type=="out"){
		for (i in 1:N){
			#look at the row of the vertex of interest to get nbhd
			nbhd<-ntwk[i,]
			k<-sum(nbhd)
			if (k>0){
				vertices<-which(nbhd==1)
				nbhd.mat<-ntwk[vertices,vertices]
				coeffs[i]<-sum(nbhd.mat)/(k-1)
			}	
		}
	}
	else if (type=="in"){
		for (i in 1:N){
			#look at the column of the vertex of interest to get nbhd
			nbhd<-ntwk[,i]
			k<-sum(nbhd)
			if (k>1){
				vertices<-which(nbhd==1)
				nbhd.mat<-ntwk[vertices,vertices]
				coeffs[i]<-sum(nbhd.mat)/(k-1)
			}
		}
	}
	else if (type=="both"){
		for (i in 1:N){
			#look at the intersection of the row and the column to get nbhd
			nbhd<-(ntwk[i,] & t(ntwk[,i]))
			k<-sum(nbhd)
			if (k>1){
				vertices<-which(nbhd==1)
				nbhd.mat<-ntwk[vertices,vertices]
				coeffs[i]<-sum(nbhd.mat)/(k-1)
			}
		}
	}
	else if (type=="either"){
		for (i in 1:N){
			#look at the intersection of the row and the column to get nbhd
			nbhd<-(ntwk[i,] | t(ntwk[,i]))
			k<-sum(nbhd)
			if (k>1){
				vertices<-which(nbhd==1)
				nbhd.mat<-ntwk[vertices,vertices]
				coeffs[i]<-sum(nbhd.mat)/(k-1)
			}
		}
	}
	return(list(coeffs,ntwk))
}

##edit so that doesn't pick env nodes
select.modules<-function(ntwk,num,type="either"){
	env.nodes<-grep("e",ntwk$genes)
	ret<-local.clustering.weighted(ntwk,type)
	coeffs<-ret[[1]]
	mat<-ret[[2]]
	coeffs.sorted<-sort(coeffs,decreasing=TRUE,index.return=TRUE)
	indices<-coeffs.sorted[[2]]
	nbhds<-vector('list',num)
	num.selected<-0
	index<-1
	already.in.module<-0
	while (num.selected<num & index<= length(coeffs)){
		# we want the selected node to not be an environment node and to not be a neighbour of something already selected
		if (length(which(env.nodes==indices[index]))==0){
			node.already.in.module<-FALSE
			if (length(which(already.in.module==indices[index]))!=0){
				node.already.in.module<-TRUE
			}
			if (node.already.in.module==FALSE){
				num.selected<-num.selected+1
				#this is complicated because we want to select the bi-directional nbhd (that's not a word, I know)
				#but we also don't want to select things that are already.in.module
				new.nbhd<-c(indices[index],which(((mat[indices[index],] | mat[,indices[index]])==1)))
				to.remove<-vector('numeric')
				rm.index<-1
				for (i in 2:length(new.nbhd)){
					if (length(which(already.in.module==new.nbhd[i]))!=0 | length(which(env.nodes==new.nbhd[i]))!=0){
						to.remove[rm.index]<-i
						rm.index<-rm.index+1
					}
				}
				if (length(to.remove)>0){
					nbhds[[num.selected]]<-new.nbhd[-to.remove]
				}
				else{
					nbhds[[num.selected]]<-new.nbhd
				}
				already.in.module<-c(already.in.module,nbhds[[num.selected]])
			}
		}
		index<-index+1
	}
	return(nbhds)
}


create.bool.ntwk<-function(N,m0,m,p,d){
	a<-sf.sw(N,m0,m,p)
	b<-make.directed(a$network,d)
	c<-make.boolean(t(b$network))
	for (i in 1:N){
		c$interactions[[i]]$expression=paste("<f(",paste(c$interactions[[i]]$input,collapse=","),"){",paste(c$interactions[[i]]$func,collapse=""),"}>",sep="")}
	d<-simplifyNetwork(c)#this is not so desirable and could potentially be removed
	return(d)
}

create.allele.ntwk<-function(N,m0,m,p,d,envNumber){
	a<-sf.sw(N,m0,m,p)
	b<-make.directed(a$network,d)
	funcs<-create.logical(t(b$network),envNumber)
	inMat<-t(b$network)
	inMat[,1:envNumber]<-matrix(0,N,envNumber)
	alspace<-funcs
	alvec<-as.character(funcs$name)
	colnames(inMat)<-alvec
	rownames(inMat)<-alvec
	return(list(alspace=alspace,alvec=alvec,inMat=inMat))#probably don't need all these returns but keep them for now
}

get.bracket.groups<-function(func){
	#func is a logical expression string i.e. "( G2 & ! G3 )"
	chars<-strsplit(func,split=" ")[[1]]
	brackets<-vector('numeric')#this will help is put brackets into pairs
	brackets[1]<-0
	brackets.loc<-vector('numeric')
	brackets.loc[1]<-0 # just for convenience to make things match up
	b.index<-2
	g.index<-1
	for (i in 1:length(chars)){
		#record where the brackets are
		if (chars[i]=="("){
			brackets[b.index]<-brackets[b.index-1]-1
			brackets.loc[b.index]<-i
			b.index<-b.index+1
		}
		else if (chars[i]==")"){
			brackets[b.index]<-brackets[b.index-1]+1
			brackets.loc[b.index]<-i
			b.index<-b.index+1
		}
	}
	pairs<-list() #this will be a list of pairs of brackets
	func.groups<-list()
	pairs.index<-1
	while (any(brackets<0)==TRUE){
		mini<-min(brackets)
		locs<-which(brackets==mini)
		for (j in 1:length(locs)){
			if (brackets[locs[j]-1]>brackets[locs[j]]){
				#need to locate the matching bracket
				locs2<-which(brackets[(locs[j]+1):length(brackets)]>mini)[1]
				pairs[[pairs.index]]<-c(brackets.loc[locs[j]],brackets.loc[locs[j]+locs2])
				func.groups[[pairs.index]]<-paste(chars[pairs[[pairs.index]][1]:pairs[[pairs.index]][2]],collapse=" ")
				pairs.index<-pairs.index+1
				brackets[locs[j]:(locs[j]+locs2-1)]<-brackets[locs[j]:(locs[j]+locs2-1)]+1
			}
		}
	}
	return(func.groups)	
}

