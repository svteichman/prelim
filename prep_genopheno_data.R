setwd('~/git/prelim/data/')


class=c('PI','NRTI','NNRTI')
datafileabbrv=c('PR','RT','RT')

drugcols=list()
drugcols[[1]]=c(7,9,11,13,15,17,19,21,23)
drugcols[[2]]=c(7,9,11,13,15,17,21,23)
drugcols[[3]]=c(19,25,27,29,31)

npos=c(99,560,560);startpos=c(25,33,33)


mutnames=c(LETTERS,'d','i') # d & i are deletion & insertion

for(i in 1:3){
	data=read.table(paste(datafileabbrv[i],'_data.txt',sep=''),fill=TRUE,header=TRUE,sep='\t')
	n=dim(data)[1]
	Y=log(data[,drugcols[[i]]])
	for(j in 1:n){
		Y[j,which(is.na(Y[j,]))]=999
	}
	write.table(Y,paste(class[i],'_Y.txt',sep=''),row.names=FALSE,col.names=FALSE)
	Ynames=unlist(strsplit(names(data)[drugcols[[i]]],'\\.Fold'))
	write.table(Ynames,paste(class[i],'_Ynames.txt',sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE)
	
	X=Xnums=NULL
	for(l in startpos[i]-1+(1:npos[i])){
		Xl=matrix(0,n,length(mutnames))
		col=data[,l]
		for(j in 1:n){
			muts=strsplit(toString(col[j]),'')[[1]]
			for(k in 1:length(muts)){
				if(any(mutnames==muts[k])){
					Xl[j,which(mutnames==muts[k])]=1
				}
			}
		}
		keep_mutnames=which(colSums(Xl)>0)
		if(length(keep_mutnames)>0){
			posnum=as.numeric(strsplit(names(data)[l],'P')[[1]][2])
			Xnums=rbind(Xnums,cbind(posnum,keep_mutnames))
			X=cbind(X,Xl[,keep_mutnames])
		}
	}
	write.table(X,paste(class[i],'_X.txt',sep=''),row.names=FALSE,col.names=FALSE)
	write.table(Xnums,paste(class[i],'_Xnums.txt',sep=''),row.names=FALSE,col.names=FALSE)
}