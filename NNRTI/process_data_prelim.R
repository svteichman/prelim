## process HIV data from http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt to send to matlab

hivdata=read.table('~/Dropbox/matlab/HIV_NNRTI/hivdata.txt',fill=TRUE,sep='\t',header=TRUE)

y_cols=4:6;x_cols=7:246 # locations in the data matrix

is.letter <- function(x) grepl("[[:alpha:]]", x)


mutnums=read.table('~/Dropbox/matlab/HIV_NNRTI/mut_to_num.txt')

X=Xnames=Xnums=NULL;n=dim(hivdata)[1]
for(i in x_cols){
	Xi=NULL
	col=hivdata[,i]
	mutnames=NULL
	for(j in 1:n){
		muts=strsplit(toString(col[j]),'')[[1]]
		for(k in 1:length(muts)){
			if(is.letter(muts[k])){
				if(any(mutnames==muts[k])){
					l=which(mutnames==muts[k])
					Xi[j,l]=1
				}else{
					vec=rep(0,n);vec[j]=1
					Xi=cbind(Xi,vec)
					mutnames=c(mutnames,muts[k])
					mutnum=mutnums[which(mutnums[,1]==muts[k]),2]
					Xnames=rbind(Xnames,c(as.numeric(strsplit(names(hivdata)[i],'P')[[1]][2]),muts[k]))
					Xnums=rbind(Xnums,c(as.numeric(strsplit(names(hivdata)[i],'P')[[1]][2]),mutnum))
				}
			}
		}
	}
	X=cbind(X,Xi)
}

Y=NULL
for(i in y_cols){
	find_na=which(is.na(hivdata[,i]))
	col=log(hivdata[,i]);col[find_na]=999
	Y=cbind(Y,col)
}


write.table(X,file='~/Dropbox/matlab/HIV_NNRTI/X.txt',row.names=FALSE,col.names=FALSE)
write.table(Y,file='~/Dropbox/matlab/HIV_NNRTI/Y.txt',row.names=FALSE,col.names=FALSE)
write.table(Xnames,file='~/Dropbox/matlab/HIV_NNRTI/Xnames.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(Xnums,file='~/Dropbox/matlab/HIV_NNRTI/Xnums.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(colnames(hivdata)[y_cols],file='~/Dropbox/matlab/HIV_NNRTI/Ynames.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)


