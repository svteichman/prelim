## process HIV data from http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006/DATA/PI_DATA.txt to send to matlab

hivdata=read.table('PI/hivdata.txt',fill=TRUE,sep='\t',header=TRUE)

y_cols=4:10;x_cols=11:109 # locations in the data matrix

is.letter <- function(x) grepl("[[:alpha:]]", x)

rmv=which(rowSums(hivdata[,x_cols]=='Krk')==1)
# remove the one sample with nonstandard mutation names
hivdata=hivdata[-rmv,]

#mutnums=read.table('~/Dropbox/matlab/HIVdata/mut_to_num.txt')
mutnums=read.table('PI/treatment_mut_table.txt')

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


write.table(X,file='~/Desktop/X.txt',row.names=FALSE,col.names=FALSE)
write.table(Y,file='~/Desktop/Y.txt',row.names=FALSE,col.names=FALSE)

write.table(X,file='~/Dropbox/matlab/HIVdata/X.txt',row.names=FALSE,col.names=FALSE)
write.table(Y,file='~/Dropbox/matlab/HIVdata/Y.txt',row.names=FALSE,col.names=FALSE)
write.table(Xnames,file='~/Dropbox/matlab/HIVdata/Xnames.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(Xnums,file='~/Dropbox/matlab/HIVdata/Xnums.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(colnames(hivdata)[y_cols],file='~/Dropbox/matlab/HIVdata/Ynames.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)



