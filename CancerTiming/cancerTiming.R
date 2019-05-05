.libPaths("Software/cancerTiming/Rlib")
library(cancerTiming)

# histories
h2 <- matrix(c(0,2,1,0),ncol=2,nrow=2,byrow=TRUE)
h3 <- matrix(c(1,3,1,0),ncol=2,nrow=2,byrow=TRUE)
h4 <- matrix(c(1,2,4,0,1,0,1,0,0),ncol=3,nrow=3,byrow=TRUE)
h5 <- matrix(c(1,2,3,5,0,0,1,0,0,1,0,0,1,0,0,0),ncol=4,nrow=4,byrow=TRUE)

# read in file
mutData <- read.table("sample1.mutData.xls", sep="\t", header=T)
cnvData <- read.table("sample1.mut_cnv.xls", sep="\t", header=T)

chrome <- as.vector(cnvData$chrome)
start <- cnvData$start
end <- cnvData$end
type <- as.vector(cnvData$type)
cn <- cnvData$tcn
piFinal <- c()
qFinal <- c()
qqFinal <- c()
mutNum <- c()

for(i in 1:length(chrome)){
	onlyMuts <- subset(mutData,chromosome==chrome[i] & position <= end[i] & position >= start[i])
	onlyMuts$t_depth <- onlyMuts$t_ref_count+onlyMuts$t_alt_count
	x <- eventTiming(x=onlyMuts$t_alt_count, m=onlyMuts$t_depth, history=get(paste("h", cn[i], sep="")), method="Bayes",
	 totalCopy=cn[i], type=type[i], normCont=0.585769)
	
	piAll <- c()
	qAll <- c()
	qqAll <- c()
	ci <- c()
	
	pi <- as.matrix(x$pi)
	piName <- rownames(pi)
	ci <- as.matrix(x$piCI)
	ciName <- rownames(ci)
	q <- as.matrix(x$q)
	qq <- as.matrix(x$call$alleleSet)
	qqName <- rownames(qq)
	
	for(j in 1:length(piName)){
		ciLH <- paste(round(ci[j,1],2),round(ci[j,2],2),sep=",")
		ciLH <- paste('[',ciLH,']',sep="")
		piAll <- union(piAll, c(paste(piName[j],round(pi[j],2),ciLH,sep=":")))
		qAll <- union(qAll, c(paste(qqName[j], round(q[j],2), sep=":")))
		qqAll <- union(qqAll, c(paste(qqName[j], round(qq[j],2), sep=":")))
	}

	piFinal[i] <- paste(piAll, collapse="__")
	qFinal[i] <- paste(qAll, collapse="__")
	qqFinal[i] <- paste(qqAll, collapse="__")
	mutNum[i] <- x$summaryTable[2,1]
}

out = data.frame(chrome=chrome,start=start,end=end,copyNumber=cn,mutNum=mutNum,type=type,pi=piFinal,q=qFinal,qNormal=qqFinal)
write.table(out, file="WGC036881.cancerTiming_result.xls",row.names=F,quote=F,sep="\t")

