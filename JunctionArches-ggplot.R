library(ggplot2)
library(GenomicRanges)

##create a GRangesList
gr<-GRanges("chr1", IRanges(c(1,20,50,60), width=c(5,10,6,7)))
combs<-combn(1:4,2)
grl<-do.call("GRangesList", apply(combs,2,function(x){
	res<-gr[x]
	values(res)<-data.frame(pvalue=runif(2))
	res
}))
values(grl)<-data.frame(counts=sample(1:100,size=6), score=rnorm(6))

#values(grl)$counts - support
#start(grl) - start values
#end(grl)  - end values


