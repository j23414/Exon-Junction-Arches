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

#Draws multiple arches.  
#startX= vector of start positions
#endX= vector of end positions
#y= single value for y position
#h= vector of hights of arches
#all vectors must be the same length or this will probably break. 
Arches<-function(startX, endX, y, h){
	#Need to add something here to make sure 
	#startX, endX, and h are the same length vectors.
	
	xx<-c()
	yy<-c()
	
	#number of points used to draw quarter of the curve
	n=500
	for(i in 1:n){
		ang<-i*pi/(2*n)
		xx[i]<-cos(ang)
		yy[i]<-sin(ang)
	}
	
	#sets point for complete curve,makes sides are even.
	xx<-c(1,xx,rev(-xx),-1)
	yy<-c(0,yy,rev(yy), 0)
	apoint<-data.frame()
	for(i in 1:length(startX)){
		temp<-data.frame(xx=xx*(abs(startX[i]-endX[i])/2)+(startX[i]+endX[i])/2,
						 yy=yy*h[i]+y,
						 junc=i)
		apoint<-rbind(apoint,temp)
	}
	ggplot(apoint, aes(xx,yy, group=junc))+geom_line()
}

#Testing code
start <-c(1,2,3)
end   <-c(5,6,7)
height<-c(2,-1.5,1)
y=2
Arches(start,end, y, height)
