> grid.lines(x,y, gp=gpar(col=rgb(0.9,0.9,0.9)))


library(ggplot2)
library(grid)

##draws the arches, half circles representing junction reads
halfcirc<-function(x,y,xr,yr){
	xx<-c()
	yy<-c()
	n<-500
	for(i in 1:n){
		ang<-i*pi/(2*n)
		xx[i]<-cos(ang)
		yy[i]<-sin(ang)
	}
	xx<-c(xx,rev(-xx))
	yy<-c(yy,rev(yy))
	xx<-xr*xx+x
	yy<-yr*yy+y
	shade<-c((8/0.75)*yr)
	qplot(xx, yy, geom="line")
	#grid.lines(xx,yy)
}

halfCircleArch<-function(x1, x2, y,h){
	xx<-c()
	yy<-c()
	n<-1000
	for(i in 1:n){
		ang<-i*pi/n
		xx[i]<-cos(ang)
		yy[i]<-sin(ang)
	}
	xx<-(((x2-x1)/2)-(0.4/(length(x2)+1)))*xx+(x2-x1)/2
	yy<-h*yy+y
	grid.lines(xx,yy,gp=gpar(col=rgb(0.9,0.9,0.9)))
}

## draws the exons as evenly spaced boxes, 
## will call the halfcirc option to draw arches
drawExons<-function(g){
	#stores number of exons in n
	#assigns an id number to each exon.
	key<-sort(unique(start(unlist(g))))
	n<-length(key)
	
	#draws the exons
	x<-c(1:n)
	x<-x/(n+1)
	y<-c(rep.int(0.2,n))
	#grid.rect(x=x, y=y, width=(0.8/(n+1)), height=0.025, just="centre", default.units="npc", draw=TRUE)
	
	#set units for arch height
	#draws in the arches
	arch<-values(g)$counts
	aunit<-c(0.78/(max(values(g)$counts)+1))
	
	#draws the arches using halfCirc() and key
	for(i in 1:length(arch)) {
		left=c()
		right=c()
		for(j in 1:length(key)){
			if(min(start(g[[i]]))==key[j]) {left <- c(j);}	
			if(max(start(g[[i]]))==key[j]) {right <- c(j);}
		}	
		halfcirc((x[right]+x[left])/2, y+0.0125, 
		((x[right]-x[left])/2)-(0.4/(n+1)),arch[i]*aunit)
		###

	}
	
	#adds the labels
	#grid.text(label=c(1:n), x=x, y=y, default.units="npc", draw=TRUE)
}

drawExons(grl)
drawExons(append(com, grl3))



## creating a GRangesList
## Thanks TengFei! for sending me the code.
####################
## for now, store whatever information into GRangesList
## suppose we have exons
## [1]---[2]---[3]---[4]
## We have splicing [1]--[2]
## [2]--[3],..., we put whatever we found as GRangesList, each entry
## contains one splicing form

library(GenomicRanges)
gr<-GRanges("chr1", IRanges(c(1,20,50,60), width=c(5,10,6,7)))
combs<-combn(1:4,2)
combs
grl<-do.call("GRangesList", apply(combs,2,function(x){
	res<-gr[x]
	values(res)<-data.frame(pvalue=runif(2))
	res
}))
values(grl)<-data.frame(counts=sample(1:100,size=6), score=rnorm(6))
grl
end(grl)
values(grl)$counts



grl[1]$seqnames

grl[1]
#first granges object
start(head(grl[1], 1))
first granges object, and returns start
unique()














####################
## for now, store whatever information into GRangesList
## suppose we have exons
## [1]---[2]---[3]---[4]
## We have splicing [1]--[2]
## [2]--[3],..., we put whatever we found as GRangesList, each entry
## contains one splicing form

library(GenomicRanges)
gr <- GRanges("chr1", IRanges(c(1, 20, 50, 60),width = c(5, 10, 6, 7)))

combs <- combn(1:4,2)
combs

grl <- do.call("GRangesList", apply(combs, 2, function(x){
  res <- gr[x]
  values(res) <- data.frame(pvalue = runif(2))
  res
}))

values(grl) <- data.frame(counts = sample(1:100, size = 6),
                                   score  = rnorm(6))


## so you could work on this grl, which contains 6 conbination of 4 exons
## use values to get the counts for each case
grl
## get first combination
grl[[1]]  # '[[' return a GRanges
grl[1] #'[' still a list 
values(grl[1])$counts # counts for first  combination
## notice the elementMetadata( or values) is not the same for list and element
values(grl)
values(grl[1])
## this is for GRanges inside the list
values(grl[[1]])
values(grl)


