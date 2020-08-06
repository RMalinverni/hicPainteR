kpLoops<-function(kp,hicPobj,loopType='Arcs'){

  if (loopType == 'Arcs'){
    cM1<-cM2<-granges(hicPobj$cMap)
    end(cM1)<-start(cM1)+2e4
    start(cM2)<-end(cM2)-2e4
    kpPlotLinks(kp,data=cM1,data2=cM2,r0=at$r0,r1=at$r1,y=0,arch.height = hicPobj$cMap$Ys*5,border=NA)
  }
  if (loopType == 'Points'){
    chr=seqlevels(hicPobj$cMap)
    kpPoints(kp,chr=chr,x = hicPobj$cMap$edge, y = hicPobj$cMap$Ys,pch = hicPobj$param$pchP,cex = hicPobj$param$cexP, col=hicPobj$param$colP,r0=at$r0,r1=at$r1)
   }
}





