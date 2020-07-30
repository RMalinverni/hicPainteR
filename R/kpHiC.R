#permit to plot an HicPainteRObj on karyoploteR

kpHiC<-function(kp,hicPobj){
  chr=seqlevels(hicPobj$cMap[1])
  kpPoints(kp,chr=chr,
           x=hicPobj$cMap$edge,
           y=hicPobj$cMap$Ys,
           col=hicPobj$cMap$colors, 
           r0=at$r0, 
           r1=at$r1,
           pch=hicPobj$param$pchP,
           cex=hicPobj$param$cexP)
}
