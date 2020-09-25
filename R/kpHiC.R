#permit to plot an HicPainteRObj on karyoploteR

kpHiC<-function(kp,hicPobj){
kpHiC<-function(kp,hicPobj,r0=0,r1=1,lim_dist=5e6){

  if  (!is.na(lim_dist)){
    hicPobj$cMap<- hicPobj$cMap %>% filter(end-start<lim_dist)
  }

  chr=seqlevels(hicPobj$cMap[1])
  kpPoints(kp,chr=chr,
           x=hicPobj$cMap$edge,
           y=hicPobj$cMap$Ys,
           col=hicPobj$cMap$colors, 
           r0=at$r0, 
           r1=at$r1,
           col=hicPobj$cMap$colors,
           r0=r0,
           r1=r1,
           pch=hicPobj$param$pchP,
           cex=hicPobj$param$cexP)
}
