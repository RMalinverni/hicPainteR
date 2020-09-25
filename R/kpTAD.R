
kpTAD<-function(kp,hicPobj,r0=0,r1=1,colorTAD="black",lim_dist = 1e6){

  if  (!is.na(lim_dist)){
    hicPobj$cMap<- hicPobj$cMap %>% filter(end-start<lim_dist)
  }

  chr=seqlevels(hicPobj$cMap[1])
  kpSegments(kp,chr=chr,x0=start(hicPobj$cMap),x1=hicPobj$cMap$edge,
             y0=0,y1=hicPobj$cMap$Ys,r0=r0,
             r1=r1,col=colorTAD,lwd=hicPobj$param$lwdP,lty=hicPobj$param$ltyP)
  kpSegments(kp,chr=chr,x0=hicPobj$cMap$edge,x1=end(hicPobj$cMap),
             y0=hicPobj$cMap$Ys,y1=0,r0=r0,
             r1=r1,col=colorTAD,lwd=hicPobj$param$lwdP,lty=hicPobj$param$ltyP)
  }

