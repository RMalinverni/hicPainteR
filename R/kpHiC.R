#permit to plot an HicPainteRObj on karyoploteR

kpHiC<-function(kp,hicPobj){

  if (hicPobj$param$mapType == 'cMap'){

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

  if (hicPobj$param$mapType == 'TAD'){

    chr=seqlevels(hicPobj$cMap[1])
    kpSegments(kp,chr=chr,x0=start(hicPobj$cMap),x1=hicPobj$cMap$edge,y0=0,y1=a$Ys,r0=at$r0,
               r1=at$r1,col=hicPobj$cMap$colors)
    kpSegments(kp,chr=chr,x0=hicPobj$cMap$edge,x1=end(hicPobj$cMap),y0=a$Ys,y1=0,r0=at$r0,
               r1=at$r1,col=hicPobj$cMap$colors)

    }


}
