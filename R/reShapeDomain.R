reShapeDomain<-function(DMN,GR,xbin=100){
  yfactor<-(-log10(1/xbin))


  selGR<- GR %>% filter(start(GR) >= start(DMN), end(GR)<= end(DMN))
  lengthTAD<-end(DMN)-start(DMN)
  #selGR<- selGR %>% mutate (y=round((end(selGR)-start(selGR))/lengthTAD,digits = yfactor))
  selGR<- selGR %>% mutate (y=(end(selGR)-start(selGR))/lengthTAD)
  selGR<- selGR %>% mutate(y=round(scales::rescale((selGR$y), to=c(1,xbin))))
  selGR<- selGR %>% mutate (x=((end(selGR)-start(selGR))/2)+start(selGR))

  if(length(selGR)!=0){
    selGR<- selGR %>% mutate (x1=round(scales::rescale(selGR$x, to=c(1,xbin)),digits =0))

    AllBinX<-c(1:xbin)
    #AllBinY<-seq(from = 0,to = 1,by=1/10^yfactor)
    AllBinY<-seq(from = 0,to = xbin,by=1)


    diffY<-setdiff(AllBinY,unique(selGR$y))
    diffX<-setdiff(AllBinX,unique(selGR$x1))

    if (length(diffY)!= 0){
      fakeGR<-selGR[rep(1,length(diffY))]
      fakeGR$y<-diffY
      fakeGR$counts<-0
      fakeGR$x1<-1
      fakeGR$Ynorm<-0
      fakeGR$Ys<-0
      selGR<-c(selGR,fakeGR)
    }

    if (length(diffX)!= 0){
      fakeGR<-selGR[rep(1,length(diffX))]
      fakeGR$x1<-diffX
      fakeGR$counts<-0
      fakeGR$y<-0
      fakeGR$Ynorm<-0
      fakeGR$Ys<-0
      selGR<-c(selGR,fakeGR)
    }
  }else{

    nelementGR<-max(length(seq(from = 0,to = xbin,by=1)),xbin)
    selGR<-toGRanges(data.frame(chr=rep(seqlevels(DMN)[1],nelementGR),
                                start=rep(start(DMN)[1],nelementGR),
                                end=rep(end(DMN)[1],nelementGR),
                                counts=rep(0,nelementGR),
                                Ynorm=rep(0,nelementGR),
                                Ys=rep(0,nelementGR),
                                y=seq(from = 0,to = 1,by=1/10^yfactor),
                                x=rep(0,nelementGR),
                                x1=c(1,rep(1:xbin))))
    selGR$x1[selGR$x1 > xbin]<-xbin
  }


  return(selGR)
}
