arrowHead2GI<-function(bedpe){
#TAD
  if(!is.data.frame(bedpe)){
    print("uploading dump file......")
    bedpe<-read.delim(bedpe,header=TRUE)
    bedpe<-bedpe[-1,]
  }
  print(nrow(bedpe))
  print(counts)
  GInt<-toGenomicInterations(bedpe,counts=NA)
  GR1<-anchors(GInt,type="first")
  end(GR1)<-start(GR1)
  GR2<-anchors(GInt,type="second")
  start(GR2)<-end(GR2)
  GInt2<-GenomicInteractions(anchor1 = GR1,anchor2 = GR2,counts = GInt$counts)
  return(GInt2)

}


