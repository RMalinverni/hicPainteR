#import a dump file from .hic and transform it in GenomicInteractions Object

hicDump2GI<-function(dump,chr){
  if(!is.data.frame(dump)){
    print("uploading dump file......")
    dump<-read.delim(dump,header=FALSE)
  }
  anchor1<-toGRanges(data.frame(chr=chr,start=dump[,1],end=dump[,1]))
  anchor2<-toGRanges(data.frame(chr=chr,start=dump[,2],end=dump[,2]))
  GInt<-GenomicInteractions(anchor1 = anchor1,anchor2 = anchor2,counts = dump[,3])
  rm(anchor1,anchor2)
  secondCenter<-start(resize(anchors(GInt,type="second"),fix="center",width=1))
  firstCenter<-start(resize(anchors(GInt,type="first"),fix="center",width=1))
  gap<-(secondCenter-firstCenter)/2

  if(sum(is.nan(GInt$counts))!=0){    # this permit to graph also the not normilized dump (hichip for example)
    GInt<-GInt[-which(is.nan(GInt$counts))]
  }
  return(GInt)
}
