createZooms<-function(zoom){

  if (class(zoom)=="GRanges"){
    GR<-zoom
    chr<-as.character(seqnames(zoom)[1])
    start<-start(zoom)[1]
    end<-end(zoom)[1]
    zoomString<-paste0(chr,":",start,"-",end)
  }

  if (is.character(zoom)){
    zoomString<-gsub(",","",zoom)
    chr<-strsplit(zoomString,":")[[1]][1]
    range<-strsplit(zoomString,":")[[1]][2]
    start<-strsplit(range,"-")[[1]][1]
    end<-strsplit(range,"-")[[1]][2]
    GR<-regioneR::toGRanges(data.frame(chr=chr,start=start,end=end))
  }

  zoomObj<-list(
    zoomString=zoomString,
    chr=chr,
    start=start,
    end=end,
    GR=GR)

  return(zoomObj)
}
