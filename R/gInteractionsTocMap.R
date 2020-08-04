#transform GenomicInteraction objects in cMap object

gInteractionsTocMap<-function(GI,zoom=NA){
  seqlevelsStyle(GI)<-'UCSC'
  Sstart<-start(resize(anchors(GI,type="first"),fix="center",width=1))
  Eend<-start(resize(anchors(GI,type="second"),fix="center",width=1))
  Cchrs<-as.character(seqnames(anchors(GI,type="first")))
  cmap<-toGRanges(data.frame(chr=Cchrs,start=Sstart,end=Eend,counts=GI$counts))
  cmap %>%  mutate( edge = start( resize ( cmap ,fix = "center", width = 1))) -> cmap


  if (!is.na(zoom)){
    zoom<-gsub(",","",zoom)
    chr<-strsplit(zoom,':')[[1]][1]
    SSEE<-strsplit(zoom,':')[[1]][2]
    SS<-as.numeric(strsplit(SSEE,'-')[[1]][1])
    EE<-as.numeric(strsplit(SSEE,'-')[[1]][2])
    cmap %>% filter(edge >= SS & edge <=EE) -> cmap
  }

  cmap %>% mutate(Ynorm=abs(width(cmap)/2)) -> cmap
  cmap %>% mutate(Ys=Ynorm/1e6) -> cmap
  cmap %>% mutate(Ynorm=(Ynorm-min(Ynorm))/(max(Ynorm)-min(Ynorm))) -> cmap

  return(cmap)
}

