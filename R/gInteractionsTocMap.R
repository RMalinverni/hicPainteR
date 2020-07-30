#transform GenomicInteraction objects in cMap object

gInteractionsTocMap<-function(GI,zoom=NULL){
  seqlevelsStyle(GI)<-'UCSC'
  Sstart<-start(resize(anchors(GI,type="first"),fix="center",width=1))
  Eend<-start(resize(anchors(GI,type="second"),fix="center",width=1))
  Cchrs<-as.character(seqnames(anchors(GI,type="first")))
  cmap<-toGRanges(data.frame(chr=Cchrs,start=Sstart,end=Eend,counts=GI$counts))
  cmap$edge<-start(resize(cmap,fix="center",width=1))
  
  if (!is.null(zoom)){
    zoom<-gsub(",","",zoom)
    chr<-strsplit(zoom,':')[[1]][1]
    SSEE<-strsplit(zoom,':')[[1]][2]
    SS<-as.numeric(strsplit(SSEE,'-')[[1]][1])
    EE<-as.numeric(strsplit(SSEE,'-')[[1]][2])
    cmap %>% filter(edge >= SS & edge <=EE) -> cmap
  }
  
  return(cmap)
}

