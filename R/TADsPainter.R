
TADsPainter<-function(tadGR,high_fact=1e6,colP="black",zoom=NA){

  tadGR %>%  mutate( edge = start( resize ( tadGR ,fix = "center", width = 1))) -> tadGR

  if (!is.na(zoom)){
    zoom<-gsub(",","",zoom)
    chr<-strsplit(zoom,':')[[1]][1]
    SSEE<-strsplit(zoom,':')[[1]][2]
    SS<-as.numeric(strsplit(SSEE,'-')[[1]][1])
    EE<-as.numeric(strsplit(SSEE,'-')[[1]][2])
    cmap %>% filter(edge >= SS & edge <=EE) -> cmap
  }

  tadGR %>% mutate(Ynorm=abs(width(tadGR)/2)) -> tadGR
  tadGR %>% mutate(Ynorm=(Ynorm-min(Ynorm))/(max(Ynorm)-min(Ynorm))) -> tadGR
  tadGR %>% mutate(color=colP)

  return(tadGR)
}
