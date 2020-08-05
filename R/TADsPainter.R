

TADsPainter<-function(tadGR,high_fact=1e6,pCol="black"){

  tadGR %>% mutate(colors=pCol)

  return(tadGR)
}
