#add color and normalization to cMap object

cMapPainter<-function(cmap,high_fact=1e6,pCol="red",alpha=0.05,enhance=NA,colBin=500,
                      ramp=NA,use_ramp=FALSE,log=FALSE,zoom=NULL,...){

  if (is.na(ramp)){
    ramp<-colorRampPalette(colors=c("blue","lightblue","yellow","orange","red","brown"),interpolate="spline",bias=1)(colBin)
  }

  cmap %>% mutate(Ynorm=abs(width(cmap)/2)) -> cmap
  cmap %>% mutate(Ynorm=(Ynorm-min(Ynorm))/(max(Ynorm)-min(Ynorm))) -> cmap

  widthF<-(end(cmap)-start(cmap))/2
  cmap$Ys<-abs((widthF)/high_fact)

  cmap %>% mutate(orig.counts=counts) -> cmap

  if (log ==TRUE) {

    cmap %>%  mutate(counts=log2(counts))-> cmap
  }

  if ( is.numeric(enhance )){
    cmap %>%  mutate(counts=counts^enhance)-> cmap
  }

  if ( use_ramp==FALSE ){
    cmap %>% mutate(counts=round(((counts-min(counts))/(max(counts)-min(counts)))*colBin,digits = 0)+1) %>%
      mutate(colors=alpha(pCol,alpha=counts*alpha)) -> cmap
  }

  if ( use_ramp==TRUE ){
    cmap %>% mutate(counts=round(((counts-min(counts))/(max(counts)-min(counts)))*colBin,digits = 0)+1) %>%
      mutate(colors=alpha(ramp[counts],alpha=counts*alpha)) -> cmap
  }

  cmap<-sort(cmap)
  return(cmap)
}
