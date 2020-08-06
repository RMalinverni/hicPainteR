loopsPainter<-function(cmap,high_fact=1e6,pCol="red",alpha=1,enhance=NA,colBin=500,
                       scores=NA,use.scores=NA,log=FALSE,...){

  cmap %>% mutate(orig.counts=counts) -> cmap

  if (log ==TRUE) {

    cmap %>%  mutate(counts=log2(counts))-> cmap

  }

  if ( is.numeric(enhance )){
    cmap %>%  mutate(counts=counts^enhance)-> cmap
  }

  if ( !is.na (scores) ){
    cmap %>% mutate(score=scores) -> cmap
    if ( use_score==TRUE){

      cmap %>% mutate(counts=counts*score) -> cmap

    }
  }

  cmap %>% mutate(counts=round(((counts-min(counts))/(max(counts)-min(counts)))*colBin,digits = 0)+1) %>%
    mutate(colors=alpha(pCol,alpha=counts*alpha)) -> cmap


  cmap<-sort(cmap)

  return(cmap)
}
