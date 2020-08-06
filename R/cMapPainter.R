#' Genomic Map Painter
#'
#' add color and graphic parameters to a GenomicRangers object created with gInteractions function
#'
#'
#' @usage cMapPainter(cmap,high_fact=1e6,pCol="red",alpha=0.05,enhance=NA,colBin=500,ramp=NA,use_ramp=FALSE,log=FALSE,zoom=NULL,...)
#'
#' @param cmap GenomicRanges object created using gInteractions function
#' @param zoom character - genomic coordinate in wich the GRanges will be contruct, if not privide it will use all the contact of GI.
#'             The format need to be "chromosome:start-end" ( default = NA )
#' @param pCol charachter - color used if  use_ramp == FALSE ( default "red" )
#' @param alpha numeric (0-1) - is the value used for the transparency, mutuated from \code{\link{ggplot2} (default 0.05)
#' @param enhance numeric - transform factor for counts on cMap ( default = NA )
#' @param colBin numeric - used to calculate the number of the palette when use_ramp == TRUE ( default  500 )
#' @param use_ramp logic - select to switch from flat color to colorRampPalette ( default FALSE )
#' @param ramp colorRampPalette object - palette used when use_ramp == TRUE, when is set to NA ramp= ( default NA )
#' @param log logic - transform the counts of the contact object in logaritmic ( default NA )
#'
#'
#'
#' @seealso ...
#'
#' @examples  ...
#'
#' @references ...
#' @importFrom ...
#' @export cMapPainter

cMapPainter<-function(cmap,high_fact=1e6,pCol="red",alpha=0.05,enhance=NA,colBin=500,
                      ramp=NA,use_ramp=FALSE,log=FALSE,zoom=NULL,...){

  if (is.na(ramp)){
    ramp<-colorRampPalette(colors=c("blue","lightblue","yellow","orange","red","brown"),interpolate="spline",bias=1)(colBin)
  }


  #widthF<-(end(cmap)-start(cmap))/2
  #cmap$Ys<-abs((widthF)/high_fact)

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
