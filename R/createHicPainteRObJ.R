#' Create an hicPainteR object
#'
#' create an hicPainteR object starting from a GenomicInteractions object
#'
#' @usage createHicPainteRObj((contact,mapType,Name,zoom=NA,high_fact=1e6,pCol="red",alpha=0.05,enhance=NA,colBin=500,
#'        use_ramp=FALSE,ramp=NA,log=FALSE,pchP=".",lwdP=2,ltyP=1,cexP=2,scores=NA,use.scores=FALSE,...))
#'
#' @param contact GenomicInteractions object
#' @param mapType character - possible choises are cMAP, TAD, Loops
#' @param Name character - is the name of the object ( default NA )
#' @param high_fact numeric - value for manage the proportion of graph ( default 1e6 ) (TODO control if making it automatic)
#' @param pCol charachter - color used if  use_ramp == FALSE ( default "red" )
#' @param alpha numeric (0-1) - is the value used for the transparency, mutuated from \code{\link{ggplot2} (default 0.05)
#' @param enhance numeric - transform factor for counts on cMap ( default = NA )
#' @param colBin numeric - used to calculate the number of the palette when use_ramp == TRUE ( default  500 )
#' @param use_ramp logic - select to switch from flat color to colorRampPalette ( default FALSE )
#' @param ramp colorRampPalette object - palette used when use_ramp == TRUE, when is set to NA ramp= ( default NA )
#' @param log logic - transform the counts of the contact object in logaritmic ( default NA )
#' @param pchP character - pch code for points drawing
#' @param lwdP numeric - lwd thick for line drawing
#' @param ltyP numeric - lty code for line drawing
#' @param cexP numeric - cex value for point drawing
#'
#'
#'
#' @seealso ...
#'
#' @examples  ...
#'
#' @references ...
#' @importFrom ...
#' @export createHicPainteRObj

createHicPainteRObj<-function(contact, mapType, Name=NA, zoom=NA, high_fact=1e6, pCol="red", alpha=0.05, enhance=NA,
                              colBin=500, use_ramp=FALSE, ramp=NA, log=FALSE, pchP=".", cexP=2, lwdP=2, ltyP=1,
                              scores=NA, use.scores=FALSE, ...){



  if (mapType %in% c('cMap','TAD','Loop') == FALSE ){
    stop ('"Type" need to be one of the follows: Hic, TAD or Loop')
  }

  if( is.na(Name) ){
    Name<-paste0(mapType,"_hP")
  }
  #contact add an error
  #ramp
  #pchP
  #cexP
  #enhance
  if (!is.numeric(high_fact)){
    stop ('"high_fact" need to be numeric')
  }

  # if ((!is.color(pCol)) & use_ramp==FALSE){  #change this for same results
  #   stop ('"pCol" need to be a valid color or select a valid "ramp" value')
  # }

  if (!is.logical(use_ramp)){
    stop ('"use_ramp" need to be a logical value')
  }

  if (!is.logical(use_ramp)){
    stop ('"use_ramp" need to be a logical value')
  }

  if (!is.logical(log)){
    stop ('"log" need to be a logical value')
  }

  if (mapType=='cMap'){

    if (class(contact)!= "GenomicInteractions") {
      stop('for cMap format, contact value need to be Genomicinteraction object')
    }

    print("reading the contact matrix")
    GI<-gInteractionsTocMap(contact,zoom=zoom)
    print("painting contact map matrix")
    GI<-cMapPainter(GI,high_fact=high_fact,pCol=pCol,alpha=alpha,enhance=enhance,colBin=colBin,
                    use_ramp=use_ramp,log=log,ramp=ramp)


    HPobj<-list(cMap=GI,param=list(Name=Name,
                                   mapType=mapType,
                                   zoom=zoom,
                                   high_fact=high_fact,
                                   pCol=pCol,
                                   alpha=alpha,
                                   enhance=enhance,
                                   colBin=colBin,
                                   use_ramp=use_ramp,
                                   log=log,ramp=ramp,
                                   pchP=pchP,cexP=cexP))

  }



  if (mapType =='TAD'){

    if (class(contact)!= "GenomicInteractions") {
      stop('for TAD format, contact value need to be GenomicInteractions object')
    }
    print("reading the contact matrix")
    GI<-gInteractionsTocMap(contact,zoom=zoom)
    print("painting contact map matrix")
    GI<-TADsPainter(GI,high_fact=high_fact,pCol=pCol)
    HPobj<-list(cMap=GI,param=list(Name=Name,
                                   mapType=mapType,
                                   zoom=zoom,
                                   high_fact=high_fact,
                                   pCol=pCol,
                                   lwdP=lwdP,ltyP=ltyP))

  }

  if (mapType =='Loop'){

    if (class(contact)!= "GenomicInteractions") {
      stop('for Loop format, contact value need to be GenomicRanges object')
    }

    print("reading the contact matrix")
    GI<-gInteractionsTocMap(contact,zoom=zoom)
    print("painting the loops ")
    GI<-loopsPainter(GI,high_fact=high_fact,pCol=pCol,alpha=alpha,enhance=enhance,
                     scores=scores,use.scores=use.scores,log=log,pchP=pchP,cexP=pchP)


    HPobj<-list(cMap=GI,param=list(Name=Name,
                                   mapType=mapType,
                                   zoom=zoom,
                                   high_fact=high_fact,
                                   pCol=pCol,
                                   alpha=alpha,
                                   enhance=enhance,
                                   scores=scores,
                                   use.scores=use.scores,
                                   log=log))

  }


  class(HPobj)<-"HicPainteRObj"


  return(HPobj)
}
