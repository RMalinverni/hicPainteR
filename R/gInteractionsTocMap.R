#' Genomics Interaction to Genomic Map
#'
#' transform a GenomicsInteraction object to a Granges with mCol "edge" and "Ys" that will be use for create a hicPainter Object,
#' if zoom is provided it will extract only this portion of genome
#'
#' @usage gInteractionsTocMap( GI, zoom = NA )
#'
#' @param GI GenomicInteractions object
#' @param zoom character - genomic coordinate in wich the GRanges will be contruct, if not privide it will use all the contact of GI.
#'             The format need to be "chromosome:start-end" ( default = NA )
#'
#'
#'
#' @seealso ...
#'
#' @examples  ...
#'
#' @references ...
#' @importFrom ...
#' @export gInteractionsTocMap


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
  cmap %>% mutate(Ys=Ynorm/5e5) -> cmap
  cmap %>% mutate(Ynorm=(Ynorm-min(Ynorm))/(max(Ynorm)-min(Ynorm))) -> cmap

  return(cmap)
}
