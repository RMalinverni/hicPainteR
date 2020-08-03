#create an "HicPainteRObj" object that will be recognize using a from karyoploteR

createHicPainteRObj<-function(Name,Type=c('HiC','TAD','Loop'),contact,
                             zoom=NA,high_fact=1e6,pCol="red",alpha=0.05,enhance=NA,colBin=500,
                             use_ramp=FALSE,log=FALSE,ramp=NA,pchP=".",cexP=2,...){
  if (class(contact)!= "GenomicInteractions") {
    stop('contact value neet to be Genomicinteraction object')
  }

  if (Type %in% c('HiC','TAD','Loop') ==FALSE){
    stop ('"Type" need to be one of the follows: Hic, TAD or Loop')
  }

  #contact add an error
  #ramp
  #pchP
  #cexP
  #enhance
  if (!is.numeric(high_fact)){
    stop ('"high_fact" need to be numeric')
  }

  if ((!is.color(pCol) & use_ramp=FALSE)){
    stop ('"pCol" need to be a valid color or select a valid "ramp" value')
  }

  if (!is.logical(use_ramp)){
    stop ('"use_ramp" need to be a locical value')
  }

  if (!is.logical(use_ramp)){
    stop ('"use_ramp" need to be a logical value')
  }

  if (!is.logical(log)){
    stop ('"log" need to be a logical value')
  }

  if (Type=='HiC'){
    print("reading the contact matrix")
    GI<-gInteractionTocMap(contact,zoom=zoom)
    print("painting contact map matrix")
    GI<-cMapPainter(GI,high_fact=high_fact,pCol=pCol,alpha=alpha,enanche=enanche,colBin=colBin,
                    use_ramp=use_ramp,log=log,ramp=ramp)
  }
  HPobj<-list(cMap=GI,param=list(Name=Name,
                                 Type=Type,
                                 zoom=zoom,
                                 high_fact=high_fact,
                                 pCol=pCol,
                                 alpha=alpha,
                                 enhance=enhance,
                                 colBin=colBin,
                                 use_ramp=use_ramp,
                                 log=log,ramp=ramp,
                                 pchP=pchP,cexP=cexP))
  class(HPobj)<-"HicPainteRObj"

  return(HPobj)
}
