#create an "HicPainteRObj" object that will be recognize from karyoploteR

createHicPainteRObj<-function(Name,contact,mapType,zoom=NA,high_fact=1e6,pCol="red",alpha=0.05,enhance=NA,colBin=500,
                              use_ramp=FALSE,log=FALSE,ramp=NA,pchP=".",cexP=2,lwdP=2,ltyP=1,scores=NA,use.scores=FALSE,...){


  if (mapType %in% c('cMap','TAD','Loop') == FALSE ){
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

  if ((!is.color(pCol)) & use_ramp==FALSE){
    stop ('"pCol" need to be a valid color or select a valid "ramp" value')
  }

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
