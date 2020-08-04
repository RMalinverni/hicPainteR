#create an "HicPainteRObj" object that will be recognize using a from karyoploteR

createHicPainteRObj<-function(Name,contact, zoom=NA,high_fact=1e6,pCol="red",alpha=0.05,enhance=NA,colBin=500,
                             use_ramp=FALSE,log=FALSE,ramp=NA,pchP=".",cexP=2,lwdP=2,ltyP=1,scores=NA,use.scores=FALSE,...){




    if (Type %in% c('HiC','TAD','Loop') == FALSE ){
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

    if (Type=='HiC'){

      if (class(contact)!= "GenomicInteractions") {
        stop('for cMap format, contact value need to be Genomicinteraction object')
      }

      print("reading the contact matrix")
      GI<-gInteractionsTocMap(contact,zoom=zoom)
      print("painting contact map matrix")
      GI<-cMapPainter(GI,high_fact=high_fact,pCol=pCol,alpha=alpha,enhance=enhance,colBin=colBin,
                      use_ramp=use_ramp,log=log,ramp=ramp)


      HPobj<-list(cMap=GI,param=list(Name=Name,
                                   Type='cMap',
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



    if (Type =='TAD'){

      if (class(contact)!= "GenomicRanges") {
        stop('for TAD format, contact value need to be GenomicRanges object')
      }

      GI<-TADsPainter(cmap,high_fact=high_fact,pCol=pCol,zoom=zoom,lwdP=lwdP,ltyP=ltyP)
      HPobj<-list(cMap=GI,param=list(Name=Name,
                                     Type='TAD',
                                     zoom=zoom,
                                     high_fact=high_fact,
                                     pCol=pCol,
                                     lwdP=lwdP,ltyP=ltyP))

    }

  if (Type =='Loop'){

    if (class(contact)!= "GenomicRanges") {
      stop('for Loop format, contact value need to be GenomicRanges object')
    }

    print("reading the contact matrix")
    GI<-gInteractionsTocMap(contact,zoom=zoom)
    print("painting the loops ")
    GI<-LoopsPainter(GI,high_fact=high_fact,pCol=pCol,alpha=alpha,enhance=enhance,colBin=colBin,
                     scores=scores,use.scores=use.scores,log=log)


    HPobj<-list(cMap=GI,param=list(Name=Name,
                                   Type='Loop',
                                   zoom=zoom,
                                   high_fact=high_fact,
                                   pCol=pCol,
                                   alpha=alpha,
                                   enhance=enhance,
                                   colBin=colBin,
                                   scores=scores,
                                   use.scores=use.scores,
                                   log=log))

  }


    class(HPobj)<-"HicPainteRObj"


  return(HPobj)
}
