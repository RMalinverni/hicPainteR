# change the parameters of an HicPainterObj without recalculate it

repaintHicPainteRobj<-function(hicPobj,
                               Name=NA,Type=NA,
                               zoom=NA,high_fact=NA,pCol=NA,alpha=NA,enhance=NA,colBin=NA,
                               use_ramp=NA,log=NA,ramp=NA,pchP=NA,cexP=NA,...){

  if (class(hicPobj)!="HicPainteRObj") {
    stop ("the object to repaint need to be a 'HicPainteRObj' class object")
  }

  if ( is.na ( Name )){ Name <- hicPobj$param$Name }
  if (is.na(Type)){ Type<-hicPobj$param$Type }
  if (is.na(zoom)){ zoom<-hicPobj$param$zoom }
  if (is.na(high_fact)){ high_fact<-hicPobj$param$high_fact }
  if (is.na(pCol)){ pCol<-hicPobj$param$pCol }
  if (is.na(alpha)){ alpha<-hicPobj$param$alpha }
  if (is.na(enhance)){ enhance<-hicPobj$param$enhance }
  if (is.na(colBin)){ colBin<-hicPobj$param$colBin }
  if (is.na(use_ramp)){ use_ramp<-hicPobj$param$use_ramp }
  if (is.na(log)){ log<-hicPobj$param$log }
  if (is.na(ramp)){ ramp<-hicPobj$param$ramp }
  if (is.na(pchP)){ pchP<-hicPobj$param$pchP }
  if (is.na(cexP)){ cexP<-hicPobj$param$cexP }

  YS<-hicPobj$cMap$Ys

  GI<-granges(hicPobj$cMap)
  GI$counts<-hicPobj$cMap$orig.counts
  GI$edge<-hicPobj$cMap$edge

  print("re-painting contact map matrix")
  GI<-cMapPainter(GI,high_fact=high_fact,pCol=pCol,alpha=alpha,enhance=enhance,colBin=colBin,
                  use_ramp=use_ramp,log=log,ramp=ramp)
  GI$Ys<-YS
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
