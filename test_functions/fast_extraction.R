library(strawr)
library(regioneR)
library(GenomicInteractions)
library(plyranges)
library(tidyverse)
library(karyoploteR)

TAD<-
zoomString<-"chr2:175,849,978-178,182,899"

createZooms<-function(zoomString){
  zoom<-gsub(",","",zoomString)
  chr<-strsplit(zoom,":")[[1]][1]
  range<-strsplit(zoom,":")[[1]][2]
  start<-strsplit(range,"-")[[1]][1]
  end<-strsplit(range,"-")[[1]][2]
  GR<-regioneR::toGRanges(data.frame(chr=chr,start=start,end=end))
  zoomObj<-list(
    zoom=zoom,
    chr=chr,
    start=start,
    end=end,
    GR=GR)
  return(zoomObj)
}

reShapeTAD<-function(TAD,GR_test,xbin=10000,yfactor=2){
  selGR<- GR_test %>% filter(start(GR_test) >= start(TAD), end(GR_test)<= end(TAD))
  lengthTAD<-end(TAD)-start(TAD)
  selGR<- selGR %>% mutate (y=round((end(selGR)-start(selGR))/lengthTAD,digits = yfactor))
  selGR<- selGR %>% mutate (x=((end(selGR)-start(selGR))/2)+start(selGR))
  selGR<- selGR %>% mutate (x1=round(scales::rescale(selGR$x, to=c(1,xbin)),digits = 0))
  return(selGR)
}

comulativeMatrix<-function(GR,GI,bin=1000,yfactor=2,verbose=FALSE){
  mat<-matrix(nrow = (10^yfactor)+1,ncol=bin,data = 0)
  colnames(mat)<-seq(1,bin)
  rownames(mat)<-seq(0.0,1,10^-yfactor)

  if (verbose==TRUE){cn<-0}

  chrs<-as.vector(seqnames(anchors(GI,type="first")))
  GR_3D<-toGRanges(data.frame(chrs,start(anchors(GI,"first")),start(anchors(GI,"second"))))
  GR_3D<- GR_3D %>% mutate(counts=GI$counts)

  for( i in 1:length(GR)){
    if (verbose==TRUE){
      cn<-cn+1
      print(cn)
    }
    GRi<-GR[i]
    selGR<-reShapeTAD(TAD=GRi,GR_test = GR_3D,xbin = bin, yfactor=yfactor)
    for(j in 1:length(selGR)){
      matX<-as.character(selGR$y[j])
      matY<-as.character(selGR$x1[j])
      mat[matX,matY]<-mat[matX,matY] + selGR$counts[j]
    }
  }

  return(mat)
}


DFfromHIC<-function(hicFile,zoom,seq_style="UCSC",bin=10000,norm="KR"){

  if (is.character(zoom)){
    zoom<-createZooms(zoom)
  }
  seqlevelsStyle(zoom$GR)<-seq_style

  chr1loc<-seqlevels(zoom$GR)
  chr2loc<-seqlevels(zoom$GR)
  hic.data.frame <- strawr::straw(norm=norm,
                                  fname=hicFile,
                                  chr1loc=chr1loc,
                                  chr2loc=chr2loc,
                                  unit="BP",
                                  bin=bin)
  return(hic.data.frame)
}

hic.data.frame <- strawr::straw("KR", "/path/to/file.hic", "11", "11", "BP", 10000)

hic.file <- "/home/malli/work/insulation/hichip/K4_hichip_q_10.hic"
Dump <- DFfromHIC(hicFile = hic.file,seq_style = "NCBI",zoom = zoomString,norm = "NONE")
GI <- hicDump2GI(Dump,chr=zoom$chr)
HiCO <- createHicPainteRObj(GI,mapType = "cMap")

zoom<-createZooms(zoomString)
kp <-plotKaryotype(zoom = zoom$GR)
kpHiC(kp,HiCO)

TAD<-regioneR::toGRanges(data.frame(chr="chr2",start=min(start(gInteractionsTocMap(GI))),end=max(start(gInteractionsTocMap(GI)))))
selGR<-reShapeTAD(TAD,gInteractionsTocMap(GI), xbin=10000)

zoomString<-"chr2:175,849,978-178,182,899"
XX<-xtabs(log(selGR$counts)~selGR$x1 +selGR$y)



mat<-as.matrix(XX)
image(log(as.matrix(XX)))



