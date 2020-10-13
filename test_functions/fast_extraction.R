library(strawr)
library(regioneR)
library(GenomicInteractions)
library(plyranges)
library(tidyverse)
library(karyoploteR)
library(RColorBrewer)
library("viridis")


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

#hic.data.frame <- strawr::straw("KR", "/path/to/file.hic", "11", "11", "BP", 10000)
zoom<-createZooms(zoomString)


hic.file <- "/home/rmalinverni/shares/BDbuschbeck/rmalinverni/07_Biodata/hic/ENCODE_hic_hepg2.chr2.hic"
Dump <- DFfromHIC(hicFile = hic.file,seq_style = "NCBI",zoom = zoom$zoom, norm = "KR")
GI <- hicDump2GI(Dump,chr=zoom$chr)
HiCO <- createHicPainteRObj(GI,mapType = "cMap")


kp <-plotKaryotype(zoom = zoom$GR)
kpHiC(kp,HiCO)



TAD<-regioneR::toGRanges(data.frame(chr="chr2",start=min(start(gInteractionsTocMap(GI))),end=max(start(gInteractionsTocMap(GI)))))
selGR<-reShapeTAD(TAD,gInteractionsTocMap(GI), xbin=1000)
fakeGenome<-regioneR::toGRanges(data.frame(chr="ChrU",start=1,end=1000))
XX<-as.data.frame(xtabs(log(selGR$counts)~selGR$x1 +selGR$y))
#XX<- XX %>% filter(Freq !=0)
kp<-plotKaryotype(genome = fakeGenome)

kpPoints(kp,
         x=as.numeric(as.character(XX$selGR.x1)),
         y = as.numeric(as.character(XX$selGR.y)),
         chr = "ChrU",
         col=alpha(alpha=rescale(log1p(XX$Freq)),col="red"),pch=18,cex=3)

image(log1p(as.matrix(as.data.frame.matrix(xtabs(log(selGR$counts)~selGR$x1 +selGR$y)))))

image(xtabs(log(selGR$counts)~selGR$x1 +selGR$y), axes=FALSE, zlim=c(-2,2), col=ramp(1000))

as.data.frame.matrix(xtabs(log(selGR$counts)~selGR$x1 +selGR$y))
colsR<-round(rescale((XX$Freq),to=c(0,100)))

ramp<-colorRampPalette(c("white","cadetblue","red","orange","yellow","yellow","yellow"))

TAD<-regioneR::toGRanges(data.frame(chr="chr2",start=min(start(gInteractionsTocMap(GI))),end=max(start(gInteractionsTocMap(GI)))))
selGR<-reShapeTAD(TAD,gInteractionsTocMap(GI), xbin=100)

zoomString<-"chr2:175,849,978-178,182,899"

library(lattice)

levelplot(xtabs(log1p(selGR$counts)~selGR$x1 +selGR$y),col.regions=ramp(100))


selGR<-reShapeTAD(TAD,gInteractionsTocMap(GI), xbin=500,yfactor=2)
rsc<-(xtabs(selGR$counts~selGR$x1 +selGR$y))
rsc<-as.matrix(as.data.frame.matrix(rsc))
levelplot(log1p(rsc),col.regions=ramp(1000))





df<-data.frame(C1=c(1,1,1,0,2),
               C2=c(2,2,2,3,3),
               C3=c(10,5,1,0,0))
XT1<-xtabs(df$C3~df$C1 + df$C2)

df2<-data.frame(C1=c(1,2,3,4,5),
               C2=c(1,2,3,4,5),
               C3=c(1,5,1,0,0))
XT2<-xtabs(df2$C3~df2$C1 + df2$C2)
mat<-as.matrix(XX)
image(log(as.matrix(XX)))
XT1<-as.matrix(as.data.frame.matrix(XT1))
XT2<-as.matrix(as.data.frame.matrix(XT2))

(XT2,XT1)

