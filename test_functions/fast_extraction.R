library(strawr)
library(regioneR)
library(GenomicInteractions)
library(plyranges)
library(tidyverse)
library(karyoploteR)
library(RColorBrewer)
library("viridis")

library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("GenomicFeatures")


a
refgene <-makeTxDbFromUCSC(genome="hg19",
                                    tablename="refGene")
transcripts <- transcripts(refgene, columns=c("tx_id", "tx_name"))

tss <- resize(transcripts, width=1, fix='start')
tss<-tss[sample(length(tss),100)]
tss<-filterChromosomes(tss)
tadtss<-extendRegions(tss,extend.start = 400000,extend.end = 400000)
hg19<-filterChromosomes(getGenome("hg19"))
tadtss2<-commonRegions(tadtss[1],hg19,ignore.strand=TRUE)

tadtss<-intersect(hg19,tadtss,ignore.strand=TRUE)




TAD<-
zoomString<-"chr2:175,849,978-178,182,899"

createZooms<-function(zoom){

  if (class(zoom)=="GRanges"){
    GR<-zoom
    chr<-seqlevels(zoom)[1]
    start<-start(zoom)[1]
    end<-end(zoom)[1]
    zoomString<-paste0(chr,":",start,"-",end)
  }

  if (is.character(zoom)){
    zoomString<-gsub(",","",zoom)
    chr<-strsplit(zoom,":")[[1]][1]
    range<-strsplit(zoom,":")[[1]][2]
    start<-strsplit(range,"-")[[1]][1]
    end<-strsplit(range,"-")[[1]][2]
    GR<-regioneR::toGRanges(data.frame(chr=chr,start=start,end=end))
  }

  zoomObj<-list(
    zoomString=zoomString,
    chr=chr,
    start=start,
    end=end,
    GR=GR)

  return(zoomObj)
}

reShapeTAD<-function(TAD,GR_test,xbin=100){
  yfactor<-(-log10(1/xbin))
  selGR<- GR_test %>% filter(start(GR_test) >= start(TAD), end(GR_test)<= end(TAD))
  lengthTAD<-end(TAD)-start(TAD)
  selGR<- selGR %>% mutate (y=round((end(selGR)-start(selGR))/lengthTAD,digits = yfactor))
  selGR<- selGR %>% mutate (x=((end(selGR)-start(selGR))/2)+start(selGR))

  if(length(selGR)!=0){
    selGR<- selGR %>% mutate (x1=round(scales::rescale(selGR$x, to=c(1,xbin)),digits = 0))

    AllBinX<-c(1:xbin)
    AllBinY<-seq(from = 0,to = 1,by=1/10^yfactor)

    diffY<-setdiff(AllBinY,unique(selGR$y))
    diffX<-setdiff(AllBinX,unique(selGR$x1))

    if (length(diffY)!= 0){
      fakeGR<-selGR[1:length(diffY)]
      fakeGR$y<-diffY
      fakeGR$counts<-0
      fakeGR$x1<-1
      fakeGR$Ynorm<-0
      fakeGR$Ys<-0
      selGR<-c(selGR,fakeGR)
    }

    if (length(diffX)!= 0){
      fakeGR<-selGR[1:length(diffX)]
      fakeGR$x1<-diffX
      fakeGR$counts<-0
      fakeGR$y<-0
      fakeGR$Ynorm<-0
      fakeGR$Ys<-0
      selGR<-c(selGR,fakeGR)
    }
  }else{

    nelementGR<-max(length(seq(from = 0,to = 1,by=1/10^yfactor)),xbin)
    selGR<-toGRanges(data.frame(chr=rep(seqlevels(TAD)[1],xbin),
               start=rep(start(TAD)[1],nelementGR),
               end=rep(end(TAD)[1],nelementGR),
               counts=rep(0,nelementGR),
               Ynorm=rep(0,nelementGR),
               Ys=rep(0,nelementGR),
               y=seq(from = 0,to = 1,by=1/10^yfactor),
               x=rep(0,nelementGR),
               x1=c(1,rep(1:xbin))))
    selGR$x1[selGR$x1 > xbin]<-xbin
  }


  return(selGR)
}

matTest<-xtabs((selGR$counts)~selGR$x1 +selGR$y)

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

zoom<-zoom$zoomString
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


hic.file <- "/home/malli/work/insulation/hichip/hicup_analisis/HiChIP_K4_R1_2.hic"
Dump <- DFfromHIC(hicFile = hic.file,seq_style = "NCBI",zoom = zoom$zoom, norm = "KR")
GI <- hicDump2GI(Dump,chr=zoom$chr)
HiCO <- createHicPainteRObj(GI,mapType = "cMap")


kp <-plotKaryotype(zoom = zoom$GR)
kpHiC(kp,HiCO)



#TAD<-regioneR::toGRanges(data.frame(chr="chr2",start=min(start(gInteractionsTocMap(GI))),end=max(start(gInteractionsTocMap(GI)))))
selGR<-reShapeTAD(zoom$GR,gInteractionsTocMap(GI), xbin=100)
fakeGenome<-regioneR::toGRanges(data.frame(chr="ChrU",start=1,end=100))
XX<-as.data.frame(xtabs(log(selGR$counts)~selGR$x1 +selGR$y))
XX<- XX %>% filter(Freq !=0)
kp<-plotKaryotype(genome = fakeGenome)

kpPoints(kp,
         x=as.numeric(as.character(XX$selGR.x1)),
         y = as.numeric(as.character(XX$selGR.y)),
         chr = "ChrU",
         col=alpha(alpha=rescale(XX$Freq),col="red"),pch=18,cex=1)

image(log1p(as.matrix(as.data.frame.matrix(xtabs(log(selGR$counts)~selGR$x1 +selGR$y)))))

image(xtabs(log(selGR$counts)~selGR$x1 +selGR$y), axes=FALSE, zlim=c(-2,2), col=ramp(1000))

as.data.frame.matrix(xtabs(log(selGR$counts)~selGR$x1 +selGR$y))
colsR<-round(rescale((XX$Freq),to=c(0,100)))

ramp<-colorRampPalette(c("white","cadetblue","red","orange","yellow","yellow","yellow"))

TAD<-regioneR::toGRanges(data.frame(chr="chr2",start=min(start(gInteractionsTocMap(GI))),end=max(start(gInteractionsTocMap(GI)))))
selGR<-reShapeTAD(TAD,gInteractionsTocMap(GI), xbin=100)

zoomString<-"chr2:175,849,978-178,182,899"

library(lattice)

selGR<-reShapeTAD(zoom$GR,gInteractionsTocMap(GI), xbin=100,yfactor = 2)
fakeGenome<-regioneR::toGRanges(data.frame(chr="ChrU",start=1,end=100))
levelplot(xtabs(log1p(selGR$counts)~selGR$x1 +selGR$y),col.regions=ramp(100))


matTest<-xtabs((selGR$counts)~selGR$x1 +selGR$y)

selGR<-reShapeTAD(TAD,gInteractionsTocMap(GI), xbin=500,yfactor=2)
rsc<-(xtabs(selGR$counts~selGR$x1 +selGR$y))
rsc<-as.matrix(as.data.frame.matrix(rsc))
levelplot(log1p(matSum),col.regions=ramp(1000))

TADs<-read.delim("/home/malli/work/insulation/TADs/10000_blocks.bedpe")
TADs<-toGRanges(TADs[-1,])
seqlevelsStyle(TADs)<-"UCSC"
TADs<-intersect(TADs,hg19,ignore.strand=TRUE)
hic.file <- "/home/malli/work/insulation/hichip/hicup_analisis/HiChIP_Macro_1st_R1_2.hicup.hic"
hic.file <- "/home/malli/work/insulation/hichip/hicup_analisis/HiChIP_K4_R1_2.hic"
matList<-list()

newReg<-sample((tadtss),60)
newReg
for(i in 1:length(newReg)){
  print(i)
  zoom<-createZooms(newReg[i])
  zoom
  Dump <- DFfromHIC(hicFile = hic.file,seq_style = "NCBI",zoom = zoom$zoomString, norm = "NONE")
  GI <- hicDump2GI(Dump,chr=zoom$chr)
  selGR<-reShapeDomain(zoom$GR,gInteractionsTocMap(GI), xbin=100)
  matTest<-xtabs(selGR$counts~selGR$x1 +selGR$y)
  # if (i == 1){
  #   matSum<-matTest
  # }else{
  #   matSum<-matSum+matTest
  # }
  matList[[i]]<-matTest
}


levelplot(log1p(matSum),col.regions=ramp(1000))

require(purrr)
x = matList
v = reduce(x, `+`) / length(x)
image(v)
levelplot(log1p(v),col.regions=ramp(1000))
levelplot((log1p(v)^2),col.regions=ramp(1000))

