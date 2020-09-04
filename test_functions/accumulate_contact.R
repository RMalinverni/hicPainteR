library("scales")
setwd("/mnt/data1/hicPainteR_data_test/")
save(mat,file="matrix_1k_k4_on_k4.RData")save(mat,file="matrix_1k_k4_on_k4.RData")
save(matmacr,file="matrix_all_macro_on_macro.RData")


matmacr

GIhichip_macro
GI_hic_TAD_10K
GI_hic_TAD_5K<-arrowHead2GI(hic_TAD_5k_file)

bed_RAD_21<-toGRanges(read.delim("/home/rmalinverni/shares/BDbuschbeck/rmalinverni/01_LAB/David/ChIPseq/analysis/MACS2/RAD21_R1RD_summits.bed",
                       header=FALSE))
bed_RAD_21<- filterChromosomes(bed_RAD_21)

H3K4_peaks<-toGRanges(read.delim(gzfile("/home/rmalinverni/Downloads/wgEncodeBroadHistoneHepg2H3k4me3StdPk.broadPeak.gz"),header = FALSE))

Peaks<-H3K4_peaks
Peaks<-resize(Peaks,width=10,fix = "center")
Peaks<-extendRegions(Peaks,extend.start = 3e5,extend.end = 3e5)
Peaks<-filterChromosomes(Peaks,keep.chr = "chr2")
Peaks<-Peaks[start(Peaks)>0]


RAD21_TADlike<-toGRanges(data.frame(chr=as.character(seqnames(bed_RAD_21)),
                                    start=start(bed_RAD_21)-5e5,
                                    end=end(bed_RAD_21)+5e5))

RAD21_TADlike<-filterChromosomes(RAD21_TADlike,keep.chr = "chr2")
RAD21_TADlike<-RAD21_TADlike[-(1:3)]

GI_hic_TAD_10K
genes.pos<-toGRanges(read.delim("/mnt/data1/David_STAR/pos_DE_genes.tsv"))
genes.neg<-toGRanges(read.delim("/mnt/data1/David_STAR/neg_DE_genes.tsv"))
genesDEG<-c(genes.pos,genes.pos)


DEG_TADlike<-toGRanges(data.frame(as.character(seqnames(genes.pos)),start(genes.pos)-5e5,start(genes.pos)+5e5))

DEG_TADlike<-DEG_TADlike[seqnames(DEG_TADlike)=="chr2"]


GR_TADS<-toGRanges(data.frame("chr2",start(anchors(GI_hic_TAD_10K,"first")),start(anchors(GI_hic_TAD_10K,"second"))))
GI_hic_TAD_pos<-GI_hic_TAD_10K[overlapRegions(GR_TADS,genes.pos,only.boolean = TRUE)]
GI_hic_TAD_neg<-GI_hic_TAD_10K[overlapRegions(GR_TADS,genes.neg,only.boolean = TRUE)]
GI_hic_TAD_DEG<-c(GI_hic_TAD_pos,GI_hic_TAD_neg)

a<-GI_hic_TAD_10K[107]
GR_hichip<-toGRanges(data.frame("chr2",start(anchors(GIhichip_macro,"first")),start(anchors(GIhichip_macro,"second"))))
GR_hichip<- GR_hichip %>% mutate(counts=GIhichip_macro$counts)
GR<-toGRanges(data.frame("chr2",start(anchors(a,"first")),start(anchors(a,"second"))))

#controllare normalizzazione per length(TAD)
reShapeTAD<-function(TAD,GR_test,xbin=10000,ybin=2){
  selGR<- GR_test %>% filter(start(GR_test) >= start(TAD), end(GR_test)<= end(TAD))
  lengthTAD<-end(TAD)-start(TAD)
  selGR<- selGR %>% mutate (y=round((end(selGR)-start(selGR))/lengthTAD,digits = ybin))
  selGR<- selGR %>% mutate (x=((end(selGR)-start(selGR))/2)+start(selGR))
  selGR<- selGR %>% mutate (x1=round(scales::rescale(selGR$x, to=c(1,xbin)),digits = 0))
  return(selGR)
}

randomPeaks<-randomizeRegions(Peaks,per.chromosome = TRUE, allow.overlaps = TRUE)

  bin<-1000
  GI_TAD<-randomPeaks
  GI_test<-GIhichip_K4
  smp<-sample(length(GI_TAD),1000)
  cn<-0

  mat<-matrix(nrow = 101,ncol=1000,data = 0)
  colnames(mat)<-seq(1,1000)
  rownames(mat)<-seq(0.0,1,0.01)
  for( i in smp){
    cn<-cn+1
    print(cn)
    GR<-GI_TAD[i]
    GR_hichip<-toGRanges(data.frame("chr2",start(anchors(GI_test,"first")),start(anchors(GI_test,"second"))))
    GR_hichip<- GR_hichip %>% mutate(counts=GI_test$counts)
    GR<-toGRanges(data.frame("chr2",start(anchors(a,"first")),start(anchors(a,"second"))))

    selGR<-reShapeTAD(TAD=GR,GR_test = GR_hichip,xbin = bin)

    for(j in 1:length(selGR)){
      mat[as.character(selGR$y[j]),as.character(selGR$x1[j])]<-mat[as.character(selGR$y[j]),as.character(selGR$x1[j])] + selGR$counts[j]
    }
  }

matK4rand<-mat
save(matK4rand,file="matrix_1K_rnd_k4_on_k4.RData")
transMat<-as.data.frame(as.table(mat))
transMat<-transMat %>% arrange(Freq)
transMat<-transMat %>% filter(Freq !=0)
transMat<-transMat %>% mutate(scaleFreq=rescale(Freq))
transMat<-transMat %>% mutate(Var1=as.numeric(as.character(Var1)))
transMat<-transMat %>% mutate(Var2=as.numeric(as.character(Var2)))
transMat<-transMat %>% filter(Var1 >= 0.01)
transMat<-transMat %>% mutate(scaleFreq=rescale(Freq))

fakeGen<-toGRanges(data.frame(chr="ChrU",1,bin))
kp <- plotKaryotype(genome=fakeGen)
Ys<-as.numeric(as.character(transMat$Var1))
Xs<-as.numeric(as.character(transMat$Var2))
kpPoints(kp,chr="ChrU",x=Xs,y=Ys,col=alpha(alpha=1*transMat$scaleFreq,colour = "red"),cex=1,pch=19)


lim<-20000
mat<-tmat
transMat<-as.data.frame(as.table(mat))
transMat$Freq[transMat$Freq > lim]<-lim
transMat<-transMat %>% arrange(Freq)
transMat<-transMat %>% filter(Freq !=0)
transMat<-transMat %>% mutate(scaleFreq=rescale(Freq))
transMat<-transMat %>% mutate(Var1=as.numeric(as.character(Var1)))
transMat<-transMat %>% mutate(Var2=as.numeric(as.character(Var2)))
#transMat<-transMat %>% filter(Var1 >= 0.04)
transMat<-transMat %>% mutate(scaleLogFreq=rescale(log(Freq)))


ramp<-colorRampPalette(colors=c("white","lightblue","red","orange","yellow"),interpolate="linear",bias=2)(100)
CCols<-ramp[round(transMat$scaleFreq*100,digits = 0)+1]
fakeGen<-toGRanges(data.frame(chr="ChrU",1,bin))
kp <- plotKaryotype(genome=fakeGen)
Ys<-as.numeric(as.character(transMat$Var1))
Xs<-as.numeric(as.character(transMat$Var2))
kpPoints(kp,chr="ChrU",x=Xs,y=Ys,col=CCols,cex=1.4,pch=18)




load(file="matrix_1k_k4_on_k4.RData")
matK4<-mat
mat<-matK4
resmatK4<-rescale(matK4)
resmatK4rand<-rescale(matK4rand)
mat<-resmatK4/resmatK4rand





reShapeTAD<-function(TAD,GR_test,xbin=10000,yfactor=2){
  selGR<- GR_test %>% filter(start(GR_test) >= start(TAD), end(GR_test)<= end(TAD))
  lengthTAD<-end(TAD)-start(TAD)
  selGR<- selGR %>% mutate (y=round((end(selGR)-start(selGR))/lengthTAD,digits = yfactor))
  selGR<- selGR %>% mutate (x=((end(selGR)-start(selGR))/2)+start(selGR))
  selGR<- selGR %>% mutate (x1=round(scales::rescale(selGR$x, to=c(1,xbin)),digits = 0))
  return(selGR)
}

bin<-1000
GI_TAD<-randomPeaks
GI_test<-GIhichip_K4
smp<-sample(length(GI_TAD),1000)
if (verbose==TRUE){
  cn<-cn+1
  print(cn)
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

kp_TransMat<-function(mat,lim=40000,ramp=NULL,autoAlpha=FALSE,filter.at=1){

  if ( is.null(ramp)){
    ramp<-colorRampPalette(colors=c("white","lightblue","red","orange","yellow"),interpolate="linear",bias=2)(100)
  }

  transMat<-as.data.frame(as.table(mat))
  transMat$Freq[transMat$Freq > lim]<-lim
  transMat<-transMat %>% arrange(Freq)
  transMat<-transMat %>% filter(Freq !=0)
  transMat<-transMat %>% mutate(scaleFreq=rescale(Freq))
  transMat<-transMat %>% filter(scaleFreq <= filter.at)
  transMat<-transMat %>% mutate(Var1=as.numeric(as.character(Var1)))
  transMat<-transMat %>% mutate(Var2=as.numeric(as.character(Var2)))
  transMat<-transMat %>% filter(Var1 >= 0.04)
  transMat<-transMat %>% mutate(scaleLogFreq=rescale(log(Freq)))
  CCols<-ramp[round(transMat$scaleFreq*100,digits = 0)+1]
  Ys<-as.numeric(as.character(transMat$Var1))
  Xs<-as.numeric(as.character(transMat$Var2))
  if (autoAlpha == FALSE){
    kpPoints(kp,chr="ChrU",x=Xs,y=Ys,col=CCols,cex=1.4,pch=18)
  }else{
    kpPoints(kp,chr="ChrU",x=Xs,y=Ys,col=alpha(alpha=1*transMat$scaleFreq,CCols),cex=1.4,pch=18)
  }
}

rndPeaks<-randomizeRegions(Peaks,per.chromosome = TRUE)
bin<-1000

registerDoMC(cores=16)

TADs<-toGRanges(data.frame(
  chr=as.vector(seqnames(anchors(GI_hic_TAD_10K,type="first"))),
  start=start(anchors(GI_hic_TAD_10K,type="first")),
  end=end(anchors(GI_hic_TAD_10K,type="second"))))
rndTADs<-randomizeRegions(TADs,per.chromosome = TRUE,genome = "hg19")



macroOnTAD<-comulativeMatrix(GR=TADs, GI = GIhichip_macro,verbose = T)
k4OnTAD<-comulativeMatrix(GR=TADs, GI = GIhichip_K4,verbose = T)
macroOnTADrnd<-comulativeMatrix(GR=rndTADs, GI = GIhichip_K4,verbose = T)

fakeGen<-toGRanges(data.frame(chr="ChrU",1,bin))
kp <- plotKaryotype(genome=fakeGen)
kp_TransMat(macroOnTAD,lim = 300)
kp <- plotKaryotype(genome=fakeGen)
kp_TransMat(k4OnTAD,lim = 10000)

load("mat_k4_on_k4_all_chr2.RData")
load("matrix_1K_rnd_k4_on_k4.RData")
kp <- plotKaryotype(genome=fakeGen)
kp_TransMat(matK4rand,lim = 5000,autoAlpha = TRUE)

kp <- plotKaryotype(genome=fakeGen)
kp_TransMat(matK4rand,lim = 3000,autoAlpha = TRUE,filter.at = 1)

hist(rescale(tmat))
hist(rescale(matK4rand))


GR=Peaks[1:10]
chrs<-as.vector(seqnames(anchors(GI,type="first")))

GR_3D
bin<-1000
yfactor<-2
GI<-GIhichip_K4
mat<-matrix(nrow = (10^yfactor)+1,ncol=bin,data = 0)
colnames(mat)<-seq(1,bin)
rownames(mat)<-seq(0.0,1,10^-yfactor)
chrs<-as.vector(seqnames(anchors(GI,type="first")))
GR_3D<-toGRanges(data.frame(chrs,start(anchors(GI,"first")),start(anchors(GI,"second"))))
GR_3D<- GR_3D %>% mutate(counts=GI$counts)


test<-foreach( i = 1:3) %do%{
  mat<-Function1(GR,GR_3D,bin,yfactor,i)
}

matCreate<-function(selGR,mat){
  matX<-as.character(selGR$y[j])
  matY<-as.character(selGR$x1[j])
  mat[matX,matY]<-mat[matX,matY] + selGR$counts[j]
  return(mat)
}

Function1<-function(GR,GR_3D,bin,yfactor,i){
  GRi<-GR[i]
  selGR<-reShapeTAD(TAD=GRi,GR_test = GR_3D,xbin = bin, yfactor=yfactor)
  for(j in 1:length(selGR)){
    mat<-matCreate(selGR = selGR, mat= mat)
  }
  return(mat)
}



test<-foreach( c = 1:3 , .combine = '+') %dopar%{
  a<-sqrt(c)
  b<-sqrt(a)}



