library("scales")
GIhichip_macro
GI_hic_TAD_10K
GI_hic_TAD_5K<-arrowHead2GI(hic_TAD_5k_file)

genes.pos<-toGRanges(read.delim("/home/malli/work/insulation/David_Star/pos_DE_genes.tsv"))
genes.neg<-toGRanges(read.delim("/home/malli/work/insulation/David_Star/neg_DE_genes.tsv"))

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


bin<-1000
fakeGen<-toGRanges(data.frame(chr="ChrU",1,bin))
kp <- plotKaryotype(genome=fakeGen)

smp<-sample(length(GI_hic_TAD_10K),50)
cn<-0

mat<-matrix(nrow = 101,ncol=1000,data = 0)
colnames(mat)<-seq(1,1000)
rownames(mat)<-seq(0.0,1,0.01)

for( i in 1:length(GI_hic_TAD_DEG)){
  cn<-cn+1
  print(cn)
  a<-GI_hic_TAD_DEG[i]
  GR_hichip<-toGRanges(data.frame("chr2",start(anchors(GIhichip_macro,"first")),start(anchors(GIhichip_macro,"second"))))
  GR_hichip<- GR_hichip %>% mutate(counts=GIhichip_macro$counts)
  GR<-toGRanges(data.frame("chr2",start(anchors(a,"first")),start(anchors(a,"second"))))

  selGR<-reShapeTAD(TAD=GR,GR_test = GR_hichip,xbin = bin)

  for(j in 1:length(selGR)){
    mat[as.character(selGR$y[j]),as.character(selGR$x1[j])]<-mat[as.character(selGR$y[j]),as.character(selGR$x1[j])] + selGR$counts[j]
  }
}

mat
transMat<-as.data.frame(as.table(mat))
transMat<-transMat %>% filter(Freq !=0)
transMat<-transMat %>% mutate(scaleFreq=rescale(Freq))
transMat<-transMat %>% mutate(Var1=as.numeric(as.character(Var1)))
transMat<-transMat %>% mutate(Var2=as.numeric(as.character(Var2)))
transMat<-transMat %>% filter(Var1 >= 0.03)
transMat<-transMat %>% mutate(scaleFreq=rescale(Freq))

fakeGen<-toGRanges(data.frame(chr="ChrU",1,bin))
kp <- plotKaryotype(genome=fakeGen)
Ys<-as.numeric(as.character(transMat$Var1))
Xs<-as.numeric(as.character(transMat$Var2))

library(RColorBrewer)
ramp<-colorRampPalette(c("white","red","yellow"))(100)
kpPoints(kp,chr="ChrU",x=Xs,y=Ys,col=alpha(alpha=1*transMat$scaleFreq,colour = "red"),cex=2,pch=".")

CCols<-ramp[round(transMat$scaleFreq*100,digits = 0)]
kpPoints(kp,chr="ChrU",x=Xs,y=Ys,col=CCols,cex=1.4,pch=18)










mat<-matrix(nrow = 100,ncol=100,data = 0)
colnames(mat)<-seq(1,100)
rownames(mat)<-seq(0.0,.99,0.01)

for(i in 1:length(selGR)){
    mat[as.character(selGR$y[i]),as.character(selGR$x1[i])]<-mat[as.character(selGR$y[i]),as.character(selGR$x1[i])] + selGR$counts[i]
}


mat["0.04","1"]<-20

selGR<-reShapeTAD(TAD=GR,GR_test = GR_hichip,xbin = bin)


file <- "http://www.sr.bham.ac.uk/~ajrs/papers/sanderson09/sanderson09_table2.txt"
a <- read.table(file, header=TRUE, sep="|")


table(a$cctype)
