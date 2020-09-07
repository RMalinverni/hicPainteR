library("karyoploteR")
library("GenomicInteractions")
library("plyranges")
library("ggplot2")

#for private
setwd("~/work/insulation/script/hicPainteR/Test_data_hicPainter")

#for work
setwd("/mnt/data1/Test_data_hicPainter")

hic_file<-"HIC_hepg2_chr2_dump.txt"
hic_loops_file<-"HIC_hepg2_chr2_merged_loops.bedpe"
hic_TAD_10k_file<-"HIC_hepg2_chr2_TADs_10K.bedpe"
hic_TAD_5K_file<-"HIC_hepg2_chr2_TADs_5K.bedpe"
hichip_k4_file<-"HiChIP_hepg2_chr2_h3k4me3_dump.txt"
hichip_macroH2A_file<-"HiChIP_hepg2_chr2_macroH2A_dump.txt"

GIhic<-hicDump2GI(hic_file,chr = "chr2")
GI_hic_loops<-hiccups2GI(hic_loops_file)
GI_hic_TAD_10K<-arrowHead2GI(hic_TAD_10k_file)
GIhichip_K4<-hicDump2GI(hichip_k4_file,chr = "chr2")
GIhichip_macro<-hicDump2GI(hichip_macroH2A_file,chr = "chr2")

zoom1<-'chr2:156,000,000-163,000,000'

hPhic<-createHicPainteRObj(Name = "hic_HepG2_chr2",mapType = 'cMap',contact = GIhic,
                           zoom=zoom1,log=TRUE,alpha=0.001,use_ramp = TRUE)

hPhichip_K4<-createHicPainteRObj(Name = "hichip_HepG2_K4_chr2", mapType = 'cMap',contact = GIhichip_K4,
                              zoom=zoom1,exp=TRUE,alpha=0.01,)

hPhichip_macro<-createHicPainteRObj(Name = "hichip_HepG2_MacroH2A_chr2", mapType = 'cMap',contact = GIhichip_macro,
                                 zoom=zoom1,alpha=0.009,pCol="blue",exp=TRUE,cexP = 2)

hP_loops<-createHicPainteRObj(Name = "hic_HepG2_loops_chr2",mapType = 'Loop',contact = GI_hic_loops,
                               zoom=zoom1, pCol="darkgreen")

hP_TADs_10K<-createHicPainteRObj(Name = "hic_HepG2_loops_chr2",mapType = 'TAD',contact = GI_hic_TAD_10K,
                                 zoom=zoom1,pCol="black")


SS<-158e6
EE<-161e6
chr="chr2"
zoom<-regioneR::toGRanges(data.frame(chr=chr,start=SS,end=EE))


pp<-getDefaultPlotParams(plot.type = 1)
pp$leftmargin<-0.15
pp$topmargin <- 15
pp$bottommargin <-15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

ntracks=3

hPhic1<-hPhic

hPhic1$cMap<-hPhic1$cMap[hPhic1$cMap$Ynorm <= 0.05]

setwd("/home/malli/work/insulation/Images/expression")
pdf(file="locus_with_point_loops_genes.pdf",width = 15, height = 9)
kp<-plotKaryotype(chromosomes = "chr2", zoom=zoom,plot.params = pp, main= "locus_extended_with_genes_C")

kpAddBaseNumbers(kp, tick.dist = 10e6, minor.tick.dist = 1e6,
                 add.units = TRUE, cex=1.2, tick.len = 3)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = FALSE,
                                    plot.transcripts.structure = FALSE)
genes.data <- addGeneNames(genes.data)
GRgenes<-genes.data$genes
selGRgenes<-GRgenes[which(GRgenes$name %in% selected_genes$SYMBOL)]


at <- autotrack(c(1,6),ntracks*2,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(c(1,2),ntracks*2,r0=0,r1=1,margin=0.5)
kpHiC(kp,hicPobj = hPhic1)
kpLoops(kp,hicPobj = hP_loops, loopType = "Points")

kpAxis(kp, ymin=0, ymax=2, numticks = 2, r0=at$r0, r1=at$r1,cex=0.7)
kpAddLabels(kp, labels = hPhic1$param$Name, r0=at$r0, r1=at$r1,
            cex=0.7, label.margin = 0.015)

hPhichip_K4<-repaintHicPainteRobj(hPhichip_K4,alpha = 0.03)
at <- autotrack(c(3,6),ntracks*2,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(c(3,4),ntracks*2,r0=0,r1=1,margin=0.5)
kpHiC(kp,hicPobj = hPhichip_K4)

kpAxis(kp, ymin=0, ymax=2, numticks = 2, r0=at$r0, r1=at$r1,cex=0.7)
kpAddLabels(kp, labels = hPhichip_K4$param$Name, r0=at$r0, r1=at$r1,
            cex=0.7, label.margin = 0.015)
#kpTAD(kp,hicPobj = hP_TADs_10K)

hPhichip_macro<-repaintHicPainteRobj(hPhichip_macro,alpha = 0.03,cexP=2.4,use_ramp = TRUE)
at <- autotrack(c(4,6),ntracks*2,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(c(5,6),ntracks*2,r0=0,r1=1,margin=0.5)
kpHiC(kp,hicPobj = hPhichip_macro)
kpAxis(kp, ymin=0, ymax=2, numticks = 2, r0=at$r0, r1=at$r1,cex=0.7)
kpAddLabels(kp, labels = hPhichip_macro$param$Name, r0=at$r0, r1=at$r1,
            cex=0.7, label.margin = 0.015)
#kpTAD(kp,hicPobj = hP_TADs_10K)

#at <- autotrack(c(4),ntracks*2,r0=0,r1=1,margin=0.1)
#kp <- kpPlotBigWig(kp, data=as.character(bigwigs[[1]]), ymax="visible.region",
                  # r0=at$r0, r1=at$r1, col = col_vector[1])
#at <- autotrack(c(4),ntracks*2,r0=0,r1=1,margin=0.1)
#kp <- kpPlotBigWig(kp, data=as.character(bigwigs[[2]]), ymax="visible.region",
                   #r0=at$r0, r1=at$r1, col = col_vector[2])
at <- autotrack(c(3,4),ntracks*2,r0=0,r1=1,margin=0.5)
kpPlotMarkers(kp, data=GRgenes, labels=GRgenes$name,cex=0.7,r0=at$r0,r1=at$r1)
at <- autotrack(c(5,6),ntracks*2,r0=0,r1=1,margin=0.5)
kpPlotMarkers(kp, data=selGRgenes, labels=selGRgenes$name,cex=0.7,r0=at$r0,r1=at$r1)

dev.off()

pdf(file="locus_with_chipseq_C.pdf",width = 11, height = 13)
kp<-plotKaryotype(chromosomes = "chr2", zoom=zoom,plot.params = pp, main= "locus_with_genes_C")

kpAddBaseNumbers(kp, tick.dist = 1e6, minor.tick.dist = 1e5,
                 add.units = TRUE, cex=1.2, tick.len = 3)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = FALSE,
                                    plot.transcripts.structure = FALSE)
genes.data <- addGeneNames(genes.data)
GRgenes<-genes.data$genes
selGRgenes<-GRgenes[which(GRgenes$name %in% selected_genes$SYMBOL)]


at <- autotrack(c(1,11),11,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(c(1,3),11,r0=0,r1=1,margin=0.5)
kpHiC(kp,hicPobj = hPhic1)
#kpLoops(kp,hicPobj = hP_loops, loopType = "Points")
#kpAxis(kp, ymin=0, ymax=2, numticks = 2, r0=at$r0, r1=at$r1,cex=0.7)
kpAddLabels(kp, labels = hPhic1$param$Name, r0=at$r0, r1=at$r1,
            cex=0.7, label.margin = 0.015)


hPhichip_K4<-repaintHicPainteRobj(hPhichip_K4,alpha = 0.03)
at <- autotrack(c(3,11),11,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(c(3,5),11,r0=0,r1=1,margin=0.5)
kpHiC(kp,hicPobj = hPhichip_K4)

#kpAxis(kp, ymin=0, ymax=2, numticks = 2, r0=at$r0, r1=at$r1,cex=0.7)
kpAddLabels(kp, labels = hPhichip_K4$param$Name, r0=at$r0, r1=at$r1,
            cex=0.7, label.margin = 0.015)
#kpTAD(kp,hicPobj = hP_TADs_10K)

hPhichip_macro<-repaintHicPainteRobj(hPhichip_macro,alpha = 0.03,cexP=2.4,use_ramp = TRUE)
at <- autotrack(c(5,11),11,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(c(5,7),11,r0=0,r1=1,margin=0.5)
kpHiC(kp,hicPobj = hPhichip_macro)
#kpAxis(kp, ymin=0, ymax=2, numticks = 2, r0=at$r0, r1=at$r1,cex=0.7)
kpAddLabels(kp, labels = hPhichip_macro$param$Name, r0=at$r0, r1=at$r1,cex=0.7, label.margin = 0.015)

at <- autotrack(c(7,11),11,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)

at <- autotrack(c(7),11,r0=0,r1=1,margin=0.1)
kp <- kpPlotBigWig(kp, data=as.character(bigwigs[[1]]), ymax="visible.region", r0=at$r0, r1=at$r1, col = col_vector[1])
kpAddLabels(kp, labels = names(bigwigs)[1], r0=at$r0, r1=at$r1, cex=0.7, label.margin = 0.015)
at <- autotrack(c(8),11,r0=0,r1=1,margin=0.1)
kp <- kpPlotBigWig(kp, data=as.character(bigwigs[[2]]), ymax="visible.region", r0=at$r0, r1=at$r1, col = col_vector[2])
kpAddLabels(kp, labels = names(bigwigs)[2], r0=at$r0, r1=at$r1, cex=0.7, label.margin = 0.015)
at <- autotrack(c(9),11,r0=0,r1=1,margin=0.1)
kp <- kpPlotBigWig(kp, data=as.character(bigwigs[[3]]), ymax="visible.region", r0=at$r0, r1=at$r1, col = col_vector[3])
kpAddLabels(kp, labels = names(bigwigs)[3], r0=at$r0, r1=at$r1, cex=0.7, label.margin = 0.015)
at <- autotrack(c(10),11,r0=0,r1=1,margin=0.1)
kp <- kpPlotBigWig(kp, data=as.character(bigwigs[[4]]), ymax="visible.region", r0=at$r0, r1=at$r1, col = col_vector[4])
kpAddLabels(kp, labels = names(bigwigs)[4], r0=at$r0, r1=at$r1, cex=0.7, label.margin = 0.015)
at <- autotrack(c(11),11,r0=0,r1=1,margin=0.1)
kp <- kpPlotBigWig(kp, data=as.character(bigwigs[[5]]), ymax="visible.region", r0=at$r0, r1=at$r1, col = col_vector[5])
kpAddLabels(kp, labels = names(bigwigs)[5], r0=at$r0, r1=at$r1, cex=0.7, label.margin = 0.015)


at <- autotrack(5,11,r0=0,r1=1,margin=0.5)
kpPlotMarkers(kp, data=selGRgenes, labels=selGRgenes$name,cex=0.7,r0=at$r0,r1=at$r1)

dev.off()


kpPlotMarkers(kp, data=GRgenes, labels=GRgenes$name,cex=0.7,r0=0.0, r1=0.13)



map<-hPhichip_macro$cMap
n=10000
Xm<-10000
Ym<-5


crossCalc<-function(map,n,Xm,Ym){
  numb<-map$edge[n]
  numbY<-map$Ys[n]
  map %>% filter (Ys >= Ys-Ym & Ys <= Ys+Ym) -> map
  map %>% filter((edge >= numb-Xm) & (edge <= numb+Xm)) %>% filter ((Ys >= numbY-Ym) &(Ys <=numbY +Ym)) -> test1
  #print(length(test1))
  res<- sum(test1$counts)
  return(res)
}

vecC<-vector()
for( i in 1:length(map)){
  vecC[i]<-crossCalc(map=map,n=i,Xm=10000,Ym=5)
}












#crossCalc(test,1,20000,.1)
#param <- bpparam()
#bpworkers(param) <- 2
#vecC<-unlist(mclapply(1:10000,FUN=crossCalc,map=test,Xm=20000,Ym=.1,mc.cores=10))



