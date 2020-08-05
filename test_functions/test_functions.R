library("karyoploteR")
library("GenomicInteractions")
library("plyranges")
library("ggplot2")

#for private
setwd("~/work/insulation/script/hicPainteR/Test_data_hicPainter")
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

zoom1<-'chr2:100,000,000-105,000,000'

hPhic<-createHicPainteRObj(Name = "hic_HepG2_chr2",mapType = 'cMap',contact = GIhic,
                           zoom=zoom1)

hPhichip_K4<-createHicPainteRObj(Name = "hichip_HepG2_K4_chr2", mapType = 'cMap',contact = GIhichip_K4,
                              zoom=zoom1)

hPhichip_macro<-createHicPainteRObj(Name = "hichip_HepG2_K4_chr2", mapType = 'cMap',contact = GIhichip_macro,
                                 zoom=zoom1)

hP_loops<-createHicPainteRObj(Name = "hic_HepG2_loops_chr2",mapType = 'Loop',contact = GI_hic_loops,
                               zoom=zoom1)

hP_TADs_10K<-createHicPainteRObj(Name = "hic_HepG2_loops_chr2",mapType = 'TAD',contact = GI_hic_TAD_10K,
                                 zoom=zoom1)


SS<-100e6
EE<-110e6
chr="chr2"
zoom<-regioneR::toGRanges(data.frame(chr=chr,start=SS,end=EE))


pp<-getDefaultPlotParams(plot.type = 1)
pp$leftmargin<-0.15
pp$topmargin <- 15
pp$bottommargin <-15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

ntracks=3

kp<-plotKaryotype(chromosomes = "chr2", zoom=zoom,plot.params = pp)

at <- autotrack(c(1,3),ntracks*2,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(1,ntracks*2,r0=0,r1=1,margin=0.2)
kpHiC(kp,hicPobj = hPhic)

at <- autotrack(c(3,4),ntracks*2,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(c(3),ntracks*2,r0=0,r1=1,margin=0.2)
kpHiC(kp,hicPobj = hPhichip_K4)

at <- autotrack(c(5,6),ntracks*2,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(c(5),ntracks*2,r0=0,r1=1,margin=0.2)
kpHiC(kp,hicPobj = hPhichip_macro)


map<-hPhichip_macro$cMap
n=10000
Xm<-10000
Ym<-5


crossCalc<-function(map,n,Xm,Ym){
  numb<-map$edge[n]
  numbY<-map$Ys[n]
  map %>% filter (Ys >= Ys-Ym & Ys <= Ys+Ym) -> map
  map %>% filter((edge >= numb-Xm) & (edge <= numb+Xm)) %>% filter ((Ys >= numbY-Ym) &(Ys <=numbY +Ym)) -> test1
  print(length(test1))
  res<- sum(test1$counts)
  return(res)
}

vecC<-vector()
for( i in 5000:5100){
  vecC[i]<-crossCalc(map=map,n=i,Xm=10000,Ym=5)
}












#crossCalc(test,1,20000,.1)
#param <- bpparam()
#bpworkers(param) <- 2
#vecC<-unlist(mclapply(1:10000,FUN=crossCalc,map=test,Xm=20000,Ym=.1,mc.cores=10))



