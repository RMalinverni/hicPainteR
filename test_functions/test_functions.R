library("karyoploteR")
library("GenomicInteractions")
library("plyranges")
library("ggplot2")

hic<-read.delim("/mnt/data1/HiCUP/Macro_Insulation/HIC/ENCODE_hepg2/dumps/ENCODE_hic_hepg2.chr2.dump.txt",header=FALSE)
hichip<-read.delim("/mnt/data1/HiCUP/Macro_Insulation/HIChIP/h3k4me3_hichip/hic/hichip_h3k4me3_chr2_dump.txt",header=FALSE)
hichipM<-read.delim("/mnt/data1/HiCUP/Macro_Insulation/HIChIP/macroH2A_hichip/hic/hichip_macro_chr2_dump.txt",header=FALSE)

GIhic<-hicDump2GI(hic,chr = "chr2")
GIhichip<-hicDump2GI(hichip,chr = "chr2")
GIhichipM<-hicDump2GI("/mnt/data1/HiCUP/Macro_Insulation/HIChIP/macroH2A_hichip/hic/hichip_macro_chr2_dump.txt",chr = "chr2")


rm(hic)
rm(hichip)


hPhic<-createHicPainteRObj(Name = "HepG2_chr2",Type = 'HiC',contact = GIhic,
                           zoom='chr2:100,000,000-110,000,000')

hPhichip<-createHicPainteRObj(Name = "HepG2_chr2",Type = 'HiC',contact = GIhichip,
                              zoom='chr2:100,000,000-110,000,000')

hPhichipM<-createHicPainteRObj(Name = "HepG2_chr2",Type = 'HiC',contact = GIhichipM,
                               zoom='chr2:100,000,000-110,000,000')
file1<-"/mnt/data1/HiCUP/Macro_Insulation/HIC/ENCODE_hepg2/TADs/ENCODE_hic_hepg2.chr2.out/5000_blocks.bedpe"
file2<-"/mnt/data1/HiCUP/Macro_Insulation/HIC/ENCODE_hepg2/loops/ENCODE_hic_hepg2.chr2.hiccups_loops/enriched_pixels_5000.bedpe"
TAD<-arrowHead2GI(bedpe=file1)
Loops<-hiccups2GI(bedpe=file2)



Loops<-read.delim("/mnt/data1/HiCUP/Macro_Insulation/HIC/ENCODE_hepg2/loops/ENCODE_hic_hepg2.chr2.hiccups_loops/enriched_pixels_5000.bedpe")
Loops<-Loops[-1,]

gInteractionsTocMap(GI=Loops)
bedpe<-TAD

toGenomicInterations(TAD)

gInteractionsTocMap(TAD)
gInteractionsTocMap(Loops)



toGenomicInterations<-function(bedpe,counts,style="UCSC"){

  if (!is.na(counts)){
    if (is.character(counts)){
      bedpe %>% mutate( counts = counts) -> bedpe
    }else{
      bedpe$counts <- counts
    }
  }

  chr<-bedpe[,1]
  anchor1<-toGRanges(data.frame(chr=chr,start=bedpe[,2],end=bedpe[,3]))
  chr<-bedpe[,4]
  anchor2<-toGRanges(data.frame(chr=chr,start=bedpe[,5],end=bedpe[,6]))
  GInt<-GenomicInteractions(anchor1 = anchor1,anchor2 = anchor2,counts = counts)
  seqlevelsStyle(GInt)<-style
  return(GInt)

}

Loops<-toGenomicInterations(Loops,c.counts="observed")

gInteractionsTocMap(Loops)

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
hPhichipM<-repaintHicPainteRobj(hPhichipM,alpha = 0.05,cexP = 1,log=FALSE,use_ramp = FALSE)
kpHiC(kp,hicPobj = hPhichipM)


#fakeObj<-hPObj
#fakeObj$cMap$orig.counts[fakeObj$cMap$orig.counts > 580]<-0.001

at <- autotrack(c(3,4),ntracks*2,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(c(3),ntracks*2,r0=0,r1=1,margin=0.2)
hPhichip<-repaintHicPainteRobj(hPhichip,alpha = 0.0005,cexP = 2,log=TRUE,use_ramp = FALSE)
kpHiC(kp,hicPobj = hPhichip)


at <- autotrack(c(5,6),ntracks*2,r0=0,r1=1,margin=0)
kpRect(kp,chr=chr,x0= SS, x1=EE,y0=0, y1=1,col="white",border=NA, r0=at$r0,r1=at$r1)
at <- autotrack(c(5),ntracks*2,r0=0,r1=1,margin=0.2)
hPhic<-repaintHicPainteRobj(hPhic,alpha = 0.0005,cexP = 2,log=TRUE,use_ramp = TRUE,enhance = 1)
kpHiC(kp,hicPobj = hPhic)



test<-hPObj$cMap

test %>% filter (ynorm)

crossCalc<-function(map,n,Xm,Ym){
  numb<-map$edge[n]
  numbY<-test$Ynorm[n]
  limitYup<-2e6/(max(mapt$edge)-min(map$edge))
  limitYdown<-1e5/(max(map$edge)-min(map$edge))
  map %>% filter (Ynorm >= limitYdown & Ynorm <= limitYup)
  map %>% filter((edge >= numb-Xm) & (edge <= numb+Xm)) %>% filter ((Ynorm >= numbY-Ym) &(Ynorm <=numbY +Ym)) -> test1


  #res<- mean(test1$orig.counts)
  return(res)
}










#crossCalc(test,1,20000,.1)
#param <- bpparam()
#bpworkers(param) <- 2
#vecC<-unlist(mclapply(1:10000,FUN=crossCalc,map=test,Xm=20000,Ym=.1,mc.cores=10))



