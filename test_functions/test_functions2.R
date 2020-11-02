zoom<-'chr10:53,000,000-55,000,000'
zoom<-createZooms(zoom)

hic_TAD_10k_file<-"/mnt/data1/HiCUP/Macro_Insulation/HIC/ENCODE_hepg2/TADs/ENCODE_hic_hepg2.chr20.out/10000_blocks.bedpe"
hic_TAD_5K_file<-"/mnt/data1/HiCUP/Macro_Insulation/HIC/ENCODE_hepg2/TADs/ENCODE_hic_hepg2.chr20.out/5000_blocks.bedpe"
GI_hic_TAD_10K<-arrowHead2GI(hic_TAD_5k_file)
hP_TADs_10K<-createHicPainteRObj(Name = "hic_HepG2_loops_chr10",mapType = 'TAD',contact = GI_hic_TAD_10K,
                                 zoom=zoom$GR,pCol="black")

GI_hic_loops<-hiccups2GI(hic_loops_file)
hP_loops<-createHicPainteRObj(Name = "hic_HepG2_loops_chr10",mapType = 'Loop',contact = GI_hic_loops,
                              zoom=zoom1, pCol="darkgreen")

Hic_DF<-DFfromHIC(hicFile = "/home/malli/work/insulation/test_data/chr10/ENCODE_hic_hepg2.chr10.hic",
                  zoom = zoom,
                  seq_style="NCBI",bin=10000,norm="KR")

Hic_GI <- hicDump2GI(Hic_DF,chr=zoom$chr)
Hic_hicPO <- createHicPainteRObj(Hic_GI,mapType = "cMap",use_ramp = TRUE,ramp=NA)



ramp<-ramp<-colorRampPalette(c("cornflowerblue","yellow","orange","red"))(100)
a<-Hic_hicPO$cMap
yreal<-(end(a)-start(a))<=1.5e6
a<-a[yreal]
a$counts<-log1p(a$orig.counts)
zscores<-abs(a$counts-mean(a$counts))/sd(a$counts)
a<-a[zscores<=3]

a$counts<-log1p(a$orig.counts)
a$counts<-round(rescale(a$counts,to=c(1,100)))
a$colors<-ramp[a$counts]
a$colors<-alpha(a$colors,alpha=a$counts/100)
a %>% arrange(counts) ->a


Hic_hicPO1<-Hic_hicPO
Hic_hicPO1$cMap<-a


kp<-plotKaryotype(zoom=zoom$GR,plot.params = pp)
genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = FALSE,
                                    plot.transcripts.structure = FALSE)
genes.data <- addGeneNames(genes.data)
kpAddBaseNumbers(kp, tick.dist = 1e6, minor.tick.dist = 1e5,
                 add.units = TRUE, cex=0.7, tick.len = 3)

at <- autotrack(c(1,2),ntracks*2,r0=0,r1=1,margin=0.1)
kpHiC(kp,hicPobj = Hic_hicPO3,r0 = at$r0, r1=at$r1)

at <- autotrack(c(3,4),ntracks*2,r0=0,r1=1,margin=0.1)
Hic_hicPO2<-Hic_hicPO1
Hic_hicPO2$param$cexP<-2.5
kpRect(kp, chr=zoom$chr, x0=as.numeric(zoom$start), x1=as.numeric(zoom$end), y0=at$r0, y1=at$r1,col="white",border=NA)
kpHiC(kp,hicPobj = Hic_hicPO,r0 = at$r0, r1=at$r1,ymax="visible.region")
at <- autotrack(c(5,6),ntracks*2,r0=0,r1=1,margin=0.1)
Hic_hicPO3<-Hic_hicPO1
Hic_hicPO3$param$cexP<-3
kpRect(kp, chr=zoom$chr, x0=as.numeric(zoom$start), x1=as.numeric(zoom$end), y0=at$r0, y1=at$r1,col="white",border=NA)
kpHiC(kp,hicPobj = Hic_hicPO3,r0 = at$r0, r1=at$r1)
  kpLoops(kp,hicPobj = hP_loops,loopType = "Points")
at <- autotrack(c(1),ntracks*2,r0=0,r1=1,margin=0.1)
kpPlotGenes(kp, data=genes.data)
