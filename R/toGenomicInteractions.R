toGenomicInterations<-function(bedpe, counts=NA,style="UCSC"){

  if (is.na(counts)){
    bedpe %>% mutate( counts = 1) -> bedpe
  } else {
    if (is.character(counts)){
      bedpe %>% mutate( counts = eval(parse(text = counts))) -> bedpe
    } else {
        stop ( "counts need to be a columnname of bedpe or NA" )
      }
   }

  chr<-bedpe[,1]
  anchor1<-toGRanges(data.frame(chr=chr,start=bedpe[,2],end=bedpe[,3]))
  chr<-bedpe[,4]
  anchor2<-toGRanges(data.frame(chr=chr,start=bedpe[,5],end=bedpe[,6]))
  GInt<-GenomicInteractions(anchor1 = anchor1,anchor2 = anchor2,counts = bedpe$counts)
  seqlevelsStyle(GInt)<-style
  return(GInt)

}
