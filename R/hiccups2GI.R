hiccups2GI<-function(bedpe,counts="observed"){

  if(!is.data.frame(bedpe)){
    print("uploading dump file......")
    bedpe<-read.delim(bedpe,header=TRUE)
    bedpe<-bedpe[-1,]
  }
  print(nrow(bedpe))
  print(counts)
  GInt<-toGenomicInterations(bedpe,counts=counts)

  return(GInt)

}
