% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cMapPainter.R
\name{cMapPainter}
\alias{cMapPainter}
\title{Genomic Map Painter}
\usage{
cMapPainter(cmap,high_fact=1e6,pCol="red",alpha=0.05,enhance=NA,colBin=500,ramp=NA,use_ramp=FALSE,log=FALSE,zoom=NULL,...)
}
\arguments{
\item{cmap}{GenomicRanges object created using gInteractions function}

\item{pCol}{charachter - color used if  use_ramp == FALSE ( default "red" )}

\item{enhance}{numeric - transform factor for counts on cMap ( default = NA )}

\item{colBin}{numeric - used to calculate the number of the palette when use_ramp == TRUE ( default  500 )}

\item{ramp}{colorRampPalette object - palette used when use_ramp == TRUE, when is set to NA ramp= ( default NA )}

\item{use_ramp}{logic - select to switch from flat color to colorRampPalette ( default FALSE )}

\item{log}{logic - transform the counts of the contact object in logaritmic ( default NA )}

\item{zoom}{character - genomic coordinate in wich the GRanges will be contruct, if not privide it will use all the contact of GI.
The format need to be "chromosome:start-end" ( default = NA )}
}
\description{
add color and graphic parameters to a GenomicRangers object created with gInteractions function
}
\examples{
 ...

}
\references{
...
}
\seealso{
...
}
