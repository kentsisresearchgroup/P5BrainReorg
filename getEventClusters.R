'Get Event Clusters.

Usage:
  getEventClusters.R <SampleManifest> <VcfFolder> [<ReadCounts>]

If ReadCounts not specific then it defaults to 0

' -> doc

library(docopt)
arguments <- docopt(doc,version='Get Event Clusters 0.9')
if(is.null(arguments$ReadCounts)) {
    readCounts=0
} else {
    readCounts=as.numeric(arguments$ReadCounts)
}

