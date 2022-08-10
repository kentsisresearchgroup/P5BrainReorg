'Get Edge Clusters.

Usage:
  getEdgeClusters.R [options] <sampleManifest> <vcfFolder>

Options:
  -v --verbose   Printing debugging/logging info
  -i=<ins> --ins=<ins> Insert Size [default: 250]
  --rv=<rv>      Set RV filter [default: 0]
  --vaf=<vaf>    Set VAF {RV/(RR+RV)} filter [default: 0]

' -> doc

library(docopt)
argv <- docopt(doc,version='Get Event Clusters 1.0')

rcFilter=as.numeric(argv$rc)
rvFilter=as.numeric(argv$rv)
vafFilter=as.numeric(argv$vaf)
insertSize=as.numeric(argv$ins)
verbose=argv$verbose

suppressPackageStartupMessages({
    require(dplyr)
    require(readr)
    require(purrr)
    require(tidygenomics)
    require(igraph)
})

# tools.R
#   read_events<-function(sampleManifest,vcfFolder)
source("P5BrainReorg/tools.R")

filter_events<-function(events) {
    events %>%
        mutate(SIZE=END-POS) %>%
        filter(FILTER=="PASS" & RC/SIZE>=rcFilter & RV>=rvFilter & RV/(RV+RR)>vafFilter) %>%
        select(CHROM,POS,END,UUID,SIZE,SAMPLE,GROUP,MOUSE,BRAIN_AREA,GENOTYPE,RC,RV)
}

clique_to_event_table<-function(cliques,events_filtered) {

    clusters=list()

    for(cii in seq(cliques)) {

        if(verbose) cat("Cluster, maxC, clusterSize=",c(cii,len(cliques),len(cliques[[cii]])),"\n")

        nii=names(cliques[[cii]])
        clusters[[cii]]=events_filtered %>%
            filter(UUID %in% nii) %>%
            mutate(CLUSTER=cc("Cluster",cii),CLUSTER_LEN=len(nii)) %>%
            select(CLUSTER,CLUSTER_LEN,everything())

    }

    bind_rows(clusters) %>% filter(CLUSTER_LEN>2)

}

get_event_clusters<-function(sampleManifest,vcfFolder) {

    events = read_events(sampleManifest,vcfFolder)

    events_filtered = filter_events(events)

    if(nrow(events_filtered)==0) {
        return(NULL)
    }

    events_intersect=genome_intersect(events_filtered,events_filtered,by=c("CHROM","POS","END")) %>% filter(POS!=END)

    edges=events_intersect %>%
        filter(UUID.x<UUID.y) %>%
        mutate(SIZE_AVG=(SIZE.x+SIZE.y)/2,PCT_OVER=(END-POS)/SIZE_AVG) %>%
        filter((SIZE.x<1e5 & SIZE.y<1e5) | PCT_OVER>0.001) %>%
        select(matches("UUID"))

    if(nrow(edges)==0) {
        return(NULL)
    }

    eventGraph=graph_from_data_frame(edges,directed=F)

    event_cliques=max_cliques(eventGraph)

    clusters=clique_to_event_table(event_cliques,events_filtered)

    clusters

}

clusters=get_event_clusters(argv$sampleManifest,argv$vcfFolder)

if(!is.null(clusters) && nrow(clusters)>0) {

    outFile=cc("clusters","EDGE","Small,Large_0.001_Overlap","",
                    "RCFilter",sprintf("%05d",rcFilter),
                    "RVFilter",sprintf("%05d",rvFilter),
                    "VAFFilter",sprintf("%07.03f",round(100*vafFilter,3)),
                    ".csv")

    write_csv(clusters,outFile)

}

