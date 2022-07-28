'Get Event Clusters.

Usage:
  getEventClusters.R [options] <sampleManifest> <vcfFolder> [<readCounts>]

Options:
  -v --verbose   Printing debugging/logging info

If ReadCounts not specified then it defaults to 0

' -> doc

library(docopt)
argv <- docopt(doc,version='Get Event Clusters 1.0')
if(is.null(argv$readCounts)) {
    readCounts=0
} else {
    readCounts=as.numeric(argv$readCounts)
}

verbose=argv$verbose

suppressPackageStartupMessages({
    require(dplyr)
    require(readr)
    require(purrr)
    require(fs)
    require(readxl)
    require(bictools)
    require(tidygenomics)
    require(igraph)
})

read_events<-function(sampleManifest,vcfFolder) {

    manifest=read_xlsx(sampleManifest) %>% select(SAMPLE,GROUP,MOUSE,BRAIN_AREA,GENOTYPE)

    if(verbose) cat("\nReading in vcf files ...")

    vv=dir_ls(vcfFolder,recur=T,regex="\\.vcf") %>% map(read_vcf)

    if(verbose) cat(" done\n")

    vs=map(vv,"vs") %>%
    map(~filter(.,SAMPLE %in% .$SAMPLE[1])) %>%
    bind_rows(.id="FILE") %>%
    mutate(FILE=basename(FILE)%>%gsub(".pre.ca.*","",.))

    vm=map(vv,"vm") %>%
        bind_rows(.id="FILE") %>%
        mutate(FILE=basename(FILE)%>%gsub(".pre.ca.*","",.))

    full_join(vm,vs,by = c("FILE", "VID")) %>%
        mutate(UUID=paste0(FILE,":",ID,":",VID)) %>%
        mutate(END=as.numeric(END)) %>%
        mutate(RC=as.numeric(RC)) %>%
        arrange(factor(CHROM,levels=c(1:19,"X","Y")),POS,END) %>%
        left_join(manifest,by = "SAMPLE") %>%
        select(CHROM,POS,END,UUID,SAMPLE,GROUP,MOUSE,BRAIN_AREA,GENOTYPE,everything())

}

filter_events<-function(events,readCounts) {
    events %>%
        filter(FILTER=="PASS" & RC>=readCounts) %>%
        mutate(SIZE=END-POS) %>%
        select(CHROM,POS,END,UUID,SIZE,SAMPLE,GROUP,MOUSE,BRAIN_AREA,GENOTYPE,RC)
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

get_event_clusters<-function(sampleManifest,vcfFolder,readCounts) {

    events = read_events(sampleManifest,vcfFolder)

    events_filtered = filter_events(events,readCounts)

    events_intersect=genome_intersect(events_filtered,events_filtered,by=c("CHROM","POS","END")) %>% filter(POS!=END)

    edges=events_intersect %>%
        filter(UUID.x<UUID.y) %>%
        mutate(SIZE_AVG=(SIZE.x+SIZE.y)/2,PCT_OVER=(END-POS)/SIZE_AVG) %>%
        filter((SIZE.x<1e5 & SIZE.y<1e5) | PCT_OVER>0.001) %>%
        select(matches("UUID"))

    eventGraph=graph_from_data_frame(edges,directed=F)

    event_cliques=max_cliques(eventGraph)

    clusters=clique_to_event_table(event_cliques,events_filtered)

    clusters

}

clusters=get_event_clusters(argv$sampleManifest,argv$vcfFolder,readCounts)

outFile=cc("clusters","","Small,Large_0.001_Overlap","","RC",sprintf("%04d",readCounts),".csv")

write_csv(clusters,outFile)

