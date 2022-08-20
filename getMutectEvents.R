'Get Mutect2 Clusters.

Usage:
  getMutectEvents.R [options] <sampleManifest> <fileOfVCFFiles>

Options:
  -v --verbose   Printing debugging/logging info
  --ad=<ad>      Set allele depth filter [default: 1]
  --af=<af>      Set allele frequency filter [default: 0]
  --dp=<dp>      Set total depth filter [default: 1]

' -> doc

library(docopt)
argv <- docopt(doc,version='Get Mutect Clusters 1.0')

adFilter=as.numeric(argv$ad)
afFilter=as.numeric(argv$af)
dpFilter=as.numeric(argv$dp)
verbose=argv$verbose

suppressPackageStartupMessages({
    require(dplyr)
    require(readr)
    require(purrr)
    require(digest)
})

source("P5BrainReorg/tools.R")

vcfFiles=scan(argv$fileOfVCFFiles,"")

cacheFile=cc("vvObj",digest(c(argv$sampleManifest,vcfFiles)),".rds")

if(file.exists(cacheFile)) {

    vv=readRDS(cacheFile)

} else {

    vv=read_mutect2_events(argv$sampleManifest,vcfFiles)
    saveRDS(vv,cacheFile,compress=T)

}

pg1=vv %>% mutate(AF0=AD/(AD+RD)) %>% ggplot(aes(AF,AF0,color=AD>0)) + theme_light() + geom_point(alpha=.5) + ggtitle(paste("CHR",19)) + coord_fixed()
pg2=vv %>% mutate(AF0=AD/(AD+RD)) %>% ggplot(aes(AF,AF0,color=FILTER=="PASS")) + theme_light() + geom_point(alpha=.5) + ggtitle(paste("CHR",19)) + coord_fixed()
pg3=vv %>% filter(FILTER=="PASS") %>% mutate(AF1=AD/(AD+RD)) %>% ggplot(aes(AF,AF0,color=AD>0)) + theme_light() + geom_point(alpha=.5) + ggtitle(paste("CHR",19,"PassOnly")) + coord_fixed()

vv %>% mutate(AF0=AD/(AD+RD),DP0=AD+RD)