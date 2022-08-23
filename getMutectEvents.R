'Get Mutect2 Clusters.

Usage:
  getMutectEvents.R [options] <sampleManifest> <fileOfVCFFiles>

Options:
  -v --verbose   Printing debugging/logging info
  --ad=<ad>      Set allele depth filter [default: 1]
  --af=<af>      Set allele frequency filter [default: 0]
  --dp=<dp>      Set total depth filter [default: 1]
  --pass         Set to filter only PASS events

' -> doc

library(docopt)
argv <- docopt(doc,version='Get Mutect Clusters 1.0')

adFilter=as.numeric(argv$ad)
afFilter=as.numeric(argv$af)
dpFilter=as.numeric(argv$dp)
verbose=argv$verbose
FILTER_FLAG=ifelse(argv$pass,"PassOnly","AllEvents")

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

# pg1=vv %>% mutate(AF0=AD/(AD+RD)) %>% ggplot(aes(AF,AF0,color=AD>0)) + theme_light() + geom_point(alpha=.5) + ggtitle(paste("CHR",19)) + coord_fixed()
# pg2=vv %>% mutate(AF0=AD/(AD+RD)) %>% ggplot(aes(AF,AF0,color=FILTER=="PASS")) + theme_light() + geom_point(alpha=.5) + ggtitle(paste("CHR",19)) + coord_fixed()
# pg3=vv %>% filter(FILTER=="PASS") %>% mutate(AF1=AD/(AD+RD)) %>% ggplot(aes(AF,AF0,color=AD>0)) + theme_light() + geom_point(alpha=.5) + ggtitle(paste("CHR",19,"PassOnly")) + coord_fixed()

vv.o=vv

cat("N.events =",nrow(vv),"\n")

if(FILTER_FLAG=="PassOnly")
    vv = vv %>% filter(FILTER=="PASS")

cat("N.events =",nrow(vv),"\n")


#vv=vv.o
vv=vv %>%
    mutate(RD2=REF_F1R2+REF_F2R1,AD2=ALT_F1R2+ALT_F2R1) %>%
    filter(2*abs(RD2-RD)/(RD2+RD)<.1 & 2*abs(AD2-AD)/(AD2+AD)<.1) %>%
    mutate(AF0=AD/(AD+RD),DP0=AD+RD) %>%
    filter(AD>=adFilter) %>%
    mutate(DPr=ifelse(AF<.4,round(RD/(1-AF),0),NA),DPa=round(AD/AF,0)) %>%
    mutate(mDP=pmax(pmax(DP0,DPr),DPa)) %>%
    mutate(AF0=AD/mDP) %>%
    filter(mDP>=dpFilter & AF0 >= afFilter)

cat("N.events =",nrow(vv),"\n")

    # select(SID,VID,CHROM,POS,REF,ALT,FILTER,RD,RD2,AD,AD2,AF,DP0,DPr,DPa)

    #rowwise() %>% mutate(PCT.DIFF=(max(c(DP0,DPr,DPa),na.rm=T)-min(c(DP0,DPr,DPa),na.rm=T))/(mean(c(DP0,DPr,DPa),na.rm=T))) %>% arrange(desc(PCT.DIFF))

tbl=vv %>%
    arrange(ETAG) %>%
    group_by(ETAG) %>%
    mutate(CLUSTER_LEN=n()) %>%
    mutate(CLUSTER=cc("CLUSTER",cur_group_id())) %>%
    filter(CLUSTER_LEN>2) %>%
    ungroup %>%
    arrange(CLUSTER_LEN) %>%
    select(
        CLUSTER,CLUSTER_LEN,
        CHROM,POS,REF,ALT,
        SID,GROUP,MOUSE,BRAIN_AREA,GENOTYPE,
        mDP,AD,AF0
        )

if(nrow(tbl)>0) {
    outFile=cc("cluster","Mutect2","V2",FILTER_FLAG,"",
                "DPFilter",sprintf("%04d",dpFilter),
                "ADFilter",sprintf("%04d",adFilter),
                "VAFFilter",sprintf("%07.03f",round(100*afFilter,3)),
                ".csv")

    write_csv(tbl,outFile)
}
