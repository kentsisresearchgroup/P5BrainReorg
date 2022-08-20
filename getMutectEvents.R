suppressPackageStartupMessages({
    require(dplyr)
    require(readr)
    require(purrr)
})

source("P5BrainReorg/tools.R")

verbose=T
vcfFiles=scan("mutect2VCFs","")

vv=read_mutect2_events("raw/manifest_B-101-295__DellyVCFs_V2.xlsx")

vv %>% select(SID,ETAG) %>% group_by(ETAG) %>% summarize(SSID=paste0(SID,collapse=";"))