suppressPackageStartupMessages({
    require(bictools)
})

read_events<-function(sampleManifest,vcfFolder) {

    manifest=readxl::read_xlsx(sampleManifest) %>% select(SAMPLE,GROUP,MOUSE,BRAIN_AREA,GENOTYPE)

    if(verbose) cat("\nReading in vcf files ...")

    vv=fs::dir_ls(vcfFolder,recur=T,regex="\\.vcf") %>% map(read_vcf)

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

