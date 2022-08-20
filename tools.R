suppressPackageStartupMessages({
    require(bictools)
    require(stringr)
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

    #
    # Get VCF fields that are simple numbers by looking at header
    #
    numberFields=grep("Type=(Integer|Float)",vv[[1]]$header,value=T) %>%
        grep("Number=1",.,value=T) %>%
        str_extract("ID=[^,]+,") %>%
        gsub("^ID=","",.) %>%
        gsub(",$","",.)
    # cat("NumberFields =",numberFields,"\n")

    full_join(vm,vs,by = c("FILE", "VID")) %>%
        mutate(UUID=paste0(FILE,":",ID,":",VID)) %>%
        mutate_at(numberFields,as.numeric) %>%
        arrange(factor(CHROM,levels=c(1:19,"X","Y")),POS,END) %>%
        left_join(manifest,by = "SAMPLE") %>%
        select(CHROM,POS,END,UUID,SAMPLE,GROUP,MOUSE,BRAIN_AREA,GENOTYPE,everything())

}

read_mutect2_events<-function(sampleManifest,vcfFiles) {

    manifest=readxl::read_xlsx(sampleManifest) %>%
        select(SAMPLE,GROUP,MOUSE,BRAIN_AREA,GENOTYPE) %>%
        mutate(SID=gsub("_IGO_.*","",SAMPLE))

    if(verbose) cat("\nReading in vcf files ...")

    vv=map(vcfFiles,read_vcf)
    names(vv)=map(strsplit(gsub(".*(adult|embrionary).","",vcfFiles),"/"),1) %>% unlist

    if(verbose) cat(" done\n")

    vs=map(vv,"vs") %>%
        map(~filter(.,SAMPLE=="TUMOR")) %>%
        bind_rows(.id="SID") %>%
        separate(AD,c("RD","AD"),sep=",") %>%
        mutate(RD=as.numeric(RD),AD=as.numeric(AD),D0=RD+AD) %>%
        mutate(AF0=AD/D0)

    vm=map(vv,"vm") %>%
        map(~mutate(.,CHROM=as.character(CHROM))) %>%
        bind_rows(.id="SID")

    #mutate_at(numberFields,as.numeric) %>%

    full_join(vm,vs,by = c("SID", "VID")) %>%
        select(-SAMPLE) %>%
        mutate(UUID=paste0(SID,":",VID)) %>%
        mutate(ETAG=paste0(CHROM,":",POS,":",REF,":",ALT)) %>%
        arrange(factor(CHROM,levels=c(1:19,"X","Y")),POS) %>%
        left_join(manifest,by = "SID") %>%
        select(CHROM,POS,UUID,SAMPLE,SID,GROUP,MOUSE,BRAIN_AREA,GENOTYPE,everything())


}

