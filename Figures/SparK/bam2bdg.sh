#!/bin/bash

BAM=$1

# SCALE=$(samtools idxstats $BAM | egrep "^[1-2]" | awk '{print s+=$3/1000000}' | tail -1 | awk '{print 1/$1}')

# echo SCALE=$SCALE

#OUT=$(basename $BAM | sed 's/.bam/.bdg/')
# bedtools genomecov -scale $SCALE -ibam $BAM -bg \
#     | bedtools intersect -a - -b mouse.GRCm38.bed -wa \
#     | sort -S 1g -k1,1V -k2,2n \
#     > $OUT

OUT=$(basename $BAM | sed 's/.bam/.UnScale.bdg/')
bedtools genomecov -ibam $BAM -bg \
    | bedtools intersect -a - -b mouse.GRCm38.bed -wa \
    | sort -S 1g -k1,1V -k2,2n \
    > $OUT

bgzip $OUT
tabix -p bed ${OUT}.gz


