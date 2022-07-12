# P5BrainReorg

## Recurrence analysis of somatic neuronal DNA rearrangements

### BIC Project B-101-295

### TASKS (2022-07-12)

1) SNV Analysis

    - Take VCF files from GATK #directory with specific file names. Analyze recurrence of somatic SNVs from GATK among 3 individuals and 3 brain regions. Generate script for producing bed files for counting statistics, where the script will include user-defined parameters for DP and VAF

2) Indel analysis

    - Same as SNV analysis but GCF files from Pindel #directory
    - Script will include user-defined parameters for Pindel supporting reads

3) Two types of SV analyses from VCF files from Delly2 #directory

    - Breakpoint specific. Overlap is insert size #Luz
    - Region specific

Script is going to iterate over the number of supporting reads (RC) from 3 to max, where max yields zero somatic variants in knockout versus wildtype brains

