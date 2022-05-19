```
.
├── Diversity
│   ├── heterozygosityKb.sh
│   ├── norelMAF.sh
│   └── tassel.sh
├── MutationAge&Load
│   ├── genetic_load.sh
│   ├── mutation_age.sh
│   └── snp_annotation.sh
├── ROH,Ancestry,&Demography
│   ├── HaplotypePhasing.sh
│   ├── ROH_analysis.sh
│   ├── finestructure.sh
│   ├── gadma_2pop.sh
│   ├── gadma_3pop.sh
│   └── gadma_boostrap.sh
├── ReadPreproccesing
│   ├── adapter_removal.sh
│   ├── alignment.sh
│   ├── base_recalibration.sh
│   └── depth&breadth.sh
├── SLiM_Simulations
│   ├── 1k.expand.slim
│   ├── 1k.postbottle.slim
│   ├── 1k.prebottle.slim
│   ├── getSLiMoutput.sh
│   └── slimFst.sh
└── VariantCalling
    ├── GenotypeGVCFs.sh
    ├── HaplotypeCaller.sh
    └── vcf_filtering.sh
```

#### Scripts in the order of pipeline ####

##### Read Preproccesing #####

1. _Adapter removal_ : Trim raw sequences to remove adapters and low-quality bases (Phred < 20) using Trimmomatic
2. _Alignment.sh_ : Mapping and read processing following GATK “Best Practice Workflow”
3. _Base_recalibration.sh_ : Recaliberate mapped sequences using known variants from high quality sequnce reads
4. _Depth&Breadth.sh_ : Calculate map stats using SAMtools

##### VariantCalling #####

1. _HaplotypeCaller.sh_ : Call genome-wide variants in each sample (GVCF) using GATK v4
2. _GenotypeGVCFs.sh_: Genotype calling from combined GVCF information of all samples and SNP filtering
3. _vcf_filtering.sh_ : SNP filtering using VCFtools

##### Diversity #####

1. _HeterozygosityKb.sh_ : Divide genome-wide SNPs by chromosomes and calucate heterozygosity in 1 Mb windows
2. _norelMAF.sh_ : Comparing diversity estimates between Arizona and Texas populations after removing close relatives