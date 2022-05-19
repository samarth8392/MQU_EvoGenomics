# MQU_EvoGenomics
## Scripts for the Montezuma Quail (MQU) Evolutionary genomics project

Samarth Mathur, John Tomeček, Luis Tarango-Arámbula, Robert M. Perez, J. Andrew DeWoody (2022). 
An evolutionary perspective on genetic load in small isolated populations. *_Molecular Ecology_* (revision resubmitted).

Initial Preprint: 10.22541/au.162495929.94655412/v1

### The summarized repository structure is presented below ###
##### (See each repository for substrutures and details)

```
├── BashScript
│   ├── Diversity&Structure
│   │   ├── heterozygosityKb.sh
│   │   ├── norelMAF.sh
│   │   └── tassel.sh
│   ├── MutationAge&Load
│   │   ├── genetic_load.sh
│   │   ├── mutation_age.sh
│   │   └── snp_annotation.sh
│   ├── ROH,Ancestry,&Demography
│   │   ├── HaplotypePhasing.sh
│   │   ├── ROH_analysis.sh
│   │   ├── finestructure.sh
│   │   ├── gadma_2pop.sh
│   │   ├── gadma_3pop.sh
│   │   └── gadma_boostrap.sh
│   ├── ReadPreproccesing
│   │   ├── adapter_removal.sh
│   │   ├── alignment.sh
│   │   ├── base_recalibration.sh
│   │   └── depth&breadth.sh
│   ├── SLiM_Simulations
│   │   ├── 1k.expand.slim
│   │   ├── 1k.postbottle.slim
│   │   ├── 1k.prebottle.slim
│   │   ├── getSLiMoutput.sh
│   │   └── slimFst.sh
│   └── VariantCalling
│       ├── GenotypeGVCFs.sh
│       ├── HaplotypeCaller.sh
│       └── vcf_filtering.sh
├── RScripts
│   ├── Load_PotReal.R
│   ├── Load_segDrift.R
│   ├── ROH_tryruns.R
│   ├── ROHs.R
│   ├── Rxy.R
│   ├── admixture.R
│   ├── analyzeSLiMOutput.R
│   ├── diversity.R
│   ├── getSLiMoutput.R
│   ├── heterozygosity_norel.R
│   ├── mut_age.R
│   ├── mutation_age.R
│   ├── population_tree.R
│   ├── recomb.R
    └── simple-sensitivity-analysis-with-r
```