# MQU_EvoGenomics
## Scripts for the Montezuma Quail (MQU) Evolutionary genomics project

Samarth Mathur, John Tomeček, Luis Tarango-Arámbula, et al. 
An evolutionary perspective on contemporary genetic load in small isolated populations. *_Molecular Ecology_* (revision resubmitted).

Initial Preprint: 10.22541/au.162495929.94655412/v1

### The summarized repository structure is presented below ###
##### (See each repository for substrutures)

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

```