##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/09/20                  Last Modified: 08/23/20 ###
###########################################################################
###########################################################################
###                     deep_coverage.sh                        		###
###########################################################################

# This script idenitifies "true sites" for base recalibration using hard filtered SNPs
# from deep coverage sequences


cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/4.1.6.0
module load samtools
module load vcftools

# Align to chicken genome (change accordingly) 

cd /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step1/
bwa mem -t 50 -M -R "@RG\tID:group1\tSM:E8452\tPL:illumina\tLB:lib1\tPU:unit1" \
/scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
/scratch/snyder/m/mathur20/MQU/ch3/deep_cov/reads/025587_E8452-PE_S44_R1_filtered.fastq \
/scratch/snyder/m/mathur20/MQU/ch3/deep_cov/reads/025587_E8452-PE_S44_R2_filtered.fastq \
> deep_mqu.sam
PicardCommandLine ValidateSamFile I=deep_mqu.sam MODE=SUMMARY O=deep_mqu.sam.txt
PicardCommandLine SortSam INPUT=deep_mqu.sam OUTPUT=sorted_deep_mqu.bam SORT_ORDER=coordinate
PicardCommandLine MarkDuplicates INPUT=sorted_deep_mqu.bam OUTPUT=dedup_deep_mqu.bam METRICS_FILE=metrics_sorted_deep_mqu.bam.txt
PicardCommandLine BuildBamIndex INPUT=dedup_deep_mqu.bam

GenomeAnalysisTK -nt 50 -T RealignerTargetCreator \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-I /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step1/dedup_deep_mqu.bam \
-o /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step2/forIndelRealigner.deep_mqu.intervals

GenomeAnalysisTK -T IndelRealigner \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-I /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step1/dedup_deep_mqu.bam \
-targetIntervals /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step2/forIndelRealigner.deep_mqu.intervals \
-o /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step2/realigned_deep_mqu.bam \
&> /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step2/log/indelrealign.deep.mqu.logfile.txt

PicardCommandLine FixMateInformation \
INPUT=/scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step2/realigned_deep_mqu.bam \
OUTPUT=/scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step3/deep.mqu.sorted.dedup.realigned.fixmate.bam \
SO=coordinate \
CREATE_INDEX=true

samtools depth -a /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step3/deep.mqu.sorted.dedup.realigned.fixmate.bam \
| awk '{c++;s+=$3}END{print s/c}' \
> /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step3/stats/deep_mqu_stats.txt

samtools depth -a /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step3/deep.mqu.sorted.dedup.realigned.fixmate.bam \
| awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' \
> /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step3/stats/deep_mqu_1x_breadth.txt
	
samtools flagstat /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step3/deep.mqu.sorted.dedup.realigned.fixmate.bam \
>  /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step3/stats/deep_mqu_mapped.txt

samtools stats /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step3/deep.mqu.sorted.dedup.realigned.fixmate.bam \
> /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/align/step3/stats/deep_mqu_samtools_stats.txt

########################################################################################################################

######### Base quality score recalibration (chicken) #########

### VARINT CALLING #####
### (For high confidence and "true sites") ####

cd /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/variants

GenomeAnalysisTK -nct 10 -T HaplotypeCaller \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-I /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/bqsr/bam/E8452.chick.sorted.dedup.realigned.fixmate.recal2.bam \
--max_alternate_alleles 30 --genotyping_mode DISCOVERY -pcrModel NONE \
-stand_emit_conf 10 -stand_call_conf 20 \
-o E8452.raw.chick.vcf &> log/deep.chick.variants.log.txt

GenomeAnalysisTK -nct 10 -T HaplotypeCaller \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-I /scratch/snyder/m/mathur20/MQU/ch3/deep_cov/bqsr/bam/E8454.chick.sorted.dedup.realigned.fixmate.recal2.bam \
--max_alternate_alleles 30 --genotyping_mode DISCOVERY -pcrModel NONE \
-stand_emit_conf 10 -stand_call_conf 20 \
-o E8454.raw.chick.vcf &> log/deep.chick.variants.log.txt


gatk --java-options "-Xmx60g" CombineGVCFs \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
--variant E8452.raw.chick.vcf \
--variant E8454.raw.chick.vcf \
-O deep_coverage_rawVariants.g.vcf
 
gatk --java-options "-Xmx500g" GenotypeGVCFs \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V deep_coverage_rawVariants.g.vcf \
--allow-old-rms-mapping-quality-annotation-data \
-O deep_coverage_rawVariants.vcf

gatk SelectVariants \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V deep_coverage_rawVariants.vcf \
--select-type-to-include SNP \
-O deep_coverage_rawSNPs.vcf

## Hard filtering ##
gatk VariantFiltration \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V deep_coverage_rawSNPs.vcf \
--filter-name "QD" \
--filter-expression "QD < 2.0" \
--filter-name "FS" \
--filter-expression "FS > 40.0" \
--filter-name "SOR" \
--filter-expression "SOR > 5.0" \
--filter-name "MQ" \
--filter-expression "MQ < 20.0" \
--filter-name "MQRankSum" \
--filter-expression "-3.0 > MQRankSum || MQRankSum > 3.0" \
--filter-name "ReadPosRankSum" \
--filter-expression "-3.0 > ReadPosRankSum || ReadPosRankSum > 3.0" \
-O deep_coverage_HardFilterSNPs.vcf 

gatk SelectVariants \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V deep_coverage_HardFilterSNPs.vcf \
--select-type-to-include SNP \
-select 'vc.isNotFiltered()' \
-O deep_coverage_HardFilterSNPs.final.vcf



