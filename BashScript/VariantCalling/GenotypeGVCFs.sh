##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     genotypeGVCFs.sh                   				###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/4.1.6.0
module load samtools

# Best 66 #

gatk --java-options "-Xmx220g" GenotypeGVCFs \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6536/E6536.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6609/E6609.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6628/E6628.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6846/E6846.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6877/E6877.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7031/E7031.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7125/E7125.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7146/E7146.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7208/E7208.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7220/E7220.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7563/E7563.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7746/E7746.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7747/E7747.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7751/E7751.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7752/E7752.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7927/E7927.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7932/E7932.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7934/E7934.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7946/E7946.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7969/E7969.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8013/E8013.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8017/E8017.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8024/E8024.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8025/E8025.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8030/E8030.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8031/E8031.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8032/E8032.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8142/E8142.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8946/E8946.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8947/E8947.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8948/E8948.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8949/E8949.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8954/E8954.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9030/E9030.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9031/E9031.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9032/E9032.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9033/E9033.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9034/E9034.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9035/E9035.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9036/E9036.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9037/E9037.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9038/E9038.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9039/E9039.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9040/E9040.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9041/E9041.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9042/E9042.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9043/E9043.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9044/E9044.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9045/E9045.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9046/E9046.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9047/E9047.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9048/E9048.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9049/E9049.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9050/E9050.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9051/E9051.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9052/E9052.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9053/E9053.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9054/E9054.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9055/E9055.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9056/E9056.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9057/E9057.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9058/E9058.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9059/E9059.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9067/E9067.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9569/E9569.raw.variants.g.vcf \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9570/E9570.raw.variants.g.vcf \
-O /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcfs/mqu.best66.raw.g.vcf


# GATK 
gatk --java-options "-Xmx220g" SelectVariants \
--select-type-to-include SNP \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcfs/mqu.best66.raw.g.vcf \
-O /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcfs/mqu.best66.raw.GATK.SNPs.vcf

gatk --java-options "-Xmx100g" VariantFiltration \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcf/mqu.best66.raw.GATK.SNPs.vcf \
-O /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcf/mqu.best66.filterTag.GATK.SNPs.vcf \
--filter-name "QD" \
--filter-expression "QD < 2.0" \
--filter-name "FS" \
--filter-expression "FS > 40.0" \
--filter-name "SOR" \
--filter-expression "SOR > 5.0" \
--filter-name "MQ" \
--filter-expression "MQ< 20.0" \
--filter-name "MQRankSum" \
--filter-expression " MQRankSum < -3.0 || MQRankSum > 3.0" \
--filter-name "ReadPosRankSum" \
--filter-expression "ReadPosRankSum < -8.0" \
--filter-name "AN" \
--filter-expression "AN < 130"

gatk --java-options "-Xmx220g" SelectVariants \
-select 'vc.isNotFiltered()' \
-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
-V /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcf/mqu.best66.filterTag.GATK.SNPs.vcf \
-O /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcf/mqu.best66.final.GATK.SNPs.vcf