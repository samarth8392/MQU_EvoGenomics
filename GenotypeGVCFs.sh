##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 06/06/20                  Last Modified: 09/07/20 ###
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

# All 100 #

#GenomeAnalysisTK -T GenotypeGVCFs \
#-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6532/E6532.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6533/E6533.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6535/E6535.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6536/E6536.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6550/E6550.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6567/E6567.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6609/E6609.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6628/E6628.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6845/E6845.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6846/E6846.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6877/E6877.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6962/E6962.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6965/E6965.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7018/E7018.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7019/E7019.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7031/E7031.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7125/E7125.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7145/E7145.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7146/E7146.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7208/E7208.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7209/E7209.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7220/E7220.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7221/E7221.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7563/E7563.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7653/E7653.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7728/E7728.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7746/E7746.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7747/E7747.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7751/E7751.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7752/E7752.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7760/E7760.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7919/E7919.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7927/E7927.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7931/E7931.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7932/E7932.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7934/E7934.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7946/E7946.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7969/E7969.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7970/E7970.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7971/E7971.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7972/E7972.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7986/E7986.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8011/E8011.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8013/E8013.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8015/E8015.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8016/E8016.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8017/E8017.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8018/E8018.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8024/E8024.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8025/E8025.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8028/E8028.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8030/E8030.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8031/E8031.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8032/E8032.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8034/E8034.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8126/E8126.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8142/E8142.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8143/E8143.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8144/E8144.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8148/E8148.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8452/E8452.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8453/E8453.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8454/E8454.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8740/E8740.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8946/E8946.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8947/E8947.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8948/E8948.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8949/E8949.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8954/E8954.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9030/E9030.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9031/E9031.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9032/E9032.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9033/E9033.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9034/E9034.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9035/E9035.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9036/E9036.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9037/E9037.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9038/E9038.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9039/E9039.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9040/E9040.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9041/E9041.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9042/E9042.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9043/E9043.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9044/E9044.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9045/E9045.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9046/E9046.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9047/E9047.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9048/E9048.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9049/E9049.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9050/E9050.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9051/E9051.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9052/E9052.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9053/E9053.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9054/E9054.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9055/E9055.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9056/E9056.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9057/E9057.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9058/E9058.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9059/E9059.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9067/E9067.raw.variants.g.vcf \
#-o /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcfs/mqu.all100.raw.g.vcf


#gatk --java-options "-Xmx120g" GenotypeGVCFs \
#-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcfs/mqu.all100.raw.g.vcf \
#-stand-call-conf 0.0 --allow-old-rms-mapping-quality-annotation-data \
#-O /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcfs/mqu.all100.call0.raw.vcf

#GenomeAnalysisTK -T GenotypeGVCFs \
#-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6532/E6532.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6533/E6533.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6535/E6535.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6536/E6536.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6550/E6550.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6567/E6567.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6609/E6609.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6628/E6628.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6845/E6845.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6846/E6846.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6877/E6877.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6962/E6962.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E6965/E6965.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7018/E7018.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7019/E7019.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7031/E7031.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7125/E7125.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7145/E7145.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7146/E7146.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7208/E7208.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7209/E7209.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7220/E7220.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7221/E7221.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7563/E7563.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7653/E7653.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7728/E7728.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7746/E7746.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7747/E7747.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7751/E7751.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7752/E7752.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7760/E7760.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7919/E7919.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7927/E7927.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7931/E7931.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7932/E7932.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7934/E7934.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7946/E7946.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7969/E7969.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7970/E7970.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7971/E7971.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7972/E7972.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E7986/E7986.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8011/E8011.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8013/E8013.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8015/E8015.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8016/E8016.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8017/E8017.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8018/E8018.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8024/E8024.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8025/E8025.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8028/E8028.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8030/E8030.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8031/E8031.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8032/E8032.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8034/E8034.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8126/E8126.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8142/E8142.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8143/E8143.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8144/E8144.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8148/E8148.raw.variants.g.vcf \
#-o /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcfs/mqu.az.raw.vcf
#
# WTX #

#GenomeAnalysisTK -T GenotypeGVCFs \
#-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8740/E8740.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8954/E8954.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9030/E9030.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9031/E9031.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9032/E9032.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9033/E9033.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9034/E9034.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9035/E9035.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9036/E9036.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9037/E9037.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9038/E9038.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9039/E9039.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9040/E9040.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9041/E9041.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9042/E9042.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9043/E9043.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9044/E9044.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9045/E9045.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9046/E9046.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9047/E9047.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9048/E9048.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9049/E9049.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9050/E9050.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9051/E9051.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9052/E9052.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9053/E9053.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9054/E9054.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9055/E9055.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9056/E9056.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9057/E9057.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9058/E9058.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9059/E9059.raw.variants.g.vcf \
#-o /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcfs/mqu.tx.raw.vcf

# ETX #
#GenomeAnalysisTK -T GenotypeGVCFs \
#-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E9067/E9067.raw.variants.g.vcf \
#-o /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcfs/mqu.E9067.raw.vcf

#MX
#GenomeAnalysisTK -T GenotypeGVCFs \
#-R /scratch/snyder/m/mathur20/MQU/ch3/align/chick/ref/chicken_genome.fa \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8946/E8946.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8947/E8947.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8948/E8948.raw.variants.g.vcf \
#-V /scratch/snyder/m/mathur20/MQU/ch3/variants/chick/E8949/E8949.raw.variants.g.vcf \
#-o /scratch/snyder/m/mathur20/MQU/ch3/variants/final_vcfs/mqu.mx.raw.vcf
