##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     tassel.sh                                       ###
###########################################################################


cd $SLURM_SUBMIT_DIR
module load bioinfo
module load TASSEL/5.2.44

cd /scratch/bell/mathur20/MQU/ch3/revise/tassel/

# Step1: Sort genotypes 

run_pipeline.pl -Xmx60G -SortGenotypeFilePlugin \
-inputFile /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.noSing.auto.nomiss.recode.vcf  \
-outputFile best66.old.auto.noSing.nomiss.sorted.vcf -fileType VCF

# Step2: Create UPGMA/Neighbor joining trees

count1=1

run_pipeline.pl -Xmx60G -fork$count1 \
-vcf best66.old.auto.noSing.nomiss.sorted.vcf \
-tree Neighbor -treeSaveDistance true \
-export best66.old.auto.noSing.nomissNJ${count1}

run_pipeline.pl -Xmx60G -fork$count1 \
-vcf best66.old.auto.noSing.nomiss.sorted.vcf \
-tree Neighbor -treeSaveDistance true \
-export best66.old.auto.noSing.nomissNJ${count1} -exportType Text

# Step3: Do PCA 


run_pipeline.pl -Xmx100G \
-vcf best66.old.auto.noSing.nomiss.sorted.vcf \
-PrincipalComponentsPlugin -covariance true -limitBy total_variance \
-totalVar 1 -reportEigenvectors false -endPlugin \
-export best66.old.auto.noSing.nomissPCA
