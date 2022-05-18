##########################################################################
###                          Samarth Mathur, PhD                        ###
###                        The Ohio State University                    ###
###                                                                     ###
###########################################################################
###########################################################################
###                     mutation_age.sh                   		        ###
###########################################################################

# Get genotypes for mutations with ages

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load vcftools


cd /scratch/bell/mathur20/MQU/ch3/revise/age/geno/

# Reformat chicken recombination maps

for i in $(seq 1 1 15; seq 17 1 28)
do
    Rscript /scratch/snyder/m/mathur20/MQU/ch3/jobcodes/RScripts/recomb.R \
    -f /scratch/snyder/m/mathur20/MQU/ch3/mut_age/run4/recMap/chicken.chr$i.recombfile \
    -o /scratch/snyder/m/mathur20/MQU/ch3/mut_age/run4/recMap/chicken.chr$i.recMap

    cat recMap/chicken.chr$i.recMap | awk '{ print $2, $3, $4}' | tail -n+2 | sed '1i Position(bp)  Rate(cM/Mb) Map(cM)' \
    > recMap/chicken.chr$i.recMap.Final
done


# Rename sites
for i in del weakdel syn
do
cat /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn/sites/$i.sites \
| sed -e 's/NC_006088.5/1/' \
-e 's/NC_006089.5/2/' \
-e 's/NC_006090.5/3/' \
-e 's/NC_006091.5/4/' \
-e 's/NC_006092.5/5/' \
-e 's/NC_006093.5/6/' \
-e 's/NC_006094.5/7/' \
-e 's/NC_006095.5/8/' \
-e 's/NC_006096.5/9/' \
-e 's/NC_006097.5/10/' \
-e 's/NC_006098.5/11/' \
-e 's/NC_006099.5/12/' \
-e 's/NC_006100.5/13/' \
-e 's/NC_006101.5/14/' \
-e 's/NC_006102.5/15/' \
-e 's/NC_006103.5/16/' \
-e 's/NC_006104.5/17/' \
-e 's/NC_006105.5/18/' \
-e 's/NC_006106.5/19/' \
-e 's/NC_006107.5/20/' \
-e 's/NC_006108.5/21/' \
-e 's/NC_006109.5/22/' \
-e 's/NC_006110.5/23/' \
-e 's/NC_006111.5/24/' \
-e 's/NC_006112.4/25/' \
-e 's/NC_006113.5/26/' \
-e 's/NC_006114.5/27/' \
-e 's/NC_006115.5/28/' \
-e 's/NC_028739.2/30/' \
-e 's/NC_028740.2/31/' \
-e 's/NC_006119.4/32/' \
-e 's/NC_008465.4/33/' \
> sites/$i.sites.rename
done

for sites in del syn
do
    for i in $(seq 1 1 15; seq 17 1 28)
    do
      vcftools --vcf /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/input/beagle/best66.chr$i.beagle.phase.vcf \
     --recode --recode-INFO-all \
     --positions sites/$sites.sites.rename \
     --out vcf/$sites/best66.chr$i.$sites
    
     /scratch/snyder/m/mathur20/MQU/ch3/softwares/geva/geva_v1beta \
     --vcf /scratch/snyder/m/mathur20/MQU/ch3/structure/fineS/run2/input/beagle/best66.chr$i.beagle.phase.vcf \
     --rec 1e-8 \
     --out out/best66.chr$i
done

#Get batch files
for i in $(seq 1 1 15; seq 17 1 28)
do
    cat out/best66.chr$i.marker.txt | awk '{ print $3;}' | tail -n+2 \
    > batch/chr.$i.BATCH.txt
done

ne=(18810.3 19074.2 19436.7 22829.8 17586.6 17698.2 5.53973e+09 4365.5 22959 18282.2 30285.7 11673.2 17959.4 15753.2 22590.6 8523.17 41624.8 7917.44 16126.4 2.81728e+10 1.17484e+11 7841.53 7110.76 4294.34 49592.2 5861.5 8753.43)
count=0
for i in $(seq 1 1 15; seq 17 1 28)
do
    /scratch/snyder/m/mathur20/MQU/ch3/softwares/geva/geva_v1beta \
    -i out/best66.chr$i.bin \
    -o weakdel/best66.chr$i \
    --positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/run2/batch/del/chr.$i.BATCH.txt \
    --Ne ${ne[$count]} --mut 3.13597726022324e-09 \
    --hmm /scratch/snyder/m/mathur20/MQU/ch3/softwares/geva/hmm/hmm_initial_probs.txt \
    /scratch/snyder/m/mathur20/MQU/ch3/softwares/geva/hmm/hmm_emission_probs.txt
    count=$[$count+1]
done

for i in $(seq 1 1 15; seq 17 1 28)
do
    Rscript /scratch/snyder/m/mathur20/MQU/ch3/softwares/geva/estimate.R \
    out2/del/best66.chr$i.del.pairs.txt 100000 
done

cd /scratch/snyder/m/mathur20/MQU/ch3/mut_age/aztwx/weakdel

#Shared v Private

ne=(18810 19074 19436 22829 17586 17698 5.53973e+09 4365 22959 18282 30285 11673 17959 15753 22590 8523 41624 7917 16126 2.81728e+10 1.17484e+11 7841 7110 4294 49592 5861 8753)
count=0

#del
for i in $(seq 1 1 15; seq 17 1 28)
do
while read -a line
do
    /scratch/snyder/m/mathur20/MQU/ch3/softwares/geva/geva_v1beta \
    -i /scratch/snyder/m/mathur20/MQU/ch3/mut_age/run4/out/best66.chr${line[0]}.bin \
    -o wtx.syn.chr${line[0]} \
    --positions /scratch/snyder/m/mathur20/MQU/ch3/mut_age/aztwx/batchfiles/wtx.syn.chr${line[0]}.BATCH.txt.2 \
    --Ne ${line[1]} --mut 3.13597726022324e-09 \
    --hmm /scratch/snyder/m/mathur20/MQU/ch3/softwares/geva/hmm/hmm_initial_probs.txt \
    /scratch/snyder/m/mathur20/MQU/ch3/softwares/geva/hmm/hmm_emission_probs.txt
    count=$[$count+1]
done < /scratch/snyder/m/mathur20/MQU/ch3/mut_age/Ne_chr

while read -a line
do
    Rscript /scratch/snyder/m/mathur20/MQU/ch3/softwares/geva/estimate.R \
    az.weakdel.chr${line[0]}.pairs.txt ${line[1]}
done < /scratch/snyder/m/mathur20/MQU/ch3/mut_age/Ne_chr



