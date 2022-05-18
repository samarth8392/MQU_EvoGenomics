##########################################################################
###                          Samarth Mathur, PhD                        ###
###                        The Ohio State University                    ###
###                                                                     ###
###########################################################################
###########################################################################
###                     HaplotypePhasing.sh                   		    ###
###########################################################################


cd $SLURM_SUBMIT_DIR
module load bioinfo
module load vcftools

cd /scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/vcfBychr/


#change chr to integer codes
cat best66.noSing.auto.nomiss.recode.vcf \
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
> best66.noSing.auto.nomiss.rename


#Breaking VCF by chromosome (chromosome names =1,2,3..)

for chr in $(seq 1 1 15; seq 17 1 28)
do
	vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.noSing.auto.nomiss.rename.recode.vcf \
	--chr ${chr} \
	--recode --recode-INFO-all \
	--out best66.noSing.auto.nomiss.chr${chr}
done

# Running phase by chr usping chicken recombination data (10 times)


cd /scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/phase/

for i in $(seq 1 1 10)
do
	mkdir run${i}
	for chr in $(seq 1 1 15; seq 17 1 28)
	do
		java -Xmx90g -jar ~/beagle.28Jun21.220.jar \
		burnin=200 iterations=500 phase-states=500 \
		map=/scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/chickMap/chick.chr${chr}.map \
		nthreads=52 ne=50000 \
		gt=/scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/vcfBychr/best66.noSing.auto.nomiss.chr${chr}.recode.vcf \
		out=run${i}/best66.noSing.auto.nomiss.chr${chr}
	done
done

# Compare switch errors between runs

count1=0
count2=0
for i in $(seq 1 1 10)
do
	((count1=count1+1))
	count2=0
	for j in $(seq 1 1 10)
	do
		((count2=count2+1))
		if [[ ${count2} -gt ${count1} ]]
		then
			if [[ ${i} != ${j} ]]
			then
				for chr in $(seq 1 1 15; seq 17 1 28)
				do	
					vcftools --diff-switch-error \
					--gzvcf /scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/phase/run${i}/best66.noSing.auto.nomiss.chr${chr}.vcf.gz \
					--gzdiff /scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/phase/run${j}/best66.noSing.auto.nomiss.chr${chr}.vcf.gz \
					--out /scratch/bell/mathur20/MQU/ch3/revise/phase/beagle/switch/chr${chr}.run${i}-run${j}
				done
			fi
		fi
	done
done









