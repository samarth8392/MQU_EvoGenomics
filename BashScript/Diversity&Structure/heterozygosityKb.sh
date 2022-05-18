##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     hetperkb.sh                   		        	###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load vcftools

cd /scratch/bell/mathur20/MQU/ch3/revise/hetperkb/bypop/

# Calculate heterozygous sites in each population

for chr in $(seq 1 1 33)
do
	for pop in az mx wtx etx
	do
		mkdir chr${chr}
		vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.noSing.nomiss.renameChr \
		--hardy \
		--bed /scratch/bell/mathur20/MQU/ch3/revise/hetperkb/bins/chr${chr}_bin.txt \
		--keep /scratch/bell/mathur20/MQU/ch3/revise/lists/${pop}.best66.list \
		--out chr${chr}/${pop}.best66
	done
done

# Calculate heterozygous sites in each individual

cd /scratch/bell/mathur20/MQU/ch3/revise/hetperkb/byindiv/

for chr in $(seq 1 1 33)
do
	while read -a line
	do
		#mkdir chr${chr}
		vcftools --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.old.auto.poly.nomiss.renameChr.vcf \
		--hardy \
		--bed /scratch/bell/mathur20/MQU/ch3/revise/hetperkb/bins/perkb/chr${chr}_bin.txt \
		--indv ${line} \
		--out chr${chr}/${line}.best66
	done < /scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list
done


# Extract heterozygous sites from each population

#for chr in $(seq 1 1 33)
#do
#	for pop in az mx wtx etx
#	do
#		less chr${chr}/${pop}.best66.hwe | cut -f3 | cut -f2 -d "/" \
#		> chr${chr}/${pop}.best66.het
#	done
#done

# Extract heterozygous sites from each individual

#for chr in $(seq 1 1 33)
for chr in $(seq 1 1 1)
do
	while read -a line
	do
		less chr${chr}/${line}.best66.hwe | cut -f3 | cut -f2 -d "/" \
		> chr${chr}/${line}.best66.het
	done < /scratch/bell/mathur20/MQU/ch3/revise/lists/best66.samples.list
done


