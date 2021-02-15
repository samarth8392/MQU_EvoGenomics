##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 08/22/20                  Last Modified: 09/02/20 ###
###########################################################################
###########################################################################
###                     depth&breadth.sh                       			###
###########################################################################

# Get coverage statistics for BAM files

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load bwa
module load picard-tools
module load bedops
module load GATK/3.6.0
module load samtools

cd /scratch/snyder/m/mathur20/MQU/ch3/align/chick/final_bam/stats/

while read -a line
do
	samtools depth -a ../dedup/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
	| awk '{c++;s+=$3}END{print s/c}' \
	> ${line[0]}_chick_meandepth.txt

	samtools depth -a ../dedup/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
	| 	awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' \
	> ${line[0]}_chick_1x_breadth.txt

	samtools depth -a ../dedup/${line[0]}.chick.sorted.dedup.realigned.fixmate.bam \
	| 	awk '{c++; if($3>=5) total+=1}END{print (total/c)*100}' \
	> ${line[0]}_chick_5x_breadth.txt
done < /scratch/snyder/m/mathur20/MQU/ch3/reads/sample.list.all