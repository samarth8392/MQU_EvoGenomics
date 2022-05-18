##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     adapter_removal.sh                        		###
###########################################################################

cd $SLURM_SUBMIT_DIR

#Step0: Create folders for jobcodes and errors

#while read -a line
#do
	#mkdir /fs/scratch/PAS1533/smathur/wgr/jobcodes/per_ind/${line}
	#mkdir /fs/scratch/PAS1533/smathur/wgr/errors/per_ind/${line}
#done < /fs/scratch/PAS1533/smathur/wgr/data/sample.list

#run for each sample

while read -a line
do 
	echo "#!/bin/sh -l
#SBATCH -A gcore
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -t 300:00:00
#SBATCH --job-name=${line[0]}_adptrim
#SBATCH -e ${line[0]}_adptrim
#SBATCH -o ${line[0]}_adptrim

cd $SLURM_SUBMIT_DIR
module load bioinfo
module load fastqc/0.11.8
#module load trimmomatic

# Step1: FASTQC on raw reads

# unzip reads

#cd /scratch/snyder/m/mathur20/MQU/ch3/reads/r4/
#unzip *zip


# FASTQC on raw reads

cd /scratch/snyder/m/mathur20/MQU/ch3/align/raw
mkdir ${line}

fastqc -o ${line}/ \
/scratch/snyder/m/mathur20/MQU/ch3/reads/r4/${line[0]}_R1.fastq \
/scratch/snyder/m/mathur20/MQU/ch3/reads/r4/${line[0]}_R2.fastq \


#
### Step2: Adapter removal 
##
##java -jar /apps/trimmomatic/0.38/trimmomatic-0.38.jar \
##PE /scratch/snyder/m/mathur20/MQU/ch3/reads/r4/${line[0]}_R1.fastq /scratch/snyder/m/mathur20/MQU/ch3/reads/r4/${line[0]}_R2.fastq \
###/scratch/snyder/m/mathur20/MQU/ch3/align/trimmed/${line[0]}.paired /scratch/snyder/m/mathur20/MQU/ch3/align/trimmed/${line[0]}.unpaired \
##/scratch/snyder/m/mathur20/MQU/ch3/align/trimmed/${line[0]}.paired /scratch/snyder/m/mathur20/MQU/ch3/align/trimmed/${line[0]}.unpaired \
##LEADING:20 TRAILING:20 MINLEN:30 \
##ILLUMINACLIP:/apps/trimmomatic/0.38/adapters/NexteraPE-PE.fa:2:40:10"
> /scratch/snyder/m/mathur20/MQU/ch3/jobcodes/mapping/align/${line[0]}_chick_alignment.sh
done < /scratch/snyder/m/mathur20/MQU/ch3/reads/best66.samples.list


# Submit jobs #

while read -a line
do
	cd /scratch/snyder/m/mathur20/MQU/ch3/per_ind/${line}
	sbatch  /scratch/snyder/m/mathur20/MQU/ch3/jobcodes/mapping/align/${line[0]}_chick_alignment.sh
done < /scratch/snyder/m/mathur20/MQU/ch3/reads/best66.samples.list