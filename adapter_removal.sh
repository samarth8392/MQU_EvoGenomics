#!/bin/sh -l
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 150:00:00
#SBATCH --job-name=adapter_removal
#SBATCH -e adapter_removal
#SBATCH -o adapter_removal

##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###     Date Created: 01/26/21                  Last Modified: 01/26/21 ###
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
#SBATCH -A PAS1533
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 150:00:00
#SBATCH --job-name=${line}_adptrem
#SBATCH -e ${line}_adptrem
#SBATCH -o ${line}_adptrem

cd $SLURM_SUBMIT_DIR
module load fastqc/0.11.8
#module load trimmomatic

# Step1: FASTQC on raw reads

cd /fs/scratch/PAS1533/smathur/wgr/read_preprocess/fastqc/trimmed/
mkdir ${line}

fastqc -o ${line}/ \
/fs/scratch/PAS1533/smathur/wgr/read_preprocess/trimmomatic/trimmed/paired/${line}*R1*fastq.paired \
/fs/scratch/PAS1533/smathur/wgr/read_preprocess/trimmomatic/trimmed/paired/${line}*R2*fastq.paired

#cd /fs/scratch/PAS1533/smathur/wgr/read_preprocess/fastqc/raw/${line}
#unzip *zip
#
### Step2: Adapter removal 
##
##java -jar /apps/trimmomatic/0.38/trimmomatic-0.38.jar \
##PE /fs/scratch/PAS1533/smathur/wgr/data/raw/${line[0]} /fs/scratch/PAS1533/smathur/wgr/data/raw/${line[1]} \
###/fs/scratch/PAS1533/smathur/wgr/data/trimmed/${line[0]}.paired /fs/scratch/PAS1533/smathur/wgr/data/trimmed/${line[0]}.unpaired \
##/fs/scratch/PAS1533/smathur/wgr/data/trimmed/${line[1]}.paired /fs/scratch/PAS1533/smathur/wgr/data/trimmed/${line[1]}.unpaired \
##LEADING:20 TRAILING:20 MINLEN:30 \
##ILLUMINACLIP:/apps/trimmomatic/0.38/adapters/NexteraPE-PE.fa:2:40:10" \
> /fs/scratch/PAS1533/smathur/wgr/jobcodes/per_ind/${line}/${line}_adptrem.sh
done < /fs/scratch/PAS1533/smathur/wgr/data/sample.list
#
#
#while read -a line
#do
#	java -jar /apps/trimmomatic/0.38/trimmomatic-0.38.jar \
#	PE /fs/scratch/PAS1533/smathur/wgr/data/raw/${line[0]} /fs/scratch/PAS1533/smathur/wgr/data/raw/${line[1]} \
#	/fs/scratch/PAS1533/smathur/wgr/data/trimmed/${line[0]}.paired /fs/scratch/PAS1533/smathur/wgr/data/trimmed/${line[0]}.unpaired \
#	/fs/scratch/PAS1533/smathur/wgr/data/trimmed/${line[1]}.paired /fs/scratch/PAS1533/smathur/wgr/data/trimmed/${line[1]}.unpaired \
#	LEADING:20 TRAILING:20 MINLEN:30 -threads 20 \
#	ILLUMINACLIP:/fs/scratch/PAS1533/smathur/wgr/read_preprocess/trimmomatic/NexteraPE-PE.fa:2:40:10 
#done < /fs/scratch/PAS1533/smathur/wgr/data/readnames.list