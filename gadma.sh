#!/bin/sh -l
#SBATCH -A gcore
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 300:00:00
#SBATCH --job-name=gadma.3p2
#SBATCH -e gadma.3p2
#SBATCH -o gadma.3p2

##########################################################################
###                             Samarth Mathur                          ###
###                        PhD Candidate, DeWoody Lab                   ###
###                           Purdue University                         ###
###                                                                     ###
###     Date Created: 08/17/20                  Last Modified: 08/17/20 ###
###########################################################################
###########################################################################
###                     gadma.sh                        				###
###########################################################################

cd $SLURM_SUBMIT_DIR

module load use.own
module load conda-env/gadma-py3.7.0
#module load conda-env/easySFS-py3.7.0
#module load bioinfo
#module load vcftools

#Only using synonymous mutations
# (N=454,896 out of a  22,629,734 Sites)
#2-pop model (TX,AZ)

#Step1: To convert a VCF (.vcf) file into a SFS (.sfs) file
cd /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/2pop/

~/privatemodules/conda-env/easySFS/easySFS.py \
-p /scratch/snyder/m/mathur20/MQU/ch3/demography/popfiles/aztx.2pop.list -a --proj 25,25 -f \
-i /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn/best66.syn.recode.vcf \
-o sfs/ \
--prefix aztx.2pop.syn


 
cd /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/R3/2pop/
 # CREATE PARAM FILE #

count=1
for i in 1 2 3
do
	for j in 1 2 3 
	do
		#for k in 1 2 3
		# do
		echo \
"#GADMA_Run3 $count
# MQU  2 Pop: AZ, TX
Output directory : /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/R3/2pop/models/m$count
Input file : /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/only_syn/2pop/sfs/dadi/TX-AZ.sfs
Population labels : TX, AZ
Projections : 25, 25
Linked SNP's : True
Theta0 :  0.2838653
Time for generation : 1
Use moments or dadi : moments
Multinom : True
#Pts : 25, 35, 45
Lower bounds : 1e-2 0 0 0 
Upper bounds : 5 5 5 1
Parameter identifiers: N, T, m, s
Only sudden : False
Initial structure : $i, $j
Final Structure: $i, $j 
Relative parameters : False
No migrations : false
Name of local optimization : optimize_log
Number of repeats : 20
Number of processes : 20
Draw models every N iteration : 100
Units of time in drawing : thousand years
Upper bound of first split: 50000
#Mean mutation strength : 0.2
#Const for mutation strength : 1.0
#Mean mutation rate : 0.2
#Const for mutation rate : 1.0
Epsilon : 0.001
Stop iteration : 1000" > params/gadma.run$count.params 
		count=$[$count+1]
	#done 
done
done

# Create Job file
for i in `seq 1 1 9`
do
	mkdir models/m$i
	echo "#!/bin/sh -l
#SBATCH -A fnrquail
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 14-00:00:00
#SBATCH --job-name=gadma.2pop.run$i
#SBATCH -e gadma.2pop.run$i
#SBATCH -o gadma.2pop.run$i

module load use.own
module load conda-env/gadma-py3.7.0

gadma -p /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/R3/2pop/params/gadma.run$i.params" \
> jobs/gadma.run$i.job
done

# Submit jobs
for i in `seq 1 1 9`
do
	cd /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/R3/2pop/errors/ 
	sbatch /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/R3/2pop/jobs/gadma.run$i.job
done

#Grab results

cd /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/R3/2pop/results/
for i in `seq 1 1 9`
do 
	cp ../models/m$i/best_logLL.png ./m$i.best.logLL.png
done

