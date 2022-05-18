##########################################################################
###                          Samarth Mathur, PhD                        ###
###                        The Ohio State University                    ###
###                                                                     ###
###########################################################################
###########################################################################
###                     gadma_bootstrap.sh                        		###
###########################################################################

# bootstrapping gadma results to get CI around estimates

cd $SLURM_SUBMIT_DIR
module load use.own
#module load conda-env/easySFS-py3.7.0
module load conda-env/gadma-py3.7.0

#Only using synonymous mutations
#Do bootstrap

for i in {1..100}
do
	mkdir /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/bootstrap/run$i
	~/privatemodules/conda-env/easySFS/easySFS.py \
	-p /scratch/snyder/m/mathur20/MQU/ch3/demography/popfiles/aztx.2pop.list -a --proj 25,25 -f \
	-i /scratch/snyder/m/mathur20/MQU/ch3/synnonsyn/best66.syn.recode.vcf \
	-o /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/bootstrap/run$i \
	--prefix aztx.2pop.syn
done

for i in {1..100}
do
	cp /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/bootstrap/run$i/dadi/TX-AZ.sfs \
	/scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/bootsfs/TX-AZ_$i.sfs
done


#Get CLAIC
cd /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/R3/claic
python gadma.R3.m4.moments.py

gadma -p /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/3pop.final.m6.params

cd /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/boot

python final_model6_moments_code.py #Get CLAIC values

gadma-run_ls_on_boot_data \
-b /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/bootsfs \
-o /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/boot/ci \
-j 10 --opt powell \
-p /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/boot/ls_params \
-d /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/boot/final_model6_moments_code.forCI.py 

gadma-get_confidence_intervals /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/boot/ci/result_table.pkl --acc 3

#Using model 9 from run 2
# For every 100 spectrum from bootrstrapped data local search (Powell method) from optimal parameters was launched.


gadma-run_ls_on_boot_data -b /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/bootsfs \
-o /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/boot \
-j 1 --opt powell \
-p /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/ls_parameters \
-d /scratch/snyder/m/mathur20/MQU/ch3/demography/gadma/final/gadma.2pop.R3.m9.moments.py




