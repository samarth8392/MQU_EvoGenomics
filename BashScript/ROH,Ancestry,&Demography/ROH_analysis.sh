##########################################################################
###                          Samarth Mathur, PhD                     	###
###                        The Ohio State University                 	###
###                                                                     ###
###########################################################################
###########################################################################
###                     ROH_analysis.sh                        			###
###########################################################################

cd $SLURM_SUBMIT_DIR
module purge
module load bioinfo
module load plink/1.90b6.4
#module load vcftools


# VCF to plink format
cd /scratch/bell/mathur20/MQU/ch3/revise/roh/

plink --vcf /scratch/bell/mathur20/MQU/ch3/revise/vcfs/best66.noSing.auto.nomiss.recode.vcf \
--allow-extra-chr -out best66.auto.noSing.nomiss


# Test runs (always 100kb kb)
# see: https://www.cog-genomics.org/plink/1.9/ibd#homozyg

#--homozyg-snp <min SNP count> # Default = 100
#--homozyg-kb <min length> # Default = 1000
#--homozyg-density <max inverse density (kb/SNP)> # Default = 50
#--homozyg-gap <max internal gap kb length> # Default = 1000
#--homozyg-het <max hets> # Default = Unlimited
#--homozyg-window-snp <scanning window size> # Default = 50
#--homozyg-window-het <max hets in scanning window hit> # Default = 1 
#--homozyg-window-missing <max missing calls in scanning window hit> # Default = 5 
#--homozyg-window-threshold <min scanning window hit rate> # Default = 0.05

for snp in 50 100 500 
do
	for den in 10 50 100
	do
		for gap in 500 1000 5000
		do
			for het in 0 2 5
			do
				for wsnp in 10 50 100
				do
					for whet in 0 1 5
					do
						for wmiss in 0 5 10
						do
							for wthresh in 0.01 0.05 0.1
							do
								plink --bim best66.noSing.auto.nomiss.bim \
								--bed best66.noSing.auto.nomiss.bed \
								--fam best66.noSing.auto.nomiss.fam \
								--allow-extra-chr --no-sex --no-parents --no-pheno \
								--homozyg-snp $snp \
								--homozyg-kb 100 \
								--homozyg-density $den \
								--homozyg-gap $gap \
								--homozyg-het $het \
								--homozyg-window-snp $wsnp \
								--homozyg-window-het $whet \
								--homozyg-window-missing $wmiss \
								--homozyg-window-threshold $wthresh \
								--out tryruns1/best66.${snp}snp.100kb.${den}den.${gap}gap.${het}het.${wsnp}wsnp.${whet}whet.${wmiss}wmiss.${wthresh}wthresh
							done
						done
					done
				done
			done
		done
	done
done

# Try runs 2
# Using only the sensitive parameters (all +ive): snp, het, wsnp, whet, wmiss, (-) gap

for snp in 50 75 100
do
	for gap in 100 200 500
	do		
		for het in 1 2 3
		do
			for wsnp in 10 20 30 40 50
			do
				for whet in 1 2 3
				do	
					for wmiss in 1 2 5
					do
						plink --bim best66.noSing.auto.nomiss.bim \
						--bed best66.noSing.auto.nomiss.bed \
						--fam best66.noSing.auto.nomiss.fam \
						--allow-extra-chr --no-sex --no-parents --no-pheno \
						--homozyg-snp $snp \
						--homozyg-kb 100 \
						--homozyg-density 50 \
						--homozyg-gap $gap \
						--homozyg-het $het \
						--homozyg-window-snp $wsnp \
						--homozyg-window-het $whet \
						--homozyg-window-missing $wmiss \
						--homozyg-window-threshold 0.05 \
						--out tryruns2/best66.${snp}snp.100kb.${gap}gap.${het}het.${wsnp}wsnp.${whet}whet.${wmiss}wmiss
					done
				done
			done
		done
	done
done
#
# Try runs 3
# Using only the sensitive parameters: snp (+), het (+), wmiss (+)

for snp in 25 50 60 70 75 80
do
	for het in 0 1 2 3
	do	
		for wmiss in 1 2 3
		do
			plink --bim best66.noSing.auto.nomiss.bim \
			--bed best66.noSing.auto.nomiss.bed \
			--fam best66.noSing.auto.nomiss.fam \
			--allow-extra-chr --no-sex --no-parents --no-pheno \
			--homozyg-snp $snp \
			--homozyg-kb 100 \
			--homozyg-density 50 \
			--homozyg-gap 1000 \
			--homozyg-het $het \
			--homozyg-window-snp 25 \
			--homozyg-window-het 2 \
			--homozyg-window-missing $wmiss \
			--homozyg-window-threshold 0.05 \
			--out tryruns3/best66.${snp}snp.100kb.${het}het.${wmiss}wmiss
		done
	done
done

# Try runs 4
# Using only the sensitive parameters: snp 

for snp in 70 80 90 100 110 120 130 140 150
do
	plink --bim best66.noSing.auto.nomiss.bim \
	--bed best66.noSing.auto.nomiss.bed \
	--fam best66.noSing.auto.nomiss.fam \
	--allow-extra-chr --no-sex --no-parents --no-pheno \
	--homozyg-snp $snp \
	--homozyg-kb 100 \
	--homozyg-density 50 \
	--homozyg-gap 1000 \
	--homozyg-het 1 \
	--homozyg-window-snp 25 \
	--homozyg-window-het 2 \
	--homozyg-window-missing 2 \
	--homozyg-window-threshold 0.05 \
	--out tryruns4/best66.${snp}snp
done


# Final run
# --homozyg-window-het 2, --homozyg-snp 50, --homozyg-kb 100, --homozyg-window-snp 20 --homozyg-het 2 \
# --homozyg-density 50 --homozyg-gap 1000 --homozyg-window-missing 2 --homozyg-window-threshold 0.05

plink --bim best66.noSing.auto.nomiss.bim \
--bed best66.noSing.auto.nomiss.bed \
--fam best66.noSing.auto.nomiss.fam \
--allow-extra-chr --no-sex --no-parents --no-pheno \
--homozyg-window-het 2, --homozyg-snp 50, --homozyg-kb 100 \
--homozyg-window-snp 20 --homozyg-het 2 \
--homozyg-density 50 --homozyg-gap 1000 --homozyg-window-missing 2 --homozyg-window-threshold 0.05 \
--out best66.noSing.auto.nomiss.final






