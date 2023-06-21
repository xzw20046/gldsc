source activate ldsc

ld_ref=/Users/zewei/local/ldsc_new_ref/baseline.10
name=/Users/zewei/local/str-HDL/sim_gwas/DHS_newpolyc0.01h0.4e3N20000/chr10simgwas_test1.sumstat

python /Users/zewei/Desktop/software/ldsc/ldsc.py --h2 $name --ref-ld $ld_ref --w-ld /Users/zewei/local/ldsc/weights_hm3_no_hla/weights.10 --overlap-annot --frqfile /Users/zewei/local/ldsc/1000G_frq/1000G.mac5eur.10 --out ~/Desktop/HKUgrad/G-LDSC/gldsc/gldsc/simulation_code/ldsc_result_test --chisq-max 1000000 --print-coefficients 
