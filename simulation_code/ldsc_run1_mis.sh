source activate ldsc

ld_ref=/Users/zewei/local/ldsc_new_ref/baselineUKB.10
name1=DHS_newN

for N in {500,2500,5000,10000,20000,50000}
do
for E in {1,2,3}
do
for H in {0.04,0.4}
do
name=$name1"c0.05h"$H"e"$E"N"$N

rm /Users/zewei/local/str-HDL/result/$name/ldsc_h2_mis.txt
rm /Users/zewei/local/str-HDL/result/$name/ldsc_inter_mis.txt

out=/Users/zewei/local/str-HDL/result/$name
in=/Users/zewei/local/str-HDL/sim_gwas/$name

for simu in {1..50}
do

python /Users/zewei/Desktop/software/ldsc/ldsc.py --h2 $in/chr10simgwas_test$simu.sumstat --ref-ld $ld_ref --w-ld /Users/zewei/local/ldsc/weights_hm3_no_hla/weights.10 --overlap-annot --frqfile /Users/zewei/local/ldsc/1000G_frq/1000G.mac5eur.10 --out $out/ldsc_mis_$simu --chisq-max 1000000 --print-coefficients 

cat $out/ldsc_mis_$simu.log | grep "Total Observed scale h2:" >> $out/ldsc_h2_mis.txt

cat $out/ldsc_mis_$simu.log | grep "Intercept: " >> $out/ldsc_inter_mis.txt

done
done
done
done
