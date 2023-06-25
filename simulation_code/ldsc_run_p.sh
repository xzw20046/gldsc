source activate ldsc

ld_ref=/Users/zewei/local/ldsc_new_ref/baseline.10
name1=new_evalup

for E in {0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3}
do
name=$name1"c0.05h0.4e"$E

rm /Users/zewei/local/str-HDL/result/$name/ldsc_h2_v2.txt
rm /Users/zewei/local/str-HDL/result/$name/ldsc_inter_v2.txt

out=/Users/zewei/local/str-HDL/result/$name
in=/Users/zewei/local/str-HDL/sim_gwas/$name

for simu in {1..200}
do

python /Users/zewei/Desktop/software/ldsc/ldsc.py --h2 $in/chr10simgwas_test$simu.sumstat --ref-ld $ld_ref --w-ld /Users/zewei/local/ldsc/weights_hm3_no_hla/weights.10 --overlap-annot --frqfile /Users/zewei/local/ldsc/1000G_frq/1000G.mac5eur.10 --out $out/ldsc_v2_$simu --chisq-max 1000000 --print-coefficients 

cat $out/ldsc_v2_$simu.log | grep "Total Observed scale h2:" >> $out/ldsc_h2_v2.txt

cat $out/ldsc_v2_$simu.log | grep "Intercept: " >> $out/ldsc_inter_v2.txt

done

done

