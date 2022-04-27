sam="JK-10" #sample name
cfDNA="TDP1807015083_cfDNA" #maternal cfDNA name in family merged vcf
M_stLFR="TDP1807015084_M_stLFR" #maternal stLFR name in family merged vcf
P_stLFR="TDP1807015714_P_stLFR" #paternal stLFR name in family merged vcf
ff=0.15 #fetal fraction

for i in `seq 1 22` X ;do
        #input: <fetal fraction> <family merged vcf> <maternal haplotypes constructed by stLFR> <paternal haplotypes constructed by stLFR> <cfDNA name in vcf:TDP1807015083_cfDNA> <maternal name in vcf:TDP1807015084_M_stLFR> <paternal name in vcf:TDP1807015714_P_stLFR> 
        perl Predict_haplotype_genotype.pl  $ff example/chr${i}.filter.vcf example/${M_stLFR}_chr${i}.vcf.gz  example/${P_stLFR}_chr${i}.vcf.gz $cfDNA  $M_stLFR  $P_stLFR > chr${i}.longhap.step1
        #input: <result of last step> <length in closest informative variant algorithm> <fetal fraction> < cfDNA name in vcf>
        perl Haplotype_adjustment_flank_length_SPRT_C5_haplotype_same_100kb.pl chr${i}.longhap.step1 100000 $ff $cfDNA > chr${i}.longhap.step2
done

