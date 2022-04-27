sam="JK-10" #sample name
cfDNA="TDP1807015083_cfDNA" #maternal cfDNA name in family merged vcf
M_stLFR="TDP1807015084_M_stLFR" #maternal gDNA name in family merged vcf
P_stLFR="TDP1807015714_P_stLFR" #paternal gDNA name in family merged vcf
ff=0.15 #fetal fraction

for i in `seq 1 22` X ;do
        #input: <fetal fraction> <family merged vcf> <maternal haplotypes constructed by stLFR sequencing> <paternal haplotypes constructed by stLFR sequencing> <maternal cfDNA name in family merged vcf> <maternal gDNA name in family merged vcf> <paternal gDNA name in family merged vcf> 
        perl Predict_haplotype_genotype.pl  $ff example/chr${i}.filter.vcf example/${M_stLFR}_chr${i}.vcf.gz  example/${P_stLFR}_chr${i}.vcf.gz $cfDNA  $M_stLFR  $P_stLFR > chr${i}.longhap.step1
        #input: <result of last step> <length in closest informative variant algorithm> <fetal fraction> <cfDNA name in vcf>
        perl Haplotype_adjustment_flank_length_SPRT_C5_haplotype_same_100kb.pl chr${i}.longhap.step1 100000 $ff $cfDNA > chr${i}.longhap.step2
done

