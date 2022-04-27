# Prediction of fetal genotype using parental haplotypes

Step1:

Step2: Phase parental sequencing samples using LongHap (https://github.com/stLFR/stLFR_LongHap) and stLFR sequence data.

Stpe3: Prediction of fetal genotype using parental haplotypes


1.  For SNPs that are heterozygous in the father but homozygous in mother (AAAB loci). The paternal-specific allele should be covered with 3 or more reads.
2.  The fetal inherited alleles from the father were determined by the closest informative variant algorithm (CIVA) for ABAB loci. CIVA can be used for predicting the fetal inheritance allele from the father or mother based on the nearest inferred variant within 100 kb upstream or downstream of the allele in the same haplotype block. If the upstream and downstream closest variants are from different inherited haplotypes than the allele within a 100-kb region, these alleles are not analyzed.
3.	Determine the maternal inheritance using maternal informative SNPs by the sequential probability ratio testing (SPRT) method [18], which includes two types of variants: a), variants heterozygous in the mother but homozygous in the father (ABAA loci) and b), variants heterozygous in both parents (ABAB loci) in the blocks where the fetal inherited haplotype from the father was inferred by CIVA.SPRT method tests whether the cumulative counts for Hap I and Hap II SNPs along a haplotype block are present in maternal plasma at the same or different concentrations using the SPRT and whether the SPRT reaches sufficient statistical confidence for Hap I or Hap II to be scored. 
4.	For maternal inheritance of variants at ABAA or ABAB that cannot be inferred by the SPRT method (SPRT-unclassified), the CIVA method can also be used.
