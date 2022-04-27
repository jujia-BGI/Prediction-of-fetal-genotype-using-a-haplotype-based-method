# Non-invasive prediction of fetal genotype using parental haplotypes and the cell-free DNA sequencing of maternal plasma

The input files of the workflow are the cell-free DNA sequencing reads from maternal plasma and the stLFR library sequencing of gDNA from the parental blood samples. The output files are the prediction of fetal genotypes. The required environment of Step 3 are Perl 5. For the required environment of other steps, please refer to the corresponding website.

Step1: Alignment and variant calling using GATK4 (https://github.com/broadinstitute/gatk).

sh Step1.SOAPnuke_BWA_GATK.sh

Step2: Parental haplotypes can be constructed by using stLFR sequencing on DNA obtained from their blood samples. Then phase parental sequencing samples using LongHap (https://github.com/stLFR/stLFR_LongHap).  

sh Step2.LongHap.sh

Step3: Prediction of fetal genotype using parental haplotypes.

sh Step3.haplotype_based_method.sh

The Principle of methodolog:
a.  For SNPs that are heterozygous in the father but homozygous in mother (AAAB loci). The paternal-specific allele should be covered with 3 or more reads.
b.  The fetal inherited alleles from the father were determined by the closest informative variant algorithm (CIVA) for ABAB loci. CIVA can be used for predicting the fetal inheritance allele from the father or mother based on the nearest inferred variant within 100 kb upstream or downstream of the allele in the same haplotype block. If the upstream and downstream closest variants are from different inherited haplotypes than the allele within a 100-kb region, these alleles are not analyzed.
c.	Determine the maternal inheritance using maternal informative SNPs by the sequential probability ratio testing (SPRT) method [18], which includes two types of variants: a), variants heterozygous in the mother but homozygous in the father (ABAA loci) and b), variants heterozygous in both parents (ABAB loci) in the blocks where the fetal inherited haplotype from the father was inferred by CIVA.SPRT method tests whether the cumulative counts for Hap I and Hap II SNPs along a haplotype block are present in maternal plasma at the same or different concentrations using the SPRT and whether the SPRT reaches sufficient statistical confidence for Hap I or Hap II to be scored. 
d.	For maternal inheritance of variants at ABAA or ABAB that cannot be inferred by the SPRT method (SPRT-unclassified), the CIVA method can also be used.
