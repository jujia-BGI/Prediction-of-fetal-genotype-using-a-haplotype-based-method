SOAPnuke filter -1 sample_1.fq.gz -2 sample_2.fq.gz -l 5 -q 0.5 -n 0.1 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -t 0,0,0,0 -Q 2 -G 2 -T 6 --seqType 1 -o sample_1.clean.fq.gz -D sample_2.clean.fq.gz && echo "** SOAPnuke done **" && \
bwa  mem -t 4  -M -Y  -R  "@RG\tID:RGID\tPL:COMPLETE\tSM:sample" hg19.fasta sample_1.clean.fq.gz  sample_2.clean.fq.gz | samtools view -Sb -o  sample.bam - && echo "** BWA MEM done **" && \
samtools sort  -@ 8 -O bam -o  sample.sorted.bam sample.bam && echo "** bam sorted done **" && \
gatk MarkDuplicates REMOVE_DUPLICATES=ture -I sample.sorted.bam -O sample.sorted.rmdup.bam -M sample_metrics.metrics.txt &&  echo "** MarkDuplicates done **" && \
gatk BaseRecalibrator  -R hg19.fasta -I sample.sorted.rmdup.bam --knownSites dbsnp_138.hg19.vcf.gz --knownSites Mills_and_1000G_gold_standard.indels.hg19.vcf.gz --knownSites 1000G_phase1.indels.hg19.vcf.gz  -O sample.recal_data.table &&  echo "** BaseRecalibrator done **" && \
gatk  ApplyBQSR  -R hg19.fasta -I sample.sorted.rmdup.bam  --bqsr-recal-file  sample.recal_data.table  -O sample.sorted.rmdup.BQSR.bam && echo "** PrintReads done **" && \
samtools index sample.sorted.rmdup.BQSR.bam && echo "** index done **" && \
gatk HaplotypeCaller  --emit-ref-confidence GVCF -R hg19.fasta -I sample.sorted.rmdup.BQSR.bam -D dbsnp_138.b37.vcf.gz  -O sample.HC.g.vcf.gz  echo "** HaplotypeCaller done **"

###family samples gvcf merge
gatk CombineGVCFs  -R hg19.fasta -V sample1.HC.g.vcf.gz -V sample2.HC.g.vcf.gz -V sample3.HC.g.vcf.gz -O family.merged.g.vcf.gz && "** CombineGVCFs done **" && \
gatk GenotypeGVCFs  -R hg19.fasta -V family.Merged.g.vcf.gz  --allow-old-rms-mapping-quality-annotation-data  -O family.merged.vcf.gz && echo "** GenotypeGVCFs done **" && \
gatk VariantRecalibrator \
   -R hg19.fasta \
   -V family.merged.vcf.gz \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
   --resource:omini,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=6.0 dbsnp_146.hg38.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
   -mode SNP \
   --max-gaussians 4 \
   -tranche 100.0 -tranche 99.9 -tranche 95.0 -tranche 90.0 \
   -O family.snps.recal \
   --tranches-file family.HC.snps.tranches \
   --rscript-file family.HC.snps.plots.R && echo "** SNP VariantRecalibrator done **"
gatk ApplyVQSR \
   -R hg19.fasta \
   -V family.merged.vcf.gz \
   --truth-sensitivity-filter-level 99.5 \  
   --tranches-file family.HC.snps.tranches \
   --recal-file family.snps.recal \
   -mode SNP \
   -O family.merged.snp.vcf.gz && echo "** SNP ApplyRecalibration done **"


