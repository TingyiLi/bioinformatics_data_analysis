WORK_PATH=/data/ZX/newGSE103048
LOG_PATH=$WORK_PATH/log
REF_PATH=/data/reference_data

mkdir 07_bigwig
mkdir 01_raw_fasta
mkdir 02_fastqc
mkdir 03_bowtie2
mkdir 04_samtools_sort
mkdir 05_macs2
mkdir 06_deeptools
mkdir log

## create a filelist
cd $WORK_PATH
ls | grep 'sra$' > filelist.txt
sed -i 's/\.sra$//g' filelist.txt


## transform from SRA to fastq
cat filelist.txt | while read line
do
    echo `date`
    echo "=====start transforming SRA to fastq for $line!====="
    fastq-dump ${line}.sra --gzip -O $WORK_PATH/01_raw_fasta |tee $LOG_PATH/${line}.fastq_dump.log
done

## fastqc for fastq files
cd $WORK_PATH/01_raw_fasta
for i in *.gz
do
    echo `date`
    echo "===== start fastqc for $i! ====="
    fastqc -t 6 $i -o $WORK_PATH/02_fastqc 2>$LOG_PATH/${i}.fastqc.log
done

## map by hisat2[optional]
##cd $WORK_PATH
##cat filelist.txt | while read line
##do
#    echo `date`
#    echo "===== start mapping for $line by hisat2! ====="
#    hisat2 -p 6 -t -x $REF_PATH/index/ensembl_GRCh38/ensembl_GRCh38 -U $WORK_PATH/01_raw_fasta/${line}.fastq.gz |samtools view -bS 1>$WORK_PATH/03_hisat_map/${line}.bam 2>$LOG_PATH/${line}.hisat_map.log
#done

## step5 : alignment to hg19/ using bowtie2 to do alignment
## ~/biosoft/bowtie/bowtie2-2.2.9/bowtie2-build ~/biosoft/bowtie/hg19_index /hg19.fa ~/biosoft/bowtie/hg19_index/hg19
## cat >run_bowtie2.sh
cd $WORK_PATH
cat filelist.txt | while read line
do
    echo $id
    bowtie2 -p 8 -x $REF_PATH/index/bowtie2_hg19/hg19_dna -U $WORK_PATH/01_raw_fasta/${line}.fastq.gz -S $WORK_PATH/03_bowtie2/${line}.sam 2>$WORK_PATH/03_bowtie2/${line}.align.log;
    samtools view -bhS -q 30 $WORK_PATH/03_bowtie2/${line}.sam > $WORK_PATH/03_bowtie2/${line}.bam ## -F 1548 https://broadinstitute.github.io/picard/explain-flags.html
    -F 0x4 #remove the reads that didn't match
    #samtools sort ${id%%.*}.bam ${id%%.*}.sort ## prefix for the output
    # samtools view -bhS a.sam | samtools sort -o - ./ > a.bam
    #samtools index ${id%%.*}.sorted.bam
done
### sort file by samtools
cd $WORK_PATH/03_bowtie2
for i in *.bam
do
	echo `date`
	echo "===== start sort by samtools! ====="
	samtools sort ${i} -@ 6 -o $WORK_PATH/04_samtools_sort/sort_${i} 2> $LOG_PATH/sort_${i}.log
    samtools index ${i} ${i}.bai -@ 4 2> $LOG_PATH/${i}.index.log
done

### construct index for bam file
cd $WORK_PATH/04_samtools_sort
echo `date`
echo "===== start index for bamfiles! ====="
for i in *.bam
do
samtools index ${i} ${i}.bai -@ 4 2> $LOG_PATH/${i}.index.log
done

## call peaks by MACS2
cd $WORK_PATH/04_samtools_sort
echo `date`
echo "===== start finding peaks by MACS2! ====="
# call regular peaks for A139
macs2 callpeak -t sort_SRR5970403.bam -c sort_SRR5970404.bam --format BAM -B --name group1 -g hs -q 0.01 --outdir $WORK_PATH/05_macs2 2> $LOG_PATH/group1_peak.log
# call regular peaks for A137
#macs2 callpeak -t sort_SRR640561.bam -c sort_SRR640562.bam --format BAM -B --name T-ALL -g hs -q 0.01 --outdir $WORK_PATH/05_macs2 2> $LOG_PATH/T-ALL_peak.log
# call regular peaks for HCC70
#macs2 callpeak -t sort_SRR4029221.bam -c sort_SRR4029233.bam --format BAM -B --name H3370 -g hs -q 0.01 --outdir $WORK_PATH/05_macs2 2> $LOG_PATH/H3370_peak.log
# call regular peaks for MDAM
#macs2 callpeak -t sort_SRR4029219.bam -c sort_SRR4029231.bam --format BAM -B --name MADAM -g hs -q 0.01 --outdir $WORK_PATH/05_macs2 2> $LOG_PATH/MDAM_peak.log



## filter the bed file (filtered bed file can be annotated by R)
cd $WORK_PATH/05_macs2
for i in *.bed
do
grep '^[0-9XY]' $i | sed 's/^[0-9XY]/chr&/' > filtered_${i}
done

## draw gene profile by deeptools
cd $WORK_PATH/04_samtools_sort
echo `date`
echo "===== start draw gene profile by deeptools! ====="
for i in *.bam
do
bamCoverage -b $i -p 8 --normalizeUsingRPKM -o $WORK_PATH/06_deeptools/${i}.bw
done

cd $WORK_PATH/06_deeptools
# profile and heatmap of TSS
computeMatrix reference-point -p 6 --referencePoint TSS -b 2000 -a 2000 -S *.bw -R $REF_PATH/bed/hsa_GRCh38.bed --skipZeros -o $WORK_PATH/06_deeptools/tss.mat.gz 2> $LOG_PATH/TSS_matrix.log
plotProfile --dpi 720 -m tss.mat.gz -out $WORK_PATH/06_deeptools/tss.profile.pdf --plotFileFormat pdf --perGroup 2> $LOG_PATH/tss_profile.log
plotHeatmap --dpi 720 -m tss.mat.gz -out $WORK_PATH/06_deeptools/tss.heatmap.pdf --plotFileFormat pdf 2> $LOG_PATH/tss_heatmap.log
# profile and heatmap of gene
computeMatrix scale-regions -p 6 -S *.bw -R $REF_PATH/bed/hsa_GRCh38.bed -b 3000 -a 3000 -m 5000 --skipZeros -o $WORK_PATH/06_deeptools/gene.mat.gz 2> $LOG_PATH/gene_matrix.log
plotProfile --dpi 720 -m gene.mat.gz -out $WORK_PATH/06_deeptools/gene.profile.pdf --plotFileFormat pdf --perGroup 2> $LOG_PATH/gene_profile.log
plotHeatmap --dpi 720 -m gene.mat.gz -out $WORK_PATH/06_deeptools/gene.heatmap.pdf --plotFileFormat pdf 2>$LOG_PATH/gene_profile.log

# bam->bdg->bigwig
cd $WORK_PATH/03_bowtie2
cat $WORK_PATH/filelist.txt | while read line
do
         samtools sort -n $WORK_PATH/03_bowtie2/${line}.bam -o $WORK_PATH/07_bigwig/sorted_${line}.bam -@ 10 2> $WORK_PATH/07_bigwig/${i}.sort.log
	 macs2 pileup --extsize 200 \
		i $WORK_PATH/07_bigwig/sorted_${line}.bam \
		o $WORK_PATH/07_bigwig/${line}.bdg \
		> $WORK_PATH/07_bigwig/${line}.bgconv.log
         sort -k1,1 -k2,2n $WORK_PATH/07_bigwig/${line}.bdg > $WORK_PATH/07_bigwig/sorted_${line}.bdg
         bedGraphToBigWig $WORK_PATH/07_bigwig/sorted_${line}.bdg $WORK_PATH/hg19.chrom.sizes $WORK_PATH/07_bigwig/${line}.bigwig
done
