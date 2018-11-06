#!/bin/bash
# Define work path and set up folders
WORK_PATH=/data/ZX/GSE84515_ATAC_bc
LOG_PATH=$WORK_PATH/log
REF_PATH=/data/reference_data

 mkdir 01_fastq
 mkdir 02_cutadapt
 mkdir 03_bowtie2
 mkdir 04_samtools_sort
 mkdir 05_picard
 mkdir 06_samtools_filter
 mkdir 07_macs2
 mkdir 08_homer
 mkdir 09_results
 mkdir 10_bigwing
 mkdir log


# 1. fastq-dump extract fastq files from sra files
cd $WORK_PATH
cat $WORK_PATH/filelist.txt | while read line
do
    fastq-dump --split-3 --gzip -O $WORK_PATH/01_fastq $line.sra >$LOG_PATH/$line.dump.log 2>&1
done

# 2. trim adaptors by using cutadapt
cd $WORK_PATH/01_fastq
cat $WORK_PATH/filelist.txt | while read line
do
    cutadapt -a CTGTCTCTTATACACATCT \
             -a AGATGTGTATAAGAGACAG  \
             -A AGATGTGTATAAGAGACAG \
             -A CTGTCTCTTATACACATCT \
             -O 5  \
             -o $WORK_PATH/02_cutadapt/trimmed_${line}_1.fastq.gz \
             -p $WORK_PATH/02_cutadapt/trimmed_${line}_2.fastq.gz \
                ${line}_1.fastq.gz ${line}_2.fastq.gz \
                1>$LOG_PATH/$line.cutadapt.log 2>&1
done

# 3. bowtie2 alignment
cd $WORK_PATH/02_cutadapt
cat $WORK_PATH/filelist.txt | while read line
do
    bowtie2 -p 20 -X 2000 -x $REF_PATH/index/bowtie2_hg19/hg19_dna \
            -1 trimmed_${line}_1.fastq.gz \
            -2 trimmed_${line}_2.fastq.gz \
             2>$LOG_PATH/$line.bowtie2.log | \
             samtools view -@ 10 -bT $REF_PATH/index/bowtie2_hg19/hg19_dna.fa \
                           -o $WORK_PATH/03_bowtie2/$line.bam \
                            >$LOG_PATH/$line.bamconv.log 2>&1
done

# 4. sort bam files by using samtools
cd $WORK_PATH/03_bowtie2
cat $WORK_PATH/filelist.txt | while read line
do
    samtools sort -m 1024M -@ 20 -o $WORK_PATH/04_samtools_sort/sorted_$line.bam $line.bam \
                   >$LOG_PATH/$line.sort.log 2>&1
done

# 5. deduplication of bam files by using picard
cd $WORK_PATH/04_samtools_sort
cat $WORK_PATH/filelist.txt | while read line
do
    java -jar /home/RNAseq_tool/picard/picard.jar MarkDuplicates \
          I=sorted_$line.bam \
          O=$WORK_PATH/05_picard/dedup_${line}.bam \
          M=$LOG_PATH/${line}.dedup.log \
          REMOVE_DUPLICATES=true
done

# 6. filter out unwanted reads by using samtools
cd $WORK_PATH/05_picard
cat $WORK_PATH/filelist.txt | while read line
do
    samtools index dedup_${line}.bam
    samtools view -b -@ 20 -q 30 -f 0x2 \
                  -o $WORK_PATH/06_samtools_filter/filtered_$line.bam dedup_${line}.bam \
                     chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 \
                     chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX \
                   >$LOG_PATH/$line.filter.log 2>&1
done

# 7. peak calling by using macs2
cd $WORK_PATH/06_samtools_filter
cat $WORK_PATH/filelist.txt | while read line
do
    macs2 callpeak --nomodel --nolambda --keep-dup all --call-summits \
                   -B -n $line \
                   -t filtered_${line}.bam \
                   --outdir $WORK_PATH/07_macs2 \
                     >$LOG_PATH/$line.callpeaks.log 2>&1
done


# 8. Motif and TF finding with HOMER
cd $WORK_PATH/08_homer
cat $WORK_PATH/filelist.txt | while read line
do
	mkdir $line
	findMotifsGenome.pl $WORK_PATH/07_macs2/${line}_summits.bed \
	hg19 $WORK_PATH/08_homer/$line \
	>$LOG_PATH/$line.homer.log 2>&1
done

# 9. bedGraphToBigWig
cd $WORK_PATH/09_results
cat $WORK_PATH/filelist.txt | while read line
do
         sort -k1,1 -k2,2n $WORK_PATH/07_macs2/${line}_treat_pileup.bdg > $WORK_PATH/09_results/sorted_${line}.bdg
         bedGraphToBigWig $WORK_PATH/09_results/sorted_${line}.bdg $WORK_PATH/hg19.chrom.sizes $WORK_PATH/09_results/${line}.bigwig
done

# 10. bam->bdg->bigwig
cd $WORK_PATH/03_bowtie2
cat $WORK_PATH/filelist.txt | while read line
do
         samtools sort -n $WORK_PATH/03_bowtie2/${line}.bam -o $WORK_PATH/10_bigwing/sorted_${line}.bam -@ 10 2> $WORK_PATH/10_bigwing/${i}.sort.log
	 macs2 pileup --extsize 10 \
		-i $WORK_PATH/10_bigwing/sorted_${line}.bam \
		-o $WORK_PATH/10_bigwing/${line}.bdg \
		2> $WORK_PATH/10_bigwing/${line}.bgconv.log
         sort -k1,1 -k2,2n $WORK_PATH/10_bigwing/${line}.bdg > $WORK_PATH/10_bigwing/sorted_${line}.bdg
         bedGraphToBigWig $WORK_PATH/10_bigwing/sorted_${line}.bdg $WORK_PATH/hg19.chrom.sizes $WORK_PATH/10_bigwing/${line}.bigwig
done
