#!/bin/bash
#SBATCH -J practice1
#SBATCH -p mrcs
#SBATCH -w Zeta
#SBATCH -n 14
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#
mkdir YC6
cd YC6

name=YC6
HomeDIR=/home/arajo0707/Rat/YC6
FASTQDIR=/home/arajo0707/Rat/link/YC6
#
echo $name

STAR \
	--genomeDir /home/data/star_index/Rnor_6.0/  \
	--readFilesIn /home/arajo0707/Rat/link/YC6_sol_1.fastq.gz /home/arajo0707/Rat/link/YC6_sol_2.fastq.gz \
	--runThreadN 14 \
	--outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 1 \
        --genomeLoad NoSharedMemory \
        --readFilesCommand zcat \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --sjdbOverhang 100 \
        --outSAMstrandField intronMotif \
        --outSAMtype None \
        --outSAMmode None \
        --outFileNamePrefix $name

STAR \
        --runMode genomeGenerate \
        --genomeDir $HomeDIR \
        --genomeFastaFiles /home/data/star_index/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
        --sjdbOverhang 100 \
        --runThreadN 14 \
        --sjdbFileChrStartEnd $name\SJ.out.tab

STAR \
	--genomeDir $HomeDIR \
	--readFilesIn /home/arajo0707/Rat/link/YC6_sol_1.fastq.gz /home/arajo0707/Rat/link/YC6_sol_2.fastq.gz \
	--runThreadN 14 \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 10 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 1000000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 1 \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM 0 \
        --readFilesCommand zcat \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --sjdbOverhang 100 \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI NM MD AS XS \
        --outSAMunmapped Within \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMheaderHD @HD VN:1.4 \
        --outSAMattrRGline ID:$name LB:$name SM:$name PL:Illumina \
        --outFileNamePrefix $name

htseq-count -m intersection-nonempty -s no -f bam $name\Aligned.sortedByCoord.out.bam /home/data/star_index/Rattus_norvegicus.Rnor_6.0.103.gtf > $name\.count

#rm $name\Aligned.sortedByCoord.out.bam
/home/tools/gatk-4.1.8.1/gatk CollectInsertSizeMetrics -H $name\_histogram.pdf -I $name\Aligned.sortedByCoord.out.bam -O $name\_size.txt
cp $name\.name /home/arajo0707/Rat/count
cp /data/MRC3_home/arajo0707/Rat/YC6/YC6_size.txt /data/MRC3_home/arajo0707/Rat/fpkm
