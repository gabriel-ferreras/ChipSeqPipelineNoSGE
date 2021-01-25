#! /bin/bash
SAMPLE_DIR=$1
i=$2
NUM_SAMPLES=$3
INS_DIR=$4
ANALYSIS=$5
BROAD=$6
PAIRED=$7
UPSTREAM=$8
DOWNSTREAM=$9
MOTIFLENGTH=${10}
MOTIFSIZE=${11}
NUM_EXP=${12}
EXP_DESIGN=${13}
CHR=${14}

echo ""
echo "=========================="
echo "|   PROCESSING CHIP $i   |"
echo "=========================="
echo ""

cd $SAMPLE_DIR

## Sample quality control and read mapping to reference genome
if [ $PAIRED -eq 0 ]
then
        fastqc chip_$i.fastq.gz
        bowtie2 -x ../../genome/index -U chip_$i.fastq.gz -S chip_$i.sam
else
	fastqc chip_${i}_1.fastq.gz
        fastqc chip_${i}_2.fastq.gz
        bowtie2 -x ../../genome/index -1 chip_${i}_1.fastq.gz -2 chip_${i}_2.fastq.gz -S chip_$i.sam
fi

## Generating sorted bam file
samtools sort -o chip_$i.bam chip_$i.sam
rm chip_$i.sam
rm *.fastq.gz
samtools index chip_$i.bam

echo ""
echo "   Chip $i processing DONE!!"
echo ""
