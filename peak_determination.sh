#! /bin/bash
RES_DIR=$1
NUM_SAMPLES=$2
ANALYSIS=$3
BROAD=$4
INS_DIR=$5
UPSTREAM=$6
DOWNSTREAM=$7
MOTIFLENGTH=$8
MOTIFSIZE=$9
NUM_EXP=${10}
EXP_DESIGN=${11}
j=${12}
CHR=${13}

## Accessing results folder
cd $RES_DIR

## Creating sample result folder
mkdir sample_${j}_result
cd sample_${j}_result

## Peak determination
echo ""
echo "  --Hey I'm the peak_determination.sh script, let me process you sample $j data real quick!"
echo ""

if [ $BROAD -eq 0 ]
then
	echo ""
	echo "  --Determining NARROW peaks for sample $j"
	echo ""
        macs2 callpeak -t ../../samples/chip_$j/chip_$j.bam -c ../../samples/control_$j/control_$j.bam -f BAM --outdir . -n ${ANALYSIS}_sample_$j

else
	echo ""
	echo "  --Determining BROAD peaks for sample $j"
	echo ""
        macs2 callpeak -t ../../samples/chip_$j/chip_$j.bam -c ../../samples/control_$j/control_$j.bam -f BAM --outdir . -n ${ANALYSIS}_sample_$j --broad

fi

## Initiating R script for peak annotation

if [ $BROAD -eq 0 ]
then
	echo ""
	echo "  --Annotating NARROW peaks for sample $j"
	echo ""
	Rscript ${INS_DIR}/ChipSeqPipelineNoSGE/target_genes.R ${ANALYSIS}_sample_${j}_peaks.narrowPeak ${ANALYSIS}_sample_${j}_summits.bed $UPSTREAM $DOWNSTREAM ${ANALYSIS}_sample_${j}_peaks_targetgenes.txt ${ANALYSIS}_sample_${j}_summits_targetgenes.txt

else
	echo ""
	echo "  --Annotating BROAD peaks for sample $j"
	echo ""
	Rscript ${INS_DIR}/ChipSeqPipelineNoSGE/target_genes_broad.R ${ANALYSIS}_sample_${j}_peaks.broadPeak $UPSTREAM $DOWNSTREAM ${ANALYSIS}_sample_${j}_peaks_targetgenes.txt
fi

## Homer motif finding, only for narrow peaks.
if [ $BROAD -eq 0 ]
then
	echo ""
	echo "  --Finding motives for sample $j. OK I lied this isn't gonna be that quick..."
	echo ""
	mkdir motifs_sample_${j}
	cd motifs_sample_${j}
	findMotifsGenome.pl ../${ANALYSIS}_sample_${j}_summits.bed tair10 . -len $MOTIFLENGTH -size $MOTIFSIZE
	echo ""
	echo "  --It wasn't that bad, was it?"
	echo ""
fi
echo ""
echo "  --Peak determination for sample $j done!!"
echo ""
cd ..
cd ..

