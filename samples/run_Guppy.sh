#!/bin/bash
#
# Basecalling using guppy on GPU
#
# If Guppy cannot find your GPU first check CUDA and Nvida drivers are ready
# 	$ nvidia-smi $ nvidia drivers installation
# 	$ nvcc -V # CUDA installation
#
# If both are working tell Guppy to find your GPU card, most likely at "channel" 0 adding
#	-x "cuda:0"
# 	to the guppy command OR
#	-x "auto"
#

threads=8

basecall="GPU" # [GPU|CPU]

source ../PARAMS.sh

################################################################################

echo ">>> GUPPY BASECALLING - RTA TRIMMING OFF <<<"

samples=(
    "hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.decap.REL5.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
	"hsa.dRNASeq.HeLa.polyA.REL5OH.long.1"
    "hsa.dRNASeq.HeLa.polyA.REL5.1"
	"hsa.dRNASeq.HeLa.polyA.PNK.REL5.1"
    "hsa.dRNASeq.HeLa.total.REL3.1"
    "hsa.dRNASeq.HeLa.total.REL3.2"
    "hsa.dRNASeq.HeLa.total.REL3.3"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.4"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.5"
    "hsa.dRNASeq.HeLa.total.REL5.long.REL3.6"
    "hsa.dRNASeq.SIRV.polyA.REL5.long.2"
)

for i in "${samples[@]}"; do
	echo "Sample: $i"
	sdir=$SAMPLE_DIR/$i
	mkdir -p $sdir/guppy

	# --trim_strategy none # Turn off trimming by guppy; we can also --trim_strategy rna or --trim_strategy dna to remove adapters specific to rna or dna
	# --reverse_sequence true # Already in the RNA basecalling model config

	if [ "$basecall" == "GPU" ]; then
		echo "Using GPU Guppy basecalling"

		guppy_basecaller \
			-x "auto" \
			--flowcell FLO-MIN106 \
		 	--kit SQK-RNA002 \
		 	--hp_correct on \
		 	--records_per_fastq 0 \
		 	--u_substitution off \
			--trim_strategy none \
		 	--input_path $sdir/raw/ \
			--save_path $sdir/guppy/ \
		 	--recursive \
		 	--compress_fastq \
		 	--fast5_out

	elif [ "$basecall" == "CPU" ]; then
		echo "Using CPU Guppy basecalling"

		guppy_basecaller \
			--num_callers $threads \
			--cpu_threads_per_caller 1 \
			--flowcell FLO-MIN106 \
		 	--kit SQK-RNA002 \
		 	--hp_correct on \
		 	--records_per_fastq 0 \
		 	--u_substitution off \
			--trim_strategy none \
		 	--input_path $sdir/raw/ \
			--save_path $sdir/guppy/ \
		 	--recursive \
		 	--compress_fastq \
		 	--fast5_out
	else
		echo "Please specify 'GPU' or 'CPU' in \$basecall variable"
		exit
	fi

	cat $(ls $sdir/guppy/fastq_runid_*.fastq.gz) > $sdir/guppy/reads.fastq.gz # Merge all fastq files from separate workers
	mkdir $sdir/fastq
	ln -s ../guppy/reads.fastq.gz $sdir/fastq/reads.1.fastq.gz # Link basecalled data
	ln -s $(ls -d guppy/workspace/*/fast5) $sdir/fast5 # Link basecalled fast5
done
