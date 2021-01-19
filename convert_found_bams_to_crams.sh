#!/usr/bin/env bash

# 37
#genome_dir="$HOME/ref/hg19_cDNA_genome.fa"
# 38
genome_dir="/lustre/scratch119/realdata/mdt1/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh38_ERCC92/genome.fa"

bamout="38_cram_output"



bam_list="$(find "3_reAligned_to_38" -name '*.bam' |sort -u)"
mkdir -p "$bamout/logs"

module load samtools
for fq_file in $bam_list;do
	base_extensionless="$(basename "$fq_file" | sed 's/.bam//g')"

	bsub -e "$bamout/logs/cram_$base_extensionless.e" -o "$bamout/logs/cram_$base_extensionless.o" \
		-J "cram_$base_extensionless" \
		-R'select[mem>6000] rusage[mem=6000]' -M6000 -n 6 -R'span[hosts=1]' \
		samtools view -O cram,seqs_per_slice=100000,level=8,use_lzma \
			-@ 6 -T "$genome_dir" -o "$bamout/$base_extensionless.cram" \
			"$fq_file"
done



