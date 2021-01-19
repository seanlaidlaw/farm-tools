#!/usr/bin/env bash
# first arg should be input folder, second arg output folder

# 37
genome_dir="$HOME/ref/hg19_cDNA_genome.fa"
# 38
#genome_dir="/lustre/scratch119/realdata/mdt1/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh38_ERCC92/genome.fa"

bamout="$2"



cram_list="$(find "$1" -name '*.cram' |sort -u)"
mkdir -p "$bamout/logs"

module load samtools
for cram_file in $cram_list;do
	base_extensionless="$(basename "$cram_file" | sed 's/.cram//g')"

	bsub -e "$bamout/logs/bam_$base_extensionless.e" -o "$bamout/logs/bam_$base_extensionless.o" \
		-J "bam_$base_extensionless" \
		-R'select[mem>1000] rusage[mem=1000]' -M1000 -n 2 -R'span[hosts=1]' \
		-q small \
		samtools view -O bam \
			-@ 2 -T "$genome_dir" -o "$bamout/$base_extensionless.bam" \
			"$cram_file"
done



