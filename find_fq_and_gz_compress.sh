#!/usr/bin/env bash


# first argument should be a folder in which to look for fastq files
input_folder="$1"


fastq_list="$(find "$input_folder" -name '*.fq' | sort -u)"


for fq_file in $fastq_list;do
	echo ""
	basename_file="$(basename "$fq_file")"
	output_folder="$(dirname "$fq_file")"
	mkdir -p "$output_folder/logs"

	cd "$output_folder"
	bsub -e "$output_folder/logs/gzip_$basename_file.e" -o "$output_folder/logs/gzip_$basename_file.o" \
		-R'select[mem>6000] rusage[mem=6000]' -M6000 -n 6 -R'span[hosts=1]' \
		pigz --best --keep --processes 6 "$fq_file"
	cd -
done



