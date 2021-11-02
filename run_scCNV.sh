#!/usr/bin/env bash

# Run specific paths
while getopts ":B:O:S:" opt; do
  case $opt in
    B) inputbams="$OPTARG"
    ;;
    O) outputfolder="$OPTARG"
    ;;
    S) speciesarg="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# make sure minimum arguments are provided
if [ "$inputbams" == "" ]; then
	echo "ERROR: missing input bams '-B' argument"
	exit 1
fi
if [ "$outputfolder" == "" ]; then
	echo "ERROR: missing output directory '-O' argument"
	exit 1
fi




if [ "$speciesarg" == "human" ]; then
	# Farm 5 Paths for human
	mappableBinsDir="/lustre/scratch117/casm/team176/sl31/scCNV_mappability/mappability_files_organized/human/GRCh37d5/1Mb"

	chr_list='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X'
	codedir="/lustre/scratch119/realdata/mdt1/team176/sl31/sccnv_seans_upgrades/scCNV/"
	badbins_path="/lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/scCNV/bad_bins.bed.gz"
	cytoband_path="/lustre/scratch119/realdata/mdt1/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCH37d5/brass/cytoband_hg19.txt"

elif [ "$speciesarg" == "mouse" ]; then
	# Farm 5 Paths for mouse
	mappableBinsDir="/lustre/scratch117/casm/team176/sl31/scCNV_mappability/mappability_files_organized/mouse/GRCm38/100kb"

	chr_list='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19'
	cytoband_path="/nfs/cancer_ref02/Mus_musculus/GRCm38/brass/cytoband.txt"
	badbins_path="/nfs/cancer_ref02/Mus_musculus/GRCm38/scCNV/bad_bins.bed.gz"
else
	echo "Unsupported species, please manually add to wrapper script"
	exit 1
fi

sccnv_singularity_image="/software/CASM/singularity/sccnv/sccnv_v2.0.0.sif"

if [[ ! $(tree "$inputbams" | grep ".bai$") ]]; then
	cd "$inputbams"

	module load samtools

	for unsorted_bam in $(find -type f -name '*.mem.bam');do
		bsub -e samtools_sort_index_${unsorted_bam%.mem.*}.e -o samtools_sort_index_${unsorted_bam%.mem.*}.o -R'select[mem>1000] rusage[mem=1000]' -M1000 \
			"samtools sort $unsorted_bam > ${unsorted_bam%.mem.*}.sorted.bam && samtools index ${unsorted_bam%.mem.*}.sorted.bam && rm $unsorted_bam"
	done
	cd -
	echo "> No index files found for bam files in input directory, samtools sort and index jobs are launched. When jobs finish rerun this command to run scCNV"
	exit 1
else
	# make output folder if doesnt exist
	if [ ! -d "$outputfolder" ]; then
		mkdir -p "$outputfolder"
	fi

	# copy bad_bins_file to expected location if not already
	if [ ! -f "$outputfolder/bad_bins_file.gz" ]; then
		cp -vn "$badbins_path" "$outputfolder/bad_bins_file.gz"
	fi

	if [ ! -f "$outputfolder/cytoband.txt.gz" ]; then
		cp -vn "$cytoband_path" "$outputfolder/cytoband.txt"
		cd "$outputfolder"
		gzip "cytoband.txt"
		cd -
	fi

	# Create sample csv with filename and local file path by writing the same filename without .bam and with .bam
	sampleList="$outputfolder/unfiltered_sample.list.csv"
	filteredSampleList="$outputfolder/filtered_sample.list.csv"
	if [ ! -f "$sampleList" ]; then
		echo 'sample,bam_file' > "$sampleList"
		ls "$inputbams" | grep 'bam$' | sed -E 's/(.*)\.bam/\1,\1.bam/g' >> "$sampleList"
	fi


	# remove samplename from sample list if the bam is empty (as scCNV aborts on first empty bam)
	if [ ! -f "$filteredSampleList" ]; then
		module load samtools
		echo "filtering sample list to $filteredSampleList"
		echo 'sample,bam_file' > "$filteredSampleList"
		for bam in $(cat $sampleList | cut -d ',' -f 2 | grep 'bam$');do
			bam_path="$inputbams/$bam"
			bam_contents_lines=$(samtools view $bam_path | head | wc -l)
			if (($bam_contents_lines != 0)); then
				cat $sampleList | grep ",$bam" >> $filteredSampleList
			else
				echo "empty bam: $bam_path"
			fi
		done
	fi
	threads=$(cat $filteredSampleList | wc -l)


	# generate manifest file from the provided mappability filename as mappability files have
	# standardised names
	# parse from mappability file data for manifest
	map_filename=$(ls "$mappableBinsDir" | sort -u | grep 'GCperc_INPUT.txt' | head -n 1 | xargs basename)
	map_species=$(echo "$map_filename" | cut -d '_' -f 2)
	map_assembly=$(echo "$map_filename" | cut -d '_' -f 3)
	map_binsize=$(echo "$map_filename" | cut -d '_' -f 4)

	echo "> Generating mappability manifest for chrs: $chr_list"
	manifest_chrlist=$(echo "$chr_list" | tr ',' '\n' | sed "s/'//g" | sed 's/"//g' \
		| sed 's/^/"/g' | sed 's/$/",/g' | tr '\n' ' ' | sed 's/, $//g' | sed 's/ //g')
	echo "$manifest_chrlist"

	echo "{\"_version\":\"1.0.0\",\"metadata\":{\"species\":\"$map_species\",\"assembly\":\"$map_assembly\",\"read_length\":126,\"chromosomes\":[$manifest_chrlist],\"bin_size\":$map_binsize},\"files\":{\"bed\":\"Combined_${map_species}_${map_assembly}_${map_binsize}_126bases_mappable_bins.txt\",\"input\":\"Combined_${map_species}_${map_assembly}_${map_binsize}_126bases_mappable_bins_GCperc_INPUT.txt\"}}" > "$outputfolder/manifest.json"


	sampleList=$filteredSampleList
	threads=24
	echo "> running scCNV for bin size: $map_binsize"
	bsub -e "$outputfolder/scCNV.e" -o "$outputfolder/scCNV.o" -G team176 -R'select[mem>8000] rusage[mem=8000]' -M8000 -n $threads -R'span[hosts=1]' \
	"singularity exec \
	--cleanenv \
	-B ${inputbams}:/data:ro \
	-B ${mappableBinsDir}:/mapBins:ro \
	-B ${sampleList}:/sample.list.txt:ro \
	-B ${badbinsdir}:/output/bad_bins_file.gz:ro \
	-B ${outputfolder}:/output \
	$sccnv_singularity_image scCNV.pl \
		--sample_names_file /sample.list.txt \
		--threads $threads \
		--bamdir /data \
		--outdir /output \
		--bin $map_binsize \
		--readlength 126 \
		--gamma 15 \
		--species $map_species \
		--assembly $map_assembly \
		--cytoband /output/cytoband.txt.gz \
		--ploidy 2 \
		--threshold 0.6 \
		--bad_bins_file /output/bad_bins_file.gz \
		--mappable_bins_dir /mapBins \
		--mappable_bins_manifest /output/manifest.json \
		--chr_list $chr_list"

fi
