#!/usr/bin/env bash
run="$1"
lane="$2"
library_type="$3"

run_lane="$run"_"$lane"

# create the data directory where we will run all the following commands
projectdir="$(pwd)"
rundir="$projectdir/Data/iRODS_Data/1.1_Sorted_by_Run_Lane"
irodsdir="$rundir"
mkdir -p "$rundir" && cd "$rundir"

mkdir -p "$run_lane"

cd "$run_lane"
# donwload CRAM into folders organised by run_lane
if [[ ! $(tree | grep ".cram$") ]];then
	# this must be run interactive as kinit requires password
	echo "> Authenticating on IRODS"
	iinit

	echo "> Downloading from IRODS:"
	echo "imeta qu -z seq -d id_run = $run and lane = $lane and type = cram"
	imeta qu -z seq -d id_run = $run and lane = $lane and type = cram \
		| grep ":" | awk '{ print $2 }' | paste - - -d/ \
		| grep -v '_phix.cram' > "cram_download_list.txt"

	mkdir -p logs
	for dl_link in $(cat "cram_download_list.txt");do
		bname="$(basename $dl_link)"
		bsub -e "logs/dl_$bname.e" -o "logs/dl_$bname.o" -J "dl_irodslink_$bname" -R'select[mem>2000] rusage[mem=2000]' -M2000 -n 1 -R'span[hosts=1]' iget -K "$dl_link" "./"
	done


fi

# donwload imeta into folders organised by run_lane
bwait -w 'ended(dl_irodslink_*)' && echo "> Downloading has finished"
if [[ ! $(tree | grep ".imeta$") ]]; then
	chmod 664 ./* #set permissions for downloaded crams to non writable

	echo "> Downloading imeta files for downloaded crams"
	for cram in $(find . | grep "cram$" | sed 's/.*\///g' )
	do
		imeta ls -d /seq/$run/$cram > $cram.imeta
		if [[ $(cat "$cram.imeta" | grep 'phi') ]]; then
			rm -v "$cram.imeta"
			rm -v "$cram"
		fi

	done
fi
chmod -R 775 logs


fastqout="$projectdir/Data/2_extracted_fastqs/"
mkdir -p "$fastqout/logs"

echo "> converting crams to fastq"
module load samtools

for cram_file in $(ls *.cram);do
	# get filename without .cram extension
	extension_less="$(echo $cram_file | sed 's/.cram$//g')"

	# extract fastq information from each cram and compress to gz
	bsub -e "$fastqout/logs/cram2fastq_$cram_file.e" -o "$fastqout/logs/cram2fastq_$cram_file.o" -J "cram2fastq_$cram_file" \
		-R'select[mem>1000] rusage[mem=1000]' -M1000 -n 4 -R'span[hosts=1]' \
		samtools fastq -c 7 -@ 4 -1 "$fastqout/$extension_less".1.fq.gz -2 "$fastqout/$extension_less".2.fq.gz -0 /dev/null -s /dev/null -n  "$cram_file"

done

bwait -w 'ended(cram2fastq_*)' && echo "> Converting to fastq has finished"

bamout="$projectdir/Data/3_reAligned/"
mkdir -p "$bamout/logs"

#mkdir -p "$bamout/refs"
#cd "$bamout/refs"
#if [[ ! $(tree | grep ".fasta$") ]]; then
	#cp "$ref_genome" "./ref_genome.fasta"

	#if [ "$library_type" == "DNA" ]; then
		#bsub -e "index_ref.e" -o "index_ref.o" -J "index_ref" -R'select[mem>30000] rusage[mem=30000]' -M30000 -n 4 -R'span[hosts=1]' bwa index ref_genome.fasta
	#elif [ "$library_type" == "RNA" ]; then
		## run star to index fasta file
	#fi

	#bwait -w 'ended(index_ref)'
#fi
#ref_genome="$(pwd)/ref_genome.fasta"
#cd -


cd "$fastqout"
for fq_file in $(ls *.1.fq.gz);do
	extension_less=$(echo "$fq_file" | sed 's/.1.fq.gz//g')

	fastq_1_file="$extension_less.1.fq.gz"
	fastq_1_file="$(realpath $fastq_1_file)"
	zcat "$fastq_1_file" > "$bamout/$extension_less.1.fastq"
	fastq_2_file="$extension_less.2.fq.gz"
	fastq_2_file="$(realpath $fastq_2_file)"
	zcat "$fastq_2_file" > "$bamout/$extension_less.2.fastq"

	if [ "$library_type" == "DNA" ]; then
		module load bwa
		module load samtools

		bsub -e "$bamout/logs/bwa_$extension_less.e" -o "$bamout/logs/bwa_$extension_less.o" \
			-J "align_$fq_file" \
			-R'select[mem>30000] rusage[mem=30000]' -M30000 -n 6 -R'span[hosts=1]' \
			"bwa mem -t 3 \
				$ref_genome \
				$fastq_1_file \
				$fastq_2_file | samtools sort -@3 -l7 -o $bamout/$extension_less.sort.bam"
	elif [ "$library_type" == "RNA" ]; then

		bsub -e "$bamout/logs/star_$extension_less.e" -o "$bamout/logs/star_$extension_less.o" \
			-J "align_$extension_less" \
			-R'select[mem>50000] rusage[mem=50000]' -M50000 -n 10 -R'span[hosts=1]' \
			"/nfs/users/nfs_r/rr11/Tools/STAR-2.5.2a/bin/Linux_x86_64_static/STAR" --runThreadN 10 \
				--outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --limitBAMsortRAM 31532137230 \
				--outSAMtype BAM SortedByCoordinate --genomeDir /lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5_ERCC92/star/75/ \
				--readFilesIn "$bamout/$extension_less.1.fastq" "$bamout/$extension_less.2.fastq" \
				--outFileNamePrefix "$bamout/$extension_less.bam"
	fi
done
cd -

bwait -w 'done(align_*)' && echo "> Aligning has completed"

# set normal gtf
ref_gtf="/lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5_ERCC92/star/e75/ensembl.gtf"

# fix ERCC names not having a gene_name
cat "$ref_gtf"  | grep '^#' > "$projectdir/hg19_cDNA_ensembl__gene_name_grepd.gtf"
cat "$ref_gtf" | grep 'gene_name' | grep 'gene_id' >> "$projectdir/hg19_cDNA_ensembl__gene_name_grepd.gtf"
cat "$ref_gtf" | grep '^ERCC-' | grep -v 'gene_name' | perl -pe 's/gene_id(.*?;)/gene_id\1 gene_name\1/g' >> "$projectdir/hg19_cDNA_ensembl__gene_name_grepd.gtf"
#cat "$ref_gtf" | grep '^ERCC-' | grep -v 'gene_name' | perl -pe 's/gene_id(.*?;)/gene_id\\1; transcript_id\\1;gene_name\\1/g' >> "$projectdir/hg19"
ref_gtf="$projectdir/hg19_cDNA_ensembl__gene_name_grepd.gtf"
#echo "$ref_gtf"

countdir="$projectdir/Data/4_counts_matrix/"
mkdir -p "${countdir}"


if [ "$library_type" == "RNA" ]; then
	# add featurecounts to path
	#export PATH="/nfs/users/nfs_r/rr11/Tools/subread-1.5.1-source/bin/:$PATH"
	#export PATH="/software/team282/download/subread-2.0.0-Linux-x86_64/bin/:$PATH"

	cd "$bamout"
	bamlist="$(find . -name '*.bam')"

	bsub \
		-R"span[hosts=1]" \
		-o $countdir/featureCounts.o -e $countdir/featureCounts.e  \
		-J featureCounts \
		-n 12 \
		-q normal \
		-G team176 \
		-R 'select[mem>=5000] rusage[mem=5000]' -M5000 \
		/software/team282/download/subread-2.0.0-Linux-x86_64/bin/featureCounts -T 12 -Q 30 -p \
		-t exon -g gene_name -a "$ref_gtf" \
		-o "$countdir/FeatureCounts_matrix_cDNA_GRCh37d5_gene_name.tsv" $bamlist
	cd -
fi
