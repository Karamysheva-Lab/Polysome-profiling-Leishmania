#!/usr/bin/env zsh


####################################################################
######################    Trimming with fastp      #################
####################################################################

#ls fastq_files | cut -d_ -f 1 | uniq > samples.txt

# Make a directory to saved trimmed files and reports
#mkdir trimmed
#mkdir trimmed/reports

#while read file;
#do
	#if [ ! -f trimmed/$file\_R1.fq.gz ]; then
		#fastp -i fastq_files/$file\_*R1_001.fastq.gz -I fastq_files/$file\_*R2_001.fastq.gz -o trimmed/$file\_R1.fq.gz -O trimmed/$file\_R2.fq.gz -c --detect_adapter_for_pe -w 16 -q 30 -h trimmed/reports/$file.html -j trimmed/reports/$file.json -R $file\_report;
	#fi
#done < samples.txt

#gunzip trimmed/*.fq.gz

####################################################################
######################    Alignment with STAR     #################
####################################################################


#wget https://tritrypdb.org/common/downloads/release-67/LmajorFriedlin/fasta/data/TriTrypDB-67_LmajorFriedlin_Genome.fasta

#wget https://tritrypdb.org/common/downloads/release-67/LmajorFriedlin/gff/data/TriTrypDB-67_LmajorFriedlin.gff

export PATH=/Volumes/Backup_ZK/Bio_tools/STAR-2.7.11b/bin/MacOSX_x86_64/:$PATH

#mkdir LmajorFriedlin_67_STARindex

#STAR --runThreadN 16 --runMode genomeGenerate --genomeDir LmajorFriedlin_67_STARindex --genomeFastaFiles TriTrypDB-67_LmajorFriedlin_Genome.fasta --sjdbGTFfile TriTrypDB-67_LmajorFriedlin.gff --sjdbOverhang 143 --genomeSAindexNbases 11

#mkdir BAM_files

# Loop to align all trimmed fastq files
while read line; do
	if [ ! -d BAM_files/$line ]; then
		STAR --runThreadN 40 --runMode alignReads --outSAMtype BAM SortedByCoordinate --genomeDir LmajorFriedlin_67_STARindex --readFilesIn trimmed/$line\_R1.fq trimmed/$line\_R2.fq  --outFileNamePrefix BAM_files/$line/$line\_ --limitBAMsortRAM 1851428885;
	fi
done < samples.txt


####################################################################
######################      MarkDuplicates         #################
####################################################################

export PATH=/Volumes/Backup_ZK/Bio_tools/jdk-22.0.1.jdk/Contents/Home/bin:$PATH

mkdir MarkDup

while read file;
do
	java -jar /Volumes/Backup_ZK/Bio_tools/picard/picard.jar MarkDuplicates -I BAM_files/$file/$file\_*.bam -O MarkDup/$file\_MarkDup.bam -M MarkDup/$file\_metrics.txt --REMOVE_DUPLICATES
	echo $file MarkDuplicates done
done < samples.txt

####################################################################
######################      FeatureCounts         ##################
####################################################################

## Call featureCounts and count
export PATH=/Volumes/Backup_ZK/Bio_tools/subread-2.0.6-macOS-x86_64/bin/:$PATH

featureCounts -p --countReadPairs -T 10 -t exon -g gene_id -s 0 -a TriTrypDB-67_LmajorFriedlin.gff -F GTF -M -o raw_counts.txt  MarkDup/*.bam




