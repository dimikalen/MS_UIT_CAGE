cpu=8
name=PROJECTNAME
#OTU clustering % : 99; 97; 95; 90
OTU=97
#common name after R1/R2
RX=001.fastq
#Database for taxonomy
ref=/Users/dka018/Documents/DataBase.nosync/Mothur/Silva/132/silva.nr_v132.align
tax=/Users/dka018/Documents/DataBase.nosync/Mothur/Silva/132/silva.nr_v132.tax
#taxonomy database name withour the extension
taxname=nr_v132
#project_directory
project_home=/Users/dka018/Documents/Projets/CAGE/Amplicon_script.nosync
#uc2otulab path if not set up, otherwise just type uc2otutab.py
pyth=/Applications/python_scriptsuc2otutab.py
#min sequence length 
min=250
#max sequence length 
max=255
#forward
	#trim left
trimleftf=13
	#trunc length
trunclenf=224
#reverse
	#trim left
trimleftr=13
	#trunc length	    
trunclenr=155

cd $project_home

mkdir $name.OTUS
mkdir $name.ZOTUS

#####################################ZIP&move#########################
	gzip *.fastq
	cp *.gz $name.OTUS/

#####################################IMPORT in qiime2#########################
cp *.gz $name.ZOTUS
cd $name.ZOTUS
source activate qiime2

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $name.demux-paired-end.qza

rm *.fastq.gz

qiime demux summarize \
  --i-data $name.demux-paired-end.qza \
  --o-visualization $name.demux.qzv
  
#####################################DADA2#########################

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $name.demux-paired-end.qza \
  --p-trim-left-f $trimleftf \
  --p-trim-left-r $trimleftr \
  --p-trunc-len-f $trunclenf \
  --p-trunc-len-r $trunclenr \
  --p-n-threads $cpu \
  --p-max-ee 1 \
  --o-table $name.ZOTUS.table.qza \
  --o-representative-sequences $name.ZOTUS.rep-seqs.qza \
  --o-denoising-stats $name.ZOTUS.denoising-stats.qza

#####################################export_files#########################
qiime tools export \
  --input-path $name.ZOTUS.table.qza \
  --output-path $name.ZOTUS.table

qiime tools export \
  --input-path $name.ZOTUS.rep-seqs.qza \
  --output-path $name.ZOTUS.rep-seqs

qiime tools export \
  --input-path $name.ZOTUS.denoising-stats.qza \
  --output-path $name.ZOTUS.denoising-stats

mkdir $project_home/$name.ZOTUs.results

mv $name.ZOTUS.table/feature-table.biom $project_home//$name.ZOTUs.results/$name.ZOTUS.biom
mv $name.ZOTUS.rep-seqs/dna-sequences.fasta $project_home//$name.ZOTUs.results/$name.ZOTUS.fa
mv $name.ZOTUS.denoising-stats/stats.tsv $project_home//$name.ZOTUs.results/$name.DADA2.stats.tsv

#######################################ASSIGN_TAXONOMY###############################
cd $project_home/$name.ZOTUs.results

mothur "#classify.seqs(fasta=$name.ZOTUS.fa, reference=$ref, taxonomy=$tax, cutoff=75, processors=$cpu)"

biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp $name.*.wang.taxonomy -i $name.ZOTUS.biom -o $name.ZOTUS.taxonomy.biom

biom convert -i $name.ZOTUS.taxonomy.biom -o $name.ZOTUS.taxonomy.txt --to-tsv --header-key taxonomy --table-type "OTU table"

rm $name.*.wang.taxonomy
rm $name.*.wang.tax.summary
rm $name.ZOTUS.biom
rm mothur*

#####################################UNZIP_Rename#########################	
cd $project_home/$name.OTUS

gunzip *.gz
	
for i in *.fastq; do mv $i $(echo $i | sed 's/\_/\./g'); done

#####################################Store the first 3 column, before RX #########################

ls *\.R1\.*\.fastq |cut -d . -f 1,2,3 > list.txt

#####################################MERGE_QUALITY_FILTERING#########################	

	while read sample; do
  		bbmerge.sh in1=$sample.R1.$RX in2=$sample.R2.$RX out=$sample.merged.fq
  		vsearch -fastq_filter $sample.merged.fq -fastaout $sample.filtered.fa -fastq_maxee 1 -threads $cpu
	done < list.txt

rm *.fastq

#######################################SAMPLES_LABELING#############################

	
	ls *.filtered.fa |cut -d . -f 1,2,3 > $name.list.txt

	while read sample; do
  		sed 's/\t/;/g' $sample.filtered.fa > $sample.filtered.temp.fa
		sed 's/[[:blank:]]//g' $sample.filtered.temp.fa > $sample.filtered.clean.fa
  		sed "-es/^>\(.*\)/>\1;barcodelabel=$sample;/" < $sample.filtered.clean.fa > $sample.filtered.label.fasta
  
	done < $name.list.txt
	
	rm *.filtered.clean.fa
	rm *.filtered.temp.fa
	
#######################################SAMPLES_CONCATENATION########################

	cat *.filtered.label.fasta > $name.filtered.label.fasta

	mothur "#summary.seqs(fasta=$name.filtered.label.fasta)"


#######################################DEREPLICATION_AND _SIZE_SORTING##############

	vsearch -derep_fulllength $name.filtered.label.fasta -output $name.filtered.label.uniques.fasta -sizeout -threads $cpu -minseqlength $min -maxseqlength $max

	vsearch -sortbysize $name.filtered.label.uniques.fasta -output $name.filtered.label.uniques.sort.fasta -minsize 2


#######################################OTU_CLUSTERING_AND_CHIMERA_CHECK#############

	usearch -unoise3 $name.filtered.label.uniques.sort.fasta -zotus $name.otus100.fa

	vsearch -sortbylength $name.otus100.fa --output $name.otus100.sorted.fa

	for id in 99 97 95 90
do
   usearch -cluster_smallmem $name.otus100.sorted.fa -id 0.$id -centroids $id.fa
   sed 's/Zotu/Otu/g' $id.fa > $name.otus$id.fa
   rm $id.fa
   
done

#######################################ASSIGN_TAXONOMY###############################

	mothur "#classify.seqs(fasta=$name.otus$OTU.fa, reference=$ref, taxonomy=$tax, cutoff=75, processors=$cpu)"	

#######################################OTU_MAPPING##################################
	
	ide=$(awk "BEGIN {print ($OTU)/100}")
	
	usearch -otutab $name.filtered.label.fasta -otus $name.otus$OTU.fa -id $ide -otutabout $name.otumap -mapout $name.map.txt -threads $cpu 
	

#######################################OTU_TABLE_CONSTRUCTION########################
	
	biom convert --table-type="OTU table" -i $name.otumap -o $name.biom --to-json

	biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp $name.*.wang.taxonomy -i $name.biom -o $name.taxonomy.biom

	biom convert -i $name.taxonomy.biom -o $name.taxonomy.txt --to-tsv --header-key taxonomy --table-type "OTU table"

	mkdir $project_home/$name.OTUs.results
	mv  $name.taxonomy.txt $project_home/$name.OTUs.results
	mv  $name.taxonomy.biom $project_home/$name.OTUs.results
	mv  $name.otus$OTU.fa $project_home/$name.OTUs.results

exit