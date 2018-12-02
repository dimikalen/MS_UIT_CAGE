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

cd $project_home

#####################################Unzip_fastq#########################
	
	gunzip *.gz
	
#####################################Rename#########################	
for i in *.fastq; do mv $i $(echo $i | sed 's/\_/\./g'); done

#####################################Store the first 3 column, before RX #########################

ls *\.R1\.*\.fastq |cut -d . -f 1,2,3 > list.txt

#####################################MERGE_QUALITY_FILTERING#########################	

	while read sample; do
  		bbmerge.sh in1=$sample.R1.$RX in2=$sample.R2.$RX out=$sample.merged.fq
  		vsearch -fastq_filter $sample.merged.fq -fastaout $sample.filtered.fa -fastq_maxee 1 -threads $cpu
	done < list.txt

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

	mkdir $name.results
	mv  $name.taxonomy.txt $name.result/
	mv  $name.taxonomy.biom $name.result/
	mv  $name.otus$id.fa $name.result/

exit