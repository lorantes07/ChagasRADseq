# Chagas Disease RAD-seq Bioinformatics pipeline
Author: Lucia Orantes (lorantes@vermont.gov)
Last updated: 3/6/2018
### Available data: [http://www.ncbi.nlm.nih.gov/bioproject/389426](http://www.ncbi.nlm.nih.gov/bioproject/389426)

Full bioinformatics pipeline to extract mixed-DNA sequences from RAD sequences    

![](https://cloud.githubusercontent.com/assets/24459510/24257107/509ad8d8-0fc0-11e7-8303-0c87b4364b18.jpg)

### Trimming + filtering 
```
> module load FASTX tools
>for X in *.fq; do $FASTX -Q33 -q 10 -p 100  -i $X -o ${X%.*}\_filtered.fq; mv ${X%.*}\_filtered.fq ../filtered; done
```
### Alignment to  parasite reference genomes
```
> load bowtie 
```
*T. cruzi* genomes indexes
```
NCBI accession #
 - AAHK01 
 - ADWP02 
 - AHKC01 
 - ANOX01 
 - AODP01 
 - AQHO01
```
Run bowtie to map trimmed reads to  \*T. cruzi\* indexes in continous loop
```
> cp ../filtered/*.fq .
> cp -r ../indexes .
> mkdir unmapped\_files sam\_files Tdim\_files mapped\_files
> for X in *.fq; do echo $X; bowtie AAHK01 $X -S ${X%.*}\_AAHK01.sam --un ${X%.*}\_unmapped1.fq; mv $X ./mapped\_files; bowtie ADWP02 ${X%.*}\_unmapped1.fq -S ${X%.*}\_ADWP02.sam --un ${X%.*}\_unmapped2.fq; bowtie AHKC01 ${X%.*}\_unmapped2.fq -S ${X%.*}\_AHKC01.sam --un ${X%.*}\_unmapped3.fq; bowtie ANOX01 ${X%.*}\_unmapped3.fq -S ${X%.*}\_ANOX01.sam --un ${X%.*}\_unmapped4.fq; bowtie AODP01 ${X%.*}\_unmapped4.fq -S ${X%.*}\_AODP01.sam --un ${X%.*}\_unmapped5.fq; bowtie AQHO01 ${X%.*}\_unmapped5.fq -S ${X%.*}\_AQHO01.sam --un ${X%.*}\_unmapped6.fq; bowtie Other\_Tcruzi ${X%.*}\_unmapped6.fq -S ${X%.*}\_Other-Tcruzi.sam --un ${X%.*}\_unmapped7.fq; bowtie TcIV ${X%.*}\_unmapped7.fq -S ${X%.*}\_TcVI.sam --un ${X%.*}\_unmappedforTdim.fq; mv *.sam ./sam\_files; mv \*\_unmappedfor\* ./Tdim\_files; mv \*\_unmapped\*  ./unmapped\_files; done
> grep "XA" Dorn\_Triatomine\_S210Ab_* > S210Ab\_Tcruzi.sam #run this command for all samples infected with T. cruzi to create individual files containing only mapped sequences. 
>grep -c "XA" *\_Tcruzi.sam #counts the total number of reads for each concatenated file of the infected samples. C\[lorantes@h1 scripts\_020116\]$ qsub ref\_map\_TcruziFloragenex.sh\_
> ls | cut -d "\_" -f 3 | uniq> tcruziList.txt #create a file with the names of the unique samples, stripped of any repetitive part of the name.
```
### Denovo leg catalog assembly
```
> load stacks
```
*T. dimidiata* catalog samples 
```
-TP901sequence1filtered_unmappedforTdim.fq
-TP-904sequence1filtered_unmappedforTdim.fq
-TP-903sequence1filtered_unmappedforTdim.fq
-TP-930sequence1filtered_unmappedforTdim.fq
-TP-939sequence1filtered_unmappedforTdim.fq
-CHJ462Leg_sequence1filtered_unmappedforTdim.fq
-Codigo-*sequence1filtered_unmappedforTdim.fq
-S431Leg_sequence1filtered_unmappedforTdim.fq
```
Run a denovo assembly to create a *T. dimidiata* catalog
```
> denovo_map.pl -T 15 -m 2 -M 3 -N 5 -n 3 -o ./ -S -i 1 -b 1 -t -s ./TP901_sequence_1_filtered_unmappedforTdim.fq -s ./TP-903_sequence_1_filtered_unmappedforTdim.fq -s ./TP-904_sequence_1_filtered_unmappedforTdim.fq -s ./TP-930_sequence_1_filtered_unmappedforTdim.fq -s ./Codigo-0124_sequence_1_filtered_unmappedforTdim.fq -s ./TP-939_sequence_1_filtered_unmappedforTdim.fq -s ./Codigo-0074_sequence_1_filtered_unmappedforTdim.fq -s ./Codigo-0188_sequence_1_filtered_unmappedforTdim.fq -s ./Codigo-0225_sequence_1_filtered_unmappedforTdim.fq -s ./Codigo-0275_sequence_1_filtered_unmappedforTdim.fq -s ./S431Leg_sequence_1_filtered_unmappedforTdim.fq -s ./CHJ462Leg_sequence_1_filtered_unmappedforTdim.fq
```
Convert leg catalog to bowtie index
```
> cut -f 3,10 batch_1.catalog.tags.tsv > catalog.tab #extract the two important columns from the cat catalog tag file. 
> awk '{print ">"$1"\n"$2}' catalog.tab > catalog.fa #convert from tab to fasta 
> run bowtie-buil #convert leg catalog to index.
```
Align vector catalog
```
> cp ../unmapped_files/Dorn_Triatomine_*_sequence_1_unmappedforTdim.fq .
> cp -r ../indexes/ .
> mkdir unmapped_files sam_files mapped_files
> for X in *.fq; do bowtie tdimFloragenex $X -S ${X%.*}_tdim.sam --un ${X%.*}_unmappedforBLAST.fq; mv $X ./mapped_files; mv *.sam ./sam_files; mv *_unmappedforBLAST* ./unmapped_files; done
```
Map to *T. cruzi* reference catalog
```
> cd ../refMapTcruzi
> cp ../bowtieTcruzi/sam_files/*_Tcruzi.sam .
> mkdir outputTcruzi 
ref_map.pl -s A1215Ab_Tcruzi.sam -s A9252_Tcruzi.sam -s A9253_Tcruzi.sam -s Abdomen5_Tcruzi.sam -s B46_Tcruzi.sam -s BS230Ab_Tcruzi.sam -s CHJ462Ab_Tcruzi.sam -s CHJ52Ab_Tcruzi.sam -s S210Ab_Tcruzi.sam -o ./outputTcruzi -n 6 -m 1 -T 8 -S -i 1 -b 1
```
Map to *T. dimidiata* reference catalog
```
> cd ../refMapTdim
> cp ../bowtieTdim/sam_files/*_tdim.sam . 
> mkdir ./outputTdim
>ref_map.pl -s A1215Ab_sequence_1_filtered_unmappedforTdim_tdim.sam -s A5201_sequence_1_filtered_unmappedforTdim_tdim.sam -s A9252_sequence_1_filtered_unmappedforTdim_tdim.sam -s A9253_sequence_1_filtered_unmappedforTdim_tdim.sam -s Abdomen5_sequence_1_filtered_unmappedforTdim_tdim.sam -s B46_sequence_1_filtered_unmappedforTdim_tdim.sam -s BS230Ab_sequence_1_filtered_unmappedforTdim_tdim.sam -s C1504_sequence_1_filtered_unmappedforTdim_tdim.sam -s CHJ14_sequence_1_filtered_unmappedforTdim_tdim.sam -s Dorn_Triatomine_CHJ21_sequence_1_filtered_unmappedforTdim_tdim.sam -s CHJ3_sequence_1_filtered_unmappedforTdim_tdim.sam -s Dorn_Triatomine_CHJ462Ab_sequence_1_filtered_unmappedforTdim_tdim.sam -s CHJ462Leg_sequence_1_filtered_unmappedforTdim_tdim.sam -s CHJ506Ab_sequence_1_filtered_unmappedforTdim_tdim.sam -s CHJ52Ab_sequence_1_filtered_unmappedforTdim_tdim.sam -s Codigo-0074_sequence_1_filtered_unmappedforTdim_tdim.sam -s Codigo-0124_sequence_1_filtered_unmappedforTdim_tdim.sam -s Codigo-0188_sequence_1_filtered_unmappedforTdim_tdim.sam -s Codigo-0225_sequence_1_filtered_unmappedforTdim_tdim.sam -s Codigo-0275_sequence_1_filtered_unmappedforTdim_tdim.sam -s S134_sequence_1_filtered_unmappedforTdim_tdim.sam -s S210Ab_sequence_1_filtered_unmappedforTdim_tdim.sam -s S431Leg_sequence_1_filtered_unmappedforTdim_tdim.sam -s TP194_sequence_1_filtered_unmappedforTdim_tdim.sam -s TP838_sequence_1_filtered_unmappedforTdim_tdim.sam -s TP-900_sequence_1_filtered_unmappedforTdim_tdim.sam -s TP901_sequence_1_filtered_unmappedforTdim_tdim.sam -s TP-903_sequence_1_filtered_unmappedforTdim_tdim.sam -s TP-904_sequence_1_filtered_unmappedforTdim_tdim.sam -s TP-930_sequence_1_filtered_unmappedforTdim_tdim.sam -s TP-937_sequence_1_filtered_unmappedforTdim_tdim.sam -s TP-939_sequence_1_filtered_unmappedforTdim_tdim.sam -o ./outputTdim -n 6 -m 3 -T 8 -S -i 1 -b 1``` 
> wc -l *_sequence_1_filtered_unmappedforTdim_tdim.snps.tsv # count number of T. dimidiata SNPs per sample
```
### Loci + SNP retrival
Retrieve *T. cruz*  markers
```
> cd ../outputTcruzi
> genotypes -b 1 -P ./ -r 2 -m 1 --min_hom_seqs 2 -t GEN
> wc -l *_Tcruzi.snps.tsv #count number of T. cruzi SNPs per sample
> for X in *.tsv; do cut -f7 ${X%.*}_Tcruzi.matches.tsv | sort -n | tail -1 #cut column 7 from the xx.matches.tsv file, then sort it lowest to highest, then get the last number. This number is the highest depth of coverage for that sample
> for X in *.tsv; do cut -f7 ${X%.*}_Tcruzi.matches.tsv | sort -n | head -1 #cut column 7 from the xx.matches.tsv file, then sort it lowest to highest, then get the last number. This number is the lowest depth of coverage for that sample
> for X in *.tsv; do awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }' ${X%.*}_Tcruzi.matches.tsv # calculate de average depth of coverage for each sample by using the depth of coverage for each read described in the column N=7 of the XX.matches.tsv file
```
Retrieve *T. dimidiata*  markers
```
> cd .../outputTdim
> genotypes -b 1 -P ./ -r 3 -m 2 --min_hom_seqs 3 -t GEN
> wc -l *_Tdim.snps.tsv #count number of T. dimidiata SNPs per sample
> for X in *.tsv; do cut -f7 ${X%.*}_Tdim.matches.tsv | sort -n | tail -1 #cut column 7 from the xx.matches.tsv file, then sort it lowest to highest, then get the last number. This number is the highest depth of coverage for that sample
> for X in *.tsv; do cut -f7 ${X%.*}_Tdim.matches.tsv | sort -n | head -1 #cut column 7 from the xx.matches.tsv file, then sort it lowest to highest, then get the last number. This number is the lowest depth of coverage for that sample
>for X in *.tsv; do awk -v N=7 '{ sum += $N } END { if (NR > 0) print sum / NR }' ${X%.*}_Tdim.matches.tsv # calculate de average depth of coverage for each sample by using the depth of coverage for each read described in the column N=7 of the XX.matches.tsv file
