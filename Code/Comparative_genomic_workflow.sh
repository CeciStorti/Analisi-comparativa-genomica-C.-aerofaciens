#!/bin/bash


###############################
# STEP 1: Genomes download    #
###############################

# The following code allows you to obtain the accession numbers of genomes available on Genome,
# excluding atypical ones and metagenomes.

datasets summary genome taxon "Collinsella aerofaciens" --assembly-source genbank --exclude-atypical --mag exclude --as-json-lines > collinsella_summary.jsonl
cat collinsella_summary.jsonl | jq '.accession' > accession_numbers.txt
sed 's/"//g' accession_numbers.txt > final_accession_numbers.txt 
rm accession_numbers.txt

# We use the accession numbers just obtained to download the genomes

ncbi-genome-download --assembly-accessions final_accession_numbers.txt bacteria --section genbank --formats fasta --flat-output --output-folder genomes

# The downloaded genomes are saved in the ‘genomes’ folder.

# Check: let’s count how many files the folder contains. 

file_count=$(find genomes -type f | wc -l )
echo "The number of files contained in the folder is: $file_count"


##################
# STEP 2: CHECKM #
##################

# We evaluate the quality of the genomes just downloaded using the CheckM tool.

checkm lineage_wf -t 50 -x fna genomes  results_checkm 

# We analyze the extended_qa_summary.tsv file obtained as output from running CheckM.
# The file contains the main statistics for each analyzed genome, reporting contamination and completeness values.
# We look for the assembly accessions of genomes with:

# - completeness < 90
# - contamination > 5
# Once identified, we will proceed with their removal.

# We look into the names of the columns 
awk -F'\t' '{for(i=1; i<=NF; i++) print i, $i}' extended_qa_summary.tsv | head -n 10 

# We select the columns related to completeness and contamination:
awk -F'\t' 'NR>1 && ($6 <90 || $7 >5) {print $1}' extended_qa_summary.tsv > accession_to_remove.txt
mkdir -p checkm_genomes
cp -r genomes checkm_genomes 

# Remove low-quality genomes:
while IFS= read -r file; do
    rm "checkm_genomes/$file.fna"
done < accession_to_remove.txt

# Check: let’s count how many files the folder contains. 
file_count=$(find checkm_genomes -type f | wc -l )
echo "The number of files contained in the folder is: $file_count"


#################
# STEP 3: PyANI #
#################

# We calculate the Average Nucleotide Identity (ANI) using the PyANI tool. In this case, two genomes belong to the same species if ANI ≥ 95%.
average_nucleotide_identity.py -i checkm_genomes -o out_pyani -m ANIm --write_excel -l logfile_pyani.log -g --gformat png

# Calculation of the Percentage of Conserved Proteins (POCP). Two genomes belong to the same species if their POCP is > 50%.
# For this, we run the dedicated script 
./POCP.sh


#########################################################################
# STEP 4: Obtain metadata and descriptive statistics of interest        #
#########################################################################

# After the cleaning phase, we can download the associated metadata. 
# We are interested in obtaining information about the host from which the sample was collected, 
# the geographic area, etc. 
# We query the NCBI Biosample database.

# We create a file with the accessions of the final genomes. 
# To do this, we need to remove the rows corresponding to "accession_to_remove.txt" 
# from the accessions we had previously obtained.

grep -v "GCA_039752835.1" final_accession_numbers.txt > accession.txt
grep -v "GCA_040915745.1" accession.txt > final_removed_accession.txt
rm accession.txt 

# Check: final number of genomes
wc -l final_removed_accession.txt 


# Metadata download

datasets summary genome accession $(< final_removed_accession.txt) --assembly-source genbank --exclude-atypical --mag exclude  --as-json-lines | dataformat tsv genome --fields accession,assminfo-biosample-accession,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value > biosample_info.tsv

# From long to wide:

python long_to_wide.py

grep -v "GCA_039752835.1" extended_qa_summary.tsv > tmp_stat.tsv
grep -v "GCA_040915745.1" tmp_stat.tsv > stats_checkm.tsv
rm tmp_stat.tsv
python descriptive_statistics.py
 

########################
# STEP 5: Annotation   #
########################

# Set the path to the autoprokka.py script
path_autoprokka="replace_with_the_path_to_autoprokka.py"

# Run genome annotation
python $path_autoprokka/autoprokka.py -i checkm_genomes -o out_prokka

# Once the annotation is finished, copy the ".gff" files into
# a second folder. These files will be used in the next step.
mkdir -p gff_file

find out_prokka -type f -name "*.gff" -exec cp {} gff_file/ \;



#####################################
# STEP 6: Pangenome Analysis        #
#####################################

# Roary parameters:
# -i 95 : minimum identity threshold to consider two genes as orthologs (95% here)
# -cd 99 : prevalence cut-off to consider a gene as core (99% ≤ strains ≤ 100%)

roary -e --mafft -f result_roary -p 40 -r -i 95 -cd 99 gff_file/*.gff

# Select genes that compose the core genome:
mkdir -p gff_for_roary
cp -r "result_roary/"* "gff_for_roary/"
cd gff_file
query_pan_genome -o query_pan_genome_result_core -a intersection *.gff
mv query_pan_genome_result_core result_roary
find -maxdepth 1 -type f ! -name "*.gff" -exec rm {} +
cd ..

# Run Panaroo for core genome analysis:
panaroo -i gff_file/*.gff -o results_panaroo --clean-mode strict --remove-invalid-genes -a core --aligner mafft --core_threshold 0.99 \
        --merge_paralogs --refind_prop_match 0.5

# Run PPANGOLIN workflow for pangenome analysis
mkdir -p ppanggolin
cd ppanggolin
ppanggolin workflow --anno genomes.gbff.list --identity 0.95 --mode 2


##############################################
# STEP 7: Phylogenetic Analysis - PhyloPhlAn #
##############################################

# Create folders for phylogeny and logs
mkdir -p filogenesi/logs

# Generate a marker database for the bacterium under study
phylophlan_setup_database -g s__Collinsella_aerofaciens -o . --verbose 2>&1 | tee filogenesi/logs/phylophlan_setup_database.log

# Set up the configuration file
phylophlan_write_config_file -o filogenesi/file_config.cfg -d a --force_nucleotides --db_aa diamond --map_aa diamond --map_dna diamond --msa mafft --trim trimal --tree1 fasttree

# Run phylogenetic tree construction
phylophlan -i checkm_genomes -o filogenesi/output_phylogeny -d s__Collinsella_aerofaciens -t a -f filogenesi/file_config.cfg --nproc 20 --diversity low --force_nucleotides --accurate --verbose 2>&1 | tee filogenesi/logs/phylophlan__s__Collinsella_aerofaciens.log

# Build tree with bootstrap support using IQ-TREE
iqtree -s genomi_checkm_concatenated.aln -m MFP -bb 1000 -T AUTO | tee filogenesi/logs/iqtree.log


#######################
# STEP 8: 16S Analysis #
#######################

input_dir="checkm_genomes"
mkdir -p 16s tmp_16s
tmp_dir="tmp_16s"
out_dir="16s"

# Extract 16S rRNA sequences with barrnap
for file in "$input_dir"/*.fasta; do
    base_name=$(basename "$file" .fasta)
    barrnap --kingdom bac --outseq "$tmp_dir"/${base_name}.fasta $file
done

# Select 16S sequences from barrnap output
for file_16 in "$tmp_dir"/*.fasta; do
    genome_name=$(basename "$file_16" .fasta)
    grep -A 1 "16S" "$file_16" | grep -v "^--" > "$out_dir"/${genome_name}_16S.fasta
done

# Select the longest 16S sequence for each file
mkdir -p selected_16s tmp
for file in "16s"/*.fasta; do
    base=$(basename "$file" _16S.fasta)
    # Robust selection of longest sequence using seqtk / awk recommended
done

# Concatenate all selected sequences into a single FASTA file with standardized headers
output_file="16s/16_TOTALE.fasta"
for file in "selected_16s"/*.fasta; do
    base=$(basename "$file" .fasta | cut -d'_' -f1-2)
    sed "s/^>.*$/>$base/" "$file" >> "$output_file"
done

# Multiple sequence alignment using Clustal Omega
clustalo --infile 16s/16_TOTALE.fasta --guidetree-out clustalo-prova_auto.dnd \
         --distmat-out prova_clust_auto.matrix --full --outfmt clustal --auto --outfile aln_clustal_auto


########################
# STEP 9: ADH Analysis #
########################

# Download genome summary for Collinsella aerofaciens from RefSeq (exclude atypical and MAG)
datasets summary genome taxon "Collinsella aerofaciens" --assembly-source refseq --exclude-atypical --mag exclude --as-json-lines > collinsella_summary.jsonl

# Extract accession numbers
jq -r '.accession' collinsella_summary.jsonl > GCF.txt

# Create folder for protein sequences
mkdir -p PROTEINE

# Download protein FASTA files for the selected genomes
ncbi-genome-download --assembly-accessions GCF.txt bacteria --formats protein-fasta --flat-output --output-folder PROTEINE -v

# Merge all protein sequences into a single FASTA and create a BLAST database
cat PROTEINE/*.faa > sequenze_proteine_totali.fasta
makeblastdb -in sequenze_proteine_totali.fasta -dbtype prot -out db_blast_protein

# Run BLASTp using ADH protein as query
blastp -query ADH_protein.faa -db db_blast_protein -out risultati_blast.txt -outfmt 6

# Filter BLAST results using custom Python script
python filter_results_blast.py

# Select columns 2 and 3 from filtered results
cut -d $'\t' -f 2,3 "filtered_blast_results_ref_C_a_4.tsv" > tmp_2.txt

# Convert filtered results to CSV using custom Python script
python ADH_to_csv.py

####################
# STEP 10: CAZYmes #
####################


input_dir="genomi_checkm"

mkdir -p cazyme
output_dir="cazyme"


for file in "$input_dir"/*.fasta; do
    filename=$(basename "$file" .fasta)
    run_dbcan "$file" prok --out_dir "$output_dir/${filename}" --verbose --dbcan_thread 4
done

mkdir -p overview

for dir in cazyme/*; do
    cp "$dir/overview.txt" overview/annotations_$(basename $dir).txt
done

python cazyme_dendrogramma.py 


############################
# STEP 11: Gene resistance #
############################

mkdir abricate_results
abricate --db argannot --csv genomi_checkm/*.fasta > abricate_results/results_argannot.csv
head -n 1 abricate_results/results_argannot.csv > abricate_results/filtered_argannot.csv
awk -F, '{if (NR >1 && $10 > 80 && $11 > 75) print }' abricate_results/results_argannot.csv >> abricate_results/filtered_argannot.csv
abricate --summary --csv abricate_results/filtered_argannot.csv > abricate_results/summary_results_argannot.csv


abricate --db ncbi --csv genomi_checkm/*.fasta > abricate_results/results_ncbi.csv
head -n 1 abricate_results/results_ncbi.csv > abricate_results/filtered_ncbi.csv
awk -F, '{if (NR >1 && $10 > 80 && $11 > 75) print }' abricate_results/results_ncbi.csv >> abricate_results/filtered_ncbi.csv
abricate --summary --csv abricate_results/filtered_ncbi.csv > abricate_results/summary_results_ncbi.csv

abricate --db card --csv genomi_checkm/*.fasta > abricate_results/results_card.csv
head -n 1 abricate_results/results_card.csv > abricate_results/filtered_card.csv
awk -F, '{if (NR >1 && $10 > 80 && $11 > 75) print }' abricate_results/results_card.csv >> abricate_results/filtered_card.csv
abricate --summary --csv abricate_results/filtered_card.csv > abricate_results/summary_results_card.csv


abricate --db resfinder --csv genomi_checkm/*.fasta > abricate_results/results_resfinder.csv
head -n 1 abricate_results/results_resfinder.csv > abricate_results/filtered_resfinder.csv
awk -F, '{if (NR >1 && $10 > 80 && $11 > 75) print }' abricate_results/results_resfinder.csv >> abricate_results/filtered_resfinder.csv
abricate --summary --csv abricate_results/filtered_resfinder.csv > abricate_results/summary_results_resfinder.csv

# FATTORI DI VIRULENZA

abricate --db vfdb --csv genomi_checkm/*.fasta > abricate_results/results_vfdb.csv
abricate --summary --csv abricate_results/results_vfdb.csv > abricate_results/summary_vfdb.csv

# PLASMIDI

abricate --db plasmidfinder --csv genomi_checkm/*.fasta > abricate_results/results_plasmid.csv
abricate --summary --csv abricate_results/results_plasmid.csv > abricate_results/summary_plasmid.csv



