#!/bin/bash


###############################
# STEP 1: Download dei genomi #
###############################

# Il seguente codice permette di ottenere gli accession number dei genomi disponibili su Genome,
# esclusi gli atypical e i metagenomi

datasets summary genome taxon "Collinsella aerofaciens" --assembly-source genbank --exclude-atypical --mag exclude --as-json-lines > collinsella_summary.jsonl
cat collinsella_summary.jsonl | jq '.accession' > accession_numbers.txt
sed 's/"//g' accession_numbers.txt > final_accession_numbers.txt 
rm accession_numbers.txt

# Utilizziamo gli accession number appena ottenuti per scaricare i genomi

ncbi-genome-download --assembly-accessions final_accession_numbers.txt bacteria --section genbank --formats fasta --flat-output --output-folder genomi

# I genomi scaricati sono salvati nella cartella genomi

# check:Contiamo quanti file contiene la cartella 

file_count=$(find genomi -type f | wc -l )
echo "Il numero di file contenuti nella cartella è: $file_count"


##################
# STEP 2: CHECKM #
##################



# Valutiamo la qualità dei genomi appena scaricati, utilizzando il tool CheckM.

checkm lineage_wf -t 50 -x fna genomi  results_checkm 

# Analizziamo il file extended_qa_summary.tsv ottenuto in output eseguendo CheckM.
# Il file contiene le pricinipali statistiche per ciascun genoma analizzato, riporta il 
# valore di contaminazione e completezza. 
# Cerchiamo gli assembly accession dei genomi con:
# - completezza < 90
# - contaminazione > 5 
# Una volta individuati, procederemo alla loro rimozione.

awk -F'\t' '{for(i=1; i<=NF; i++) print i, $i}' extended_qa_summary.tsv | head -n 10 # trovo i numeri delle colonne
# Selezioniamo le colonne relative alla completezza e alla contaminazione:
awk -F'\t' 'NR>1 && ($6 <90 || $7 >5) {print $1}' extended_qa_summary.tsv > accession_da_rimuovere.txt
# copia i genomi contenuti nella prima cartella, nella seconda.
mkdir -p genomi_checkm
cp -r genomi genomi_checkm 
# Rimuoviamo i genomi non completi o contaminati:
while IFS= read -r file; do
    rm "genomi_checkm/$file.fna"
done < accession_da_rimuovere.txt

# Check: contiamo quanti file contiene la cartella 
file_count=$(find genomi_checkm -type f | wc -l )
echo "Il numero di file contenuti nella cartella è: $file_count"




#################
# STEP 3: PyANI #
#################

average_nucleotide_identity.py -i genomi_checkm -o out_pyani -m ANIm --write_excel -l logfile_pyani.log -g --gformat png




#########################################################################
# STEP 4: Otteniamo i metadati e le statistiche descrittive d'interesse #
#########################################################################



# Conclusa la fase di pulizia, possiamo scaricare i metadati associati. Siamo interessati ad ottenere informazioni 
# sull'host da cui è stato prelevato il campione, l'aerea geografica, ecc. 
# Interroghiamo il database Biosample di NCBI.

# Creiamo file con gli accession dei genomi finali, per farlo dobbiamo rimuovere le righe relative agli "accession_da_rimuovere.txt"
# dagli accession che avevamo ottenuto. 

grep -v "GCA_039752835.1" final_accession_numbers.txt > accession.txt
grep -v "GCA_040915745.1" accession.txt > final_removed_accession.txt
rm accession.txt 

# Check: 
wc -l final_removed_accession.txt # contiene 272 righe 


# Download dei metadati 

datasets summary genome accession $(< final_removed_accession.txt) --assembly-source genbank --exclude-atypical --mag exclude  --as-json-lines | dataformat tsv genome --fields accession,assminfo-biosample-accession,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value > biosample_info.tsv

# Trasformazione dei dati da long a wide:

python long_to_wide.py

grep -v "GCA_039752835.1" extended_qa_summary.tsv > tmp_stat.tsv
grep -v "GCA_040915745.1" tmp_stat.tsv > statistiche_checkm.tsv
rm tmp_stat.tsv
python statistiche_descrittive.py
 

########################
# STEP 5: Annotazione  #
########################


path_autoprokka="sostituire con il percorso del file autoprokka.py"

python $path_autoprokka/autoprokka.py -i genomi_checkm -o out_prokka

# Un volta conclusa l'annotazione, copiamo i file con estensione ".gff" in 
# una seconda cartella. Questi file verranno utilizzati nello step successivo.

mkdir -p gff_file

find out_prokka -type f -name "*.gff" -exec cp {} gff_file \;



#####################################
# STEP 6: Pangenome analysis #
#####################################



# I valori con cui viene chiamato roary sono:
#-i 95: Soglia di identità minima per considerare due geni come ortologhi (95% in questo caso).
#-cd 99: Cut-off di prevalenza per considerare un gene come core (99% ≤ ceppi ≤ 100%)

roary -e --mafft -f result_roary -p 40 -r -i 95 -cd 99 gff_file/*.gff

# Selezioniamo i geni che compongono il core-genome:
cp -r "result_roary/"* "gff_for_roary/"
cd gff_file
query_pan_genome -o query_pan_genome_result_core -a intersection *.gff
mv query_pan_genome_result_core  result_roary
find -type f ! -name "*.gff" -exec rm {} +
cd ..

panaroo -i gff_file/*.gff -o results_panaroo --clean-mode strict --remove-invalid-genes -a core --aligner mafft --core_threshold 0.99 \
        --merge_paralogs --refind_prop_match 0.5

mkdir -p ppanggolin
cd ppanggolin 
ppanggolin workflow --anno genomes.gbff.list --identity 0.95 --mode 2 


##############################################
# STEP 7: Phylogenetic analysis - Phylophlan #
##############################################

mkdir -p filogenesi 
cd filogenesi 
mkdir -p logs 
cd ..

# Generiamo un database di markers per il batterio in esame

phylophlan_setup_database -g s__Collinsella_aerofaciens -o . --verbose 2>&1 | tee filogenesi/logs/phylophlan_setup_database.log

# Impostiamo il file di configurazione

phylophlan_write_config_file -o filogenesi/file_config.cfg -d a --force_nucleotides --db_aa diamond --map_aa diamond --map_dna diamond --msa mafft --trim trimal --tree1 fasttree

# Procediamo con la costruzione dell'albero

phylophlan -i genomi_checkm -o filogenesi/output_phylogeny -d s__Collinsella_aerofaciens -t a -f filogenesi/file_config.cfg --nproc 20 --diversity low --force_nucleotides --accurate --verbose 2>&1 |tee filogenesi/logs/phylophlan__s__Collinsella_aerofaciens.log

# Costruzione dell'albero con bootstrap
iqtree -s genomi_checkm_concatenated.aln -m MFP -bb 1000 -T AUTO


#######################
# STEP 8: 16S Analysis #
######################à

input_dir="genomi_checkm"
mkdir -p 16s
mkdir -p tmp_16s
tmp_dir="tmp_16s"
out_dir="16s"

for file in "$input_dir"/*.fasta;do
    base_name=$(basename "$file" .fasta)
    barrnap --kingdom bac --outseq "$tmp_dir"/${base_name}.fasta $file
done

for file_16 in "$tmp_dir"/*.fasta;do
    genome_name=$(basename "$file_16" .fasta)
    grep -A 1 "16S" "$file_16" | grep -v "^--" > "$out_dir"/${genome_name}_16S.fasta
done

# Selezioniamo la sequenza più lunga per ciascun file 
mkdir -p selected_16s
mkdir -p tmp
for file in "16s"/*.fasta; do
    base=$(basename  "$file" _16S.fasta)
    rg -B1 -n "^.{$(wc -L < $file)}$" $file > "tmp/tmp_16s.fasta"
    head -n 2 "tmp/tmp_16s.fasta" | sed 's/^[0-9]*-//' | sed 's/^[0-9]*://' > "selected_16s/${base}.fasta"
    rm tmp/tmp_16s.fasta
done

output_file="16s/16s_file_tot.fasta"
for file in "selected_16s"/*.fasta;do
    base=$(basename "$file" .fasta | cut -d'_' -f1-2)
    sed "s/^>.*$/>$base/" "$file" >> "$output_file"
done 

clustalo --infile /home/biouserc/16_TOTALE.fasta --guidetree-out clustalo-prova_auto.dnd --distmat-out prova_clust_auto.matrix --full --outfmt clustal --auto --outfile aln_clustal_auto

########################
# Step 9: ADH Analysis #
########################

datasets summary genome taxon "Collinsella aerofaciens" --assembly-source refseq --exclude-atypical --mag exclude --as-json-lines > collinsella_summary.jsonl
cat collinsella_summary.jsonl | jq '.accession' > accession_numbers.txt
sed 's/"//g' accession_numbers.txt > GCF.txt 
rm accession_numbers.txt
mkdir -p PROTEINE
ncbi-genome-download --assembly-accessions GCF.txt bacteria --formats protein-fasta --flat-output --output-folder PROTEINE -v

cat PROTEINE/*.faa > sequenze_proteine_totali.fasta
makeblastdb -in sequenze_proteine_totali.fasta -dbtype prot -out db_blast_protein

blastp -query ADH_protein.faa -db db_blast_protein -out risultati_blast.txt -outfmt 6 

python filter_results_blast.py

cut -d $'\t' -f2,3 "filtered_blast_results_ref_C_a_4.tsv" > tmp_2.txt

python ADH_to_csv.python

####################
# STEP 10: CAZYmes #
####################


input_dir="genomi_checkm"

mkdir -p cazyme
output_dir="cazyme"


for file in "$input_dir"/*.fasta; do
    # Estrai il nome del file senza estensione
    filename=$(basename "$file" .fasta)
    run_dbcan "$file" prok --out_dir "$output_dir/${filename}" --verbose --dbcan_thread 4
done

mkdir -p overview
mkdir -p dbcan

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



