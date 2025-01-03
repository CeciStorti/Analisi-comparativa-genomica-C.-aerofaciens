#!/bin/bash
#CALCOLO DEL POCP

directory_genome=$1
num_threads=$2



# Creiamo un array dove salvare i file contenuti nella cartella genome
genomes=()
# Riempiamo l'array
for file in "$directory_genome"/*.faa; do
    if [ -f "$file" ]; then
        genomes+=("$file")
    fi
done

# Matrice PCOP
matrix_file="pcop_matrix.tsv"
echo -n -e "\t" > "$matrix_file"
for ((j=0; j<${#genomes[@]}; j++)); do
    genome=${genomes[j]}
    column_name="${genome##*/}" 
    column_name="${column_name%.faa}" 
    if [[ $j -lt $((${#genomes[@]} -1)) ]]; then
        echo -n -e "$column_name\t" >> "$matrix_file"
    else 
        echo -n -e "$column_name" >> "$matrix_file" 
    fi 
done
echo "" >> "$matrix_file"


for ((i=0; i<${#genomes[@]}; i++)); do
    genome1=${genomes[i]}
    row_name="${genome1##*/}"
    row_name="${row_name%.faa}"
    echo -n -e "$row_name\t" >> "$matrix_file"

    for ((j=i; j<${#genomes[@]}; j++)); do  # Inizia da "i" per evitare calcoli duplicati
        genome2=${genomes[j]}

        if [ $i -eq $j ]; then
            echo -n -e "100\t" >> "$matrix_file"  # 100 per la diagonale (genome vs genome)
            continue
        fi

        cp "$genome1" genes_A.faa_HH
        cp "$genome2" genes_B.faa_HH

        # Make BLAST index files.
        makeblastdb -in genes_A.faa_HH -dbtype prot -parse_seqids 1>/dev/null
        makeblastdb -in genes_B.faa_HH -dbtype prot -parse_seqids 1>/dev/null

        # Get total gene count.
        grep -c '>' genes_A.faa_HH > count_T1_HH
        grep -c '>' genes_B.faa_HH > count_T2_HH

        # BLAST genes.
        blastp -query genes_A.faa_HH -db genes_B.faa_HH -outfmt "6 std qlen" -out BLAST_A_against_B_HH -num_threads "$num_threads" -max_target_seqs 1 
        blastp -query genes_B.faa_HH -db genes_A.faa_HH -outfmt "6 std qlen" -out BLAST_B_against_A_HH -num_threads "$num_threads" -max_target_seqs 1 

        # Make sure each gene has only one hit.
        awk '{i=$1;if(i!=j){print};j=$1}' BLAST_A_against_B_HH > BLAST_A_against_B_top_hit_HH
        awk '{i=$1;if(i!=j){print};j=$1}' BLAST_B_against_A_HH > BLAST_B_against_A_top_hit_HH

        # Count genes above the cut-offs: identity = 40%; e-value = 1e-05; length of query = 50%
        awk '$3>=40 && $4/$13>=.5 && $11<1e-05' BLAST_A_against_B_top_hit_HH | wc -l > count_C1_HH
        awk '$3>=40 && $4/$13>=.5 && $11<1e-05' BLAST_B_against_A_top_hit_HH | wc -l > count_C2_HH

        # Calculate POCP value.
        pcop=$(awk 'BEGIN{print ('`cat count_C1_HH`'+'`cat count_C2_HH`')/('`cat count_T1_HH`'+'`cat count_T2_HH`')*100}')
        echo -n -e "$pcop\t" >> "$matrix_file"

        # Rimuoviamo i file temporanei
        rm count_T1_HH count_T2_HH count_C1_HH count_C2_HH 
        rm BLAST_A_against_B_HH BLAST_B_against_A_HH
        rm BLAST_A_against_B_top_hit_HH BLAST_B_against_A_top_hit_HH
        rm genes_A.faa_HH genes_B.faa_HH genes_A.faa_HH.p* genes_B.faa_HH.p*
    done

    echo "" >> "$matrix_file"
done


# Generazione del file triangolare_sup.csv
# (Il tuo codice che genera il file...)

# Creazione della matrice simmetrica usando Python
python3 << EOF
import pandas as pd
import numpy as np


mat = pd.read_csv("pcop_matrix.tsv", sep='\t')

mat.columns.values[0] = "Unnamed"
mat = mat.replace(r'^\s*$', np.nan, regex=True)
df = mat.rename(columns={"Unnamed": "index"}).set_index("index") 
print("Dataframe iniziale:")
print(df)

#  Matrice simmetrica
matrix = df.to_numpy()
matrix = np.where(np.isnan(matrix), matrix.T, matrix)  


symmetric_df = pd.DataFrame(matrix, index=df.index, columns=df.columns)


symmetric_df.to_csv("matrice_simmetrica.csv")
EOF


