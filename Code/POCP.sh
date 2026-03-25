#!/bin/bash
# CALCULATION OF POCP

directory_genome=$1
num_threads=$2

# Create an array to store the files contained in the genome folder
genomes=()
# Fill the array
for file in "$directory_genome"/*.faa; do
    if [ -f "$file" ]; then
        genomes+=("$file")
    fi
done

# PCOP matrix
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

    for ((j=i; j<${#genomes[@]}; j++)); do  # Start from "i" to avoid duplicate calculations
        genome2=${genomes[j]}

        if [ $i -eq $j ]; then
            echo -n -e "100\t" >> "$matrix_file"  # 100 for the diagonal (genome vs genome)
            continue
        fi

        cp "$genome1" genes_A.faa_HH
        cp "$genome2" genes_B.faa_HH

        # Make BLAST index files
        makeblastdb -in genes_A.faa_HH -dbtype prot -parse_seqids 1>/dev/null
        makeblastdb -in genes_B.faa_HH -dbtype prot -parse_seqids 1>/dev/null

        # Count total genes
        grep -c '>' genes_A.faa_HH > count_T1_HH
        grep -c '>' genes_B.faa_HH > count_T2_HH

        # Run BLASTP: compare genome A proteins with genome B and vice versa
        blastp -query genes_A.faa_HH -db genes_B.faa_HH -outfmt "6 std qlen" -out BLAST_A_against_B_HH -num_threads "$num_threads" -max_target_seqs 1 
        blastp -query genes_B.faa_HH -db genes_A.faa_HH -outfmt "6 std qlen" -out BLAST_B_against_A_HH -num_threads "$num_threads" -max_target_seqs 1 

        # Select the highest matches for each gene from BLAST results
        awk '{i=$1;if(i!=j){print};j=$1}' BLAST_A_against_B_HH > BLAST_A_against_B_top_hit_HH
        awk '{i=$1;if(i!=j){print};j=$1}' BLAST_B_against_A_HH > BLAST_B_against_A_top_hit_HH

        # Filter matches based on percent identity and alignment coverage
        awk '$3>=40 && $4/$13>=.5 && $11<1e-05' BLAST_A_against_B_top_hit_HH | wc -l > count_C1_HH
        awk '$3>=40 && $4/$13>=.5 && $11<1e-05' BLAST_B_against_A_top_hit_HH | wc -l > count_C2_HH

        # POCP calculation:
        # Defined as: POCP = (C1 + C2) / (T1 + T2) * 100, where:
        # C1 = genes of A found in B
        # C2 = genes of B found in A
        # T1 = total number of genes in A
        # T2 = total number of genes in B
        pcop=$(awk 'BEGIN{print ('`cat count_C1_HH`'+'`cat count_C2_HH`')/('`cat count_T1_HH`'+'`cat count_T2_HH`')*100}')
        echo -n -e "$pcop\t" >> "$matrix_file"

        # Remove temporary files
        rm count_T1_HH count_T2_HH count_C1_HH count_C2_HH 
        rm BLAST_A_against_B_HH BLAST_B_against_A_HH
        rm BLAST_A_against_B_top_hit_HH BLAST_B_against_A_top_hit_HH
        rm genes_A.faa_HH genes_B.faa_HH genes_A.faa_HH.p* genes_B.faa_HH.p*
    done

    echo "" >> "$matrix_file"
done

# Create symmetric matrix using Python
python3 << EOF
import pandas as pd
import numpy as np

mat = pd.read_csv("pcop_matrix.tsv", sep='\t')

mat.columns.values[0] = "Unnamed"
mat = mat.replace(r'^\s*$', np.nan, regex=True)
df = mat.rename(columns={"Unnamed": "index"}).set_index("index") 
print("Initial dataframe:")
print(df)

# Symmetric matrix
matrix = df.to_numpy()
matrix = np.where(np.isnan(matrix), matrix.T, matrix)  

symmetric_df = pd.DataFrame(matrix, index=df.index, columns=df.columns)

symmetric_df.to_csv("POCP_matrix.csv")
EOF
EOF


