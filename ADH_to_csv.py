
import pandas as pd

def create_presence_absence_matrix(input_file, output_file, gene_name):
    # Leggi il file di input
    df = pd.read_csv(input_file,sep='\t')

    # Crea la colonna del gene con valori 1 o 0 in base alla quantità
    df[gene_name] = df["pident"].apply(lambda x: 1 if x > 70 else 0)

    # Seleziona solo le colonne necessarie (ID e gene_name)
    result = df[["subject_id", gene_name]]

    # Salva la matrice in un nuovo file CSV
    result.to_csv(output_file, index=False)

    print(f"Matrice di presenza/assenza salvata in {output_file}")

# Esempio di utilizzo
input_file = "/home/biouserc/ADH_blast/ref_C_a_4/tmp_2.txt"  
output_file = "/home/biouserc/ADH_blast/ref_C_a_4/ADH_4.csv" 
gene_name = "ADH_4"  

create_presence_absence_matrix(input_file, output_file, gene_name)