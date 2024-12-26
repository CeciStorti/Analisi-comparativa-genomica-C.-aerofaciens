import pandas as pd


import pandas as pd


blast_results = pd.read_csv("risultati_blast_ref_C_a_4.txt", sep="\t", header=None)

# Aggiungi nomi alle colonne per semplicità (puoi cambiarli in base al tuo file)
blast_results.columns = ["query_id", "subject_id", "pident", "length", "mismatch", 
                         "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


blast_results["genome_id"] = blast_results["subject_id"].str.extract(r"(GCF_[^_]+)")


filtered_results = (
    blast_results.loc[blast_results.groupby("genome_id")["pident"].idxmax()]
)


filtered_results.to_csv("filtered_blast_results_ref_C_a_4.tsv", sep="\t", index=False)


