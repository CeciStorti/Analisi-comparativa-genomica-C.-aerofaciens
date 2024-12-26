import pandas as pd

# Leggi il file .tsv
df = pd.read_csv('/home/cstorti/qLS/Collinsella/collinsella_filtered/biosample_info.tsv', sep='\t')

# Trasforma il formato long in wide
df_wide = df.pivot(index=['Assembly Accession','Assembly BioSample Accession'], 
                   columns='Assembly BioSample Attribute Name', 
                   values='Assembly BioSample Attribute Value')

# Resetta l'indice se necessario
df_wide.reset_index(inplace=True)

# Scrivi il file .tsv di output
df_wide.to_csv('biosample_info_wide.tsv', sep='\t', index=False)

print("Trasformazione completata. File salvato come 'output_wide.tsv'.")
