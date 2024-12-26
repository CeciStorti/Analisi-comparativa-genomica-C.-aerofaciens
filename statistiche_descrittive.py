import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv("statistiche_checkm.tsv",sep="\t")

df_selected=df[['Genome size (bp)','Completeness','Contamination','GC','# scaffolds','N50 (scaffolds)','# contigs',
               'N50 (contigs)','# predicted genes']].copy()

df_selected['Genome size (bp)']=df['Genome size (bp)']/1000000
df_selected['Genome size (bp)']=df_selected['Genome size (bp)'].astype(float)
df_selected.rename(columns={'Genome size (bp)':'Genome size (Mbp)'},inplace=True)
# Dimensione della figura
fig, axes = plt.subplots(3, 3, figsize=(20, 15))


columns = df_selected.columns
y_labels = ['Genome Size (Mbp)', 'Completeness (%)', 'Contamination (%)', 'G+C content (%)',
            'Number of Scaffold Count', 'N50 Scaffold Size (Mbp)', 'Number of contigs',
            'N50 Contig Size (Mbp)', 'Number of Predicted Genes']
x_labels=['Genomes']*len(columns)

# Crea un subplot per ogni colonna del DataFrame
for i, col in enumerate(columns):
    ax = axes[i // 3, i % 3]  # Definisce la posizione del subplot
    sns.boxplot(data=df_selected[col], ax=ax,color='lightgrey',medianprops={"color": "black", "linewidth": 1.5})  # Crea il boxplot per la colonna
    #ax.set_title(f'Boxplot of {titles[i]}') 
    ax.set_xlabel(x_labels[i],fontsize=15)
    ax.set_ylabel(y_labels[i],fontsize=15) 
    sns.stripplot(data=df_selected[col], ax=ax, color="red", size=5, jitter=True)
    for spine in ax.spines.values(): 
        spine.set_visible(True)  
        spine.set_linewidth(1) 
        spine.set_linestyle('--')
    
 


plt.tight_layout()


plt.show()
plt.savefig("./Boxplots caratteristiche dei genomi.png")


# Analizziamo i dati relativi all'area geografica e alla provenienza dei campioni.

file_path = "biosample_info_wide.tsv"
biosample_data = pd.read_csv(file_path, sep='\t')

biosample_data['geo_loc_name']=biosample_data['geo_loc_name'].replace({'not provided':'na','not determined':'na','unknown':'na','not applicable':'na','NaN':'na'})
biosample_data['geo_loc_name']=biosample_data['geo_loc_name'].replace({'USA:Boston':'USA','USA: New York':'USA','USA:Chicago':'USA','USA:New York':'USA','USA: New York City':'USA','China: Shenzhen':'China'})
df=biosample_data['geo_loc_name'].value_counts()

# Barplot Area geografica
df.plot(kind='bar', color='lightgreen')

# Titolo ed etichette
plt.title('Barplot Aree geografiche')
plt.xlabel('Area geografica')
plt.ylabel('Count')


plt.show()
plt.savefig("./Istogramma_geogrphic_area.png")

# Alcuni campi non presentano lo stesso nome, ma dei sinonimi. Uniformiamo i contenuti

biosample_data['isolation_source']=biosample_data['isolation_source'].replace({'stool':'feces','human feces': 'human fecal sample','missing':'na','Korean adult feces':'human fecal sample','fecal sample':'feces'})

# Sostituiamo i valori con NaN con 'na'
biosample_data['isolation_source'] = biosample_data['isolation_source'].fillna('na')
df_isolation=biosample_data['isolation_source'].value_counts()


# Barplot Isolation source
df_isolation.plot(kind='bar', color='lightblue')


plt.xlabel('Isolation source')
plt.ylabel('Count')

plt.show()
plt.savefig("./Istogramma_isolation_source.png")