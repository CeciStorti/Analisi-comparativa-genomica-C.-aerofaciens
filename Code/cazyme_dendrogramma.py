import pandas as pd
import os
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.cluster import hierarchy


folder_path = "overview"  


file_list = [f for f in os.listdir(folder_path) if f.endswith('.txt')]


df_list = []


for file in file_list:
    
    df = pd.read_csv(os.path.join(folder_path, file), sep="\t")
    genome_name = file.rsplit('.',1)[0].replace('annotations_', '')  
    df['Genome'] = genome_name
    df_list.append(df)

# Combine all DataFrames in one
df_combined = pd.concat(df_list, ignore_index=True)

# Generate the presence/absence matrix 
df_relevant = df_combined[['Genome', 'DIAMOND']]
df_relevant = df_relevant.drop_duplicates()
df_presence_absence = df_relevant.groupby(['Genome', 'DIAMOND']).size().unstack(fill_value=0)


df_final=df_presence_absence.iloc[1:]

df_final

df_final.to_csv('dbcan/presence_absence_matrix.csv')

# DENDROGRAM:

data=pd.read_csv("dbcan/presence_absence_matrix.csv")


distance_matrix = pdist(data, metric='jaccard')

Z = linkage(distance_matrix, method='complete')  

# Cluster calculation
threshold = 0.6 * max(Z[:, 2])  
clusters = fcluster(Z, t=threshold, criterion='distance')


unique_clusters = np.unique(clusters)
cluster_colors = ["steelblue", "orange", "orangered", "violet","lightblue","lightgreen","cadetblue"]  
color_map = {cluster: cluster_colors[i % len(cluster_colors)] for i, cluster in enumerate(unique_clusters)}
hierarchy.set_link_color_palette(["steelblue", "orange", "orangered", "violet","lightblue","lightgreen","cadetblue"])
# Plot 
plt.figure(figsize=(9, 6))
dendro = hierarchy.dendrogram(Z, color_threshold=threshold, above_threshold_color='gray',no_labels=True)

# Legend
cluster_sizes = {cluster: list(clusters).count(cluster) for cluster in unique_clusters}
legend_elements = [Patch(facecolor=color_map[cluster], label=f'Cluster {cluster}: {size} elements')
                   for cluster, size in cluster_sizes.items()]
plt.legend(handles=legend_elements, loc='upper right', title="Clusters and dimensions")

plt.title('Dendrogram based on the Cazymes presence/absence matrix')
plt.xlabel('Genomes')
plt.ylabel('Distance (Jaccard)')
plt.show()

