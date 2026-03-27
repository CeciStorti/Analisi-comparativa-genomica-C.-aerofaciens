import pandas as pd

def create_presence_absence_matrix(input_file, output_file):

    df = pd.read_csv(input_file, sep='\t')

    df["presence"] = df["pident"].apply(lambda x: 1 if x > 70 else 0)

    result = df[["genome_id", "query_id", "presence"]]

    # Wide format
    matrix = result.pivot(index="genome_id", columns="query_id", values="presence")

    matrix = matrix.fillna(0)

    # 👉 salva in formato Excel
    matrix.to_excel(output_file)

    print(f"Matrix saved in {output_file}")


input_file = "filtered_blast_results.tsv"
output_file = "presence_ADH.xlsx"

create_presence_absence_matrix(input_file, output_file)
