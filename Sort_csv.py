import pandas as pd
import os
import numpy as np
import ast


def run():
    with open("Scaffold_dictionary.txt", "r", encoding="utf-8") as f:
        text = f.read()


    Scaffold_dictionary = ast.literal_eval(text)

    print("Loaded dictionary with", len(Scaffold_dictionary), "PDB entries")


    with open("Assembly_dictionary.txt", "r", encoding="utf8") as f:
        text = f.read()

    assembly_dictionary = ast.literal_eval(text)

    print("Loaded", len(assembly_dictionary), "entries.")

    with open("coordinates.txt", "r", encoding="utf-8") as f:
        text = f.read()
    coordinates = ast.literal_eval(text)
    print("Loaded", len(coordinates), "coordinate entries.")



    #### File Names ###
    raw_csv = "final_gene_scaffold_assembly_coordinates.csv"
    sorted_csv = "sorted_gene_table.csv"



    gene_to_coords = {}

    for pdb_file, genes in coordinates.items():
        for gene_id, coord_tuple in genes.items():
            # convert gene_id to int for consistency
            gene_to_coords[int(gene_id)] = coord_tuple


    if os.path.exists(raw_csv):
        print(f"Found existing file: {raw_csv}")
        print("Skipping DataFrame creation.")
        df = pd.read_csv(raw_csv)

    else:
        print("CSV not found — creating new full table...")

        rows = []

        # Loop through each PDB file
        for pdb_file in Scaffold_dictionary:

            scaffold_map = Scaffold_dictionary[pdb_file]        # gene → scaffold
            assembly_map = assembly_dictionary.get(pdb_file, {})  # gene → assembly
            coord_map = coordinates.get(pdb_file, {})           # gene → (start, end)

            for gene_id in scaffold_map:

                scaffold_id = scaffold_map.get(gene_id, "")
                assembly_id = assembly_map.get(gene_id, "")
                start, end = coord_map.get(gene_id, ("", ""))

                rows.append({
                    "Gene ID": gene_id,
                    "PDB File": pdb_file,
                    "Scaffold ID": scaffold_id,
                    "Assembly ID": assembly_id,
                    "Start": start,
                    "End": end
                })

        # Build DataFrame
        df = pd.DataFrame(rows)

        # Save
        output_file = "final_gene_scaffold_assembly_coordinates.csv"
        df.to_csv(output_file, index=False)

        print(f"\nSaved CSV file: {output_file}")
        print(df.head())
        print("\nTotal rows:", len(df))

    ######## Sorting Dataframe ############

    df = pd.read_csv("final_gene_scaffold_assembly_coordinates.csv")

    df_sorted = df.sort_values(
        by=["Assembly ID", "Scaffold ID", "Start"],
        ascending=[True, True, True]
    ).reset_index(drop=True)

    ############## Genomic Distance ################
    df_sorted["next_start"] = df_sorted["Start"].shift(-1)
    df_sorted["next_scaffold"] = df_sorted["Scaffold ID"].shift(-1)
    df_sorted["next_assembly"] = df_sorted["Assembly ID"].shift(-1)

    def compute_distance(row):
        # Only compute if next gene is on same scaffold & same assembly
        if (
            row["Scaffold ID"] == row["next_scaffold"] and
            row["Assembly ID"] == row["next_assembly"]
        ):
            return row["next_start"] - row["End"]
        else:
            return np.nan

    df_sorted["genomic_distance"] = df_sorted.apply(compute_distance, axis=1)

    # remove helper columns
    df_sorted = df_sorted.drop(columns=["next_start", "next_scaffold", "next_assembly"])



    ################ Saving###############

    df_sorted.to_csv(sorted_csv, index=False)

    print(f"✔ Sorted CSV saved as {sorted_csv}")
    print(df_sorted.head(15))

if __name__ == "__main__":
    run()