import pandas as pd
import os

# -----------------------------
# PARAMETERS
# -----------------------------
SORTED_CSV = "sorted_gene_table.csv"
OUTPUT_CSV = "identified_clusters.csv"

def run(output, distance, length, numberOfGenes):
    distance = int(distance)
    length = int(length)
    numberOfGenes = int(numberOfGenes)
    # -----------------------------
    # LOAD SORTED CSV
    # -----------------------------
    df = pd.read_csv(SORTED_CSV)
    print(f"Loaded sorted CSV with {len(df)} rows.")


    # -----------------------------
    # CLUSTER IDENTIFICATION
    # -----------------------------
    clusters = []   # list of list-of-rows (cluster)
    prev_row = None

    current_cluster = []
    cluster_start = None
    cluster_end = None


    for idx, row in df.iterrows():

        gene_start = row["Start"]
        gene_end = row["End"]

        if prev_row is None:
            # Start first cluster
            current_cluster = [row]
            cluster_start = gene_start
            cluster_end = gene_end

        else:
            same_scaffold = (row["Scaffold ID"] == prev_row["Scaffold ID"])
            same_assembly = (row["Assembly ID"] == prev_row["Assembly ID"])

            if same_scaffold and same_assembly:

                # distance between this gene & previous
                adj_dist = gene_start - prev_row["End"]

                # cluster span if we add this gene
                new_span = gene_end - cluster_start

                if adj_dist <= distance and new_span <= length:
                    # Add gene
                    current_cluster.append(row)
                    cluster_end = gene_end
                else:
                    # Save old cluster
                    clusters.append(current_cluster)
                    # start a new one
                    current_cluster = [row]
                    cluster_start = gene_start
                    cluster_end = gene_end

            else:
                # scaffold/assembly changed
                clusters.append(current_cluster)
                current_cluster = [row]
                cluster_start = gene_start
                cluster_end = gene_end

        prev_row = row

    # append last cluster
    clusters.append(current_cluster)


    # -----------------------------
    # BUILD OUTPUT WITH DISTANCES
    # -----------------------------
    cluster_rows = []
    cluster_id = 1

    for cluster in clusters:

        # Filter out tiny clusters
        if len(cluster) < numberOfGenes:
            continue

        # Compute total span for the cluster
        cluster_first = cluster[0]["Start"]
        cluster_last = cluster[-1]["End"]
        cluster_span = cluster_last - cluster_first

        prev_gene_end = None

        for r in cluster:
            # compute distance to previous gene in the same cluster
            if prev_gene_end is None:
                adj_dist = None   # first gene → no distance
            else:
                adj_dist = r["Start"] - prev_gene_end

            prev_gene_end = r["End"]

            cluster_rows.append({
                "Cluster ID": cluster_id,
                "Gene ID": r["Gene ID"],
                "Assembly ID": r["Assembly ID"],
                "Scaffold ID": r["Scaffold ID"],
                "Start": r["Start"],
                "End": r["End"],
                "Adjacent Distance": adj_dist,
                "Cluster Span": cluster_span
            })

        cluster_id += 1


    # Save output
    output_dir = output
    output_file = "results.csv"
    path = os.path.join(output_dir, output_file)
    df_out = pd.DataFrame(cluster_rows)
    df_out.to_csv(path, index=False)

    print("✔ Cluster identification complete!")
    print("✔ Removed clusters with < " + str(numberOfGenes) + "genes.")
    print("✔ Output saved to:", path)
    print(df_out.head(15))

if __name__ == "__main__":
    run()