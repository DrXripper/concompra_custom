import os
import subprocess
import tempfile
import shutil
import pandas as pd

# Root directory
root_dir = os.getcwd()

taxonomy_file = os.path.join(root_dir, "taxonomy.tsv")
species_db = os.path.join(root_dir, "species_taxid.fasta")

# Traverse directories
for dirpath, dirnames, filenames in os.walk(root_dir):
    if "otu_table.csv" in filenames and "clustered_consensus.fasta" in filenames:
        results_dir = dirpath
        print(f"\n🔍 Found results folder: {results_dir}")

        otu_path = os.path.join(results_dir, "otu_table.csv")
        fasta_path = os.path.join(results_dir, "clustered_consensus.fasta")

        output_dir = os.path.join(results_dir, "output")
        output_unmapped_dir = os.path.join(results_dir, "output_with_unmapped")
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(output_unmapped_dir, exist_ok=True)

        parent_of_results = os.path.dirname(results_dir)
        unmapped_dir = os.path.join(parent_of_results, "unmapped")

        with tempfile.TemporaryDirectory() as tmpdir:
            shutil.copy(otu_path, os.path.join(tmpdir, "otu_table.csv"))
            shutil.copy(fasta_path, os.path.join(tmpdir, "clustered_consensus.fasta"))
            shutil.copy(taxonomy_file, os.path.join(tmpdir, "taxonomy.tsv"))
            shutil.copy(species_db, os.path.join(tmpdir, "species_taxid.fasta"))

            bash_script = """
            vsearch --usearch_global clustered_consensus.fasta \\
                --db species_taxid.fasta \\
                --id 0.9 \\
                --userout hits.txt \\
                --userfields query+target+id

            awk '
                BEGIN {
                    FS=OFS="\\t"
                    while ((getline < "taxonomy.tsv") > 0) {
                        if (NR==1) {
                            header = $0
                            continue
                        }
                        tax[$1] = $0
                    }
                }
                {
                    split($2, parts, ":")
                    tid = parts[1]
                    identity = $3
                    if (tid in tax) {
                        print $1, identity, tax[tid]
                    } else {
                        print $1, identity, tid, "UNKNOWN"
                    }
                }
            ' hits.txt > classified_full_tmp.tsv

            {
                printf "OTU\\tidentity\\t%s\\n", "$(head -n1 taxonomy.tsv)"
                cat classified_full_tmp.tsv
            } > classified_taxonomy.tsv

            rm hits.txt classified_full_tmp.tsv
            """
            subprocess.run(bash_script, shell=True, check=True, cwd=tmpdir)

            classified_df = pd.read_csv(os.path.join(tmpdir, "classified_taxonomy.tsv"), sep="\t")
            otu_df = pd.read_csv(os.path.join(tmpdir, "otu_table.csv"))

            for sample in otu_df.columns[1:]:
                print(f"📁 Processing sample: {sample}")
                sample_df = otu_df[["#OTU ID", sample]].copy()
                sample_df.columns = ["OTU", "abundance"]

                mapped_total = sample_df["abundance"].sum()
                sample_df["abundance"] = sample_df["abundance"] / mapped_total

                merged = pd.merge(classified_df, sample_df, on="OTU", how="inner")
                merged_clean = (
                    merged.drop_duplicates(subset=["species"])
                    if "species" in merged.columns else merged
                )

                sample_name = sample.replace(".CONCOMPRA", "")
                output_path = os.path.join(output_dir, f"{sample_name}.tsv")
                merged_clean.to_csv(output_path, sep="\t", index=False)
                print(f"    ✅ Saved mapped-only: {output_path}")

                # =========================
                # ADD UNMAPPED ROW SECTION
                # =========================

                df = pd.read_csv(output_path, sep="\t")

                # ✅ Get actual raw counts per OTU
                raw_abundance = otu_df.set_index("#OTU ID")[sample]
                df["raw_count"] = df["OTU"].map(raw_abundance.to_dict())
                df = df.dropna(subset=["raw_count"])

                # ✅ Count unmapped reads
                unmapped_fastq = os.path.join(unmapped_dir, f"{sample_name}.fastq")
                if os.path.exists(unmapped_fastq):
                    with open(unmapped_fastq, "r") as f:
                        unmapped_lines = sum(1 for _ in f)
                    unmapped_read_count = unmapped_lines // 4
                else:
                    unmapped_read_count = 0
                    print(f"⚠️  No unmapped fastq for {sample_name}, using 0")

                mapped_read_count = df["raw_count"].sum()
                total_reads = mapped_read_count + unmapped_read_count

                df["abundance"] = df["raw_count"] / total_reads
                df.drop(columns=["raw_count"], inplace=True)

                # ✅ Add unmapped row
                unmapped_row = pd.Series({col: "" for col in df.columns})
                unmapped_row["species"] = "unmapped"
                unmapped_row["abundance"] = unmapped_read_count / total_reads
                df = pd.concat([df, pd.DataFrame([unmapped_row])], ignore_index=True)

                # ✅ Final check and save
                final_sum = df["abundance"].astype(float).sum()
                print(f"    ✅ With unmapped: normalized sum = {final_sum:.10f}")

                output_unmapped_path = os.path.join(output_unmapped_dir, f"{sample_name}.tsv")
                df.to_csv(output_unmapped_path, sep="\t", index=False)
                print(f"    📂 Saved with unmapped: {output_unmapped_path}")

print("\n🎉 All ConCompRA results processed successfully with and without unmapped rows.")
