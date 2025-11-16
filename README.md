## ConCompRA Post-Processing Script (`concompra.py`)

This script processes **ConCompRA OTU results** to produce clean, normalized species abundance tables, including unmapped reads if available. 

> ⚠️ **Important:** This script is designed to run **after you have run ConCompRA**. It expects each dataset folder to contain `otu_table.csv` and `clustered_consensus.fasta`.

---

### Requirements

1. **Input files per dataset folder**
   - `otu_table.csv` → OTU table produced by ConCompRA.
   - `clustered_consensus.fasta` → Consensus sequences produced by ConCompRA.
   - `taxonomy.tsv` → Reference taxonomy table (place in the root directory).
   - `species_taxid.fasta` → Reference species database (place in the root directory).
   - Optionally, `unmapped/*.fastq` → Unmapped reads for each sample, if available.

   Both `taxonomy.tsv` and `species_taxid.fasta` can be downloaded from EMU’s GitLab repository:  
   [EMU Database on GitLab](https://gitlab.com/treangenlab/emu/-/tree/master/emu_database?ref_type=heads)

2. **Software dependencies**
   - Python 3 with `pandas` installed
   - `vsearch` → Performs sequence matching
   - Standard Unix tools like `awk` (usually available on Linux/macOS)

---

### How the Script Works

1. Searches recursively for folders containing `otu_table.csv` and `clustered_consensus.fasta`.
2. Copies required files into a **temporary working directory** (original files remain unchanged).
3. Uses `vsearch` to map OTUs to `species_taxid.fasta`.
4. Combines vsearch results with `taxonomy.tsv` to classify sequences.
5. For each sample:
   - Normalizes OTU abundances.
   - Removes duplicate species.
   - Produces a **mapped-only TSV** in `output/`.
6. Optionally adds unmapped reads:
   - Counts unmapped reads from `unmapped/*.fastq`.
   - Recalculates relative abundances including unmapped reads.
   - Produces a **TSV with unmapped reads** in `output_with_unmapped/`.

---

### Output

- **Mapped-only TSVs:** Normalized abundances per species for each sample.
- **Mapped + Unmapped TSVs:** Includes a row for unmapped reads, showing their relative abundance.
- Output files are organized in `output/` and `output_with_unmapped/` subfolders per dataset folder.

---

### Summary

- **Input:** ConCompRA outputs (`otu_table.csv`, `clustered_consensus.fasta`), `taxonomy.tsv`, `species_taxid.fasta`
- **Processing:** vsearch classification, merge taxonomy, normalize abundances, handle unmapped reads
- **Output:** Clean TSV tables ready for downstream analysis

✅ **Tip:** Do not modify original ConCompRA outputs. The script safely works in a temporary directory and preserves your original data.
