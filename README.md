# Metadata for Viral RefSeq Genomes

This script creates a metadata file with taxonomy annotations for each genome in RefSeq.

## Dependencies

The script requires the following Python packages:

- Biopython
- pytaxonkit
- pandas
- tqdm

You can install the dependencies using `mamba`, `conda`, or `pip`.

### Mamba

```bash
mamba install -c conda-forge -c bioconda biopython pytaxonkit pandas tqdm
```

### Conda

```bash
conda install -c conda-forge -c bioconda biopython pytaxonkit pandas tqdm
```

### Pip

```bash
pip install biopython pytaxonkit pandas tqdm
```

## Usage

The script takes the following arguments:

- `--refseq_nucl_file`: The RefSeq nucleotide file (required).
- `--refseq_version`: The RefSeq version (e.g., 211, 216, etc., required).
- `--output_folder`: The name of the output folder. Defaults to the same folder as the input fasta file.
- `--taxonomy_database`: Path to the NCBI taxonomy database folder (required).