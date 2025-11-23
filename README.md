
-----

# Genomic Language Model Pre-training Pipeline with Ladderpath Tokenization

## Project Overview

This repository provides a comprehensive end-to-end workflow for pre-training a genomic language model utilizing the **Ladderpath** tokenization algorithm. The pipeline begins with raw genomic sequences (FASTA format) and gene annotation files, proceeds through data preprocessing, parallelized Ladderpath computation, result merging, and vocabulary construction, and concludes with training a BERT-based Masked Language Model (MLM) using the Hugging Face Transformers library.

## Workflow

The pipeline consists of five primary steps, each handled by a dedicated Python script. Please execute them in the following order:

1.  **`split_chrom.py`**: **Data Preprocessing**

      * Reads a FASTA reference genome file and a CSV gene annotation file.
      * Extracts all gene sequences and intergenic non-coding sequences.
      * Handles overlapping genes and extracts sequences from chromosome ends.
      * Consolidates all extracted DNA sequences into a unified CSV file for the next step.

2.  **`Ladderpath_multiprocess_DNASequence.py`**: **Parallel Ladderpath Computation**

      * Reads the CSV file containing DNA sequences generated in the previous step.
      * Splits long sequences into manageable batches (default: max 1 million characters per batch).
      * Runs the Ladderpath algorithm in parallel using multiprocessing to analyze the hierarchical structure of each batch.
      * Saves the resulting Poset-Multiset (POM) data for each batch as individual JSON files.

3.  **`Merge_the_multiplicities.py`**: **Merging Parallel Results**

      * Reads all JSON files generated in Step 2.
      * Accumulates the occurrence counts (multiplicities) for each "Ladder unit" (Token).
      * Saves the final aggregated results into a single JSON file containing global statistics (`merge_multiplicities.json`).

4.  **`build_vocab.py`**: **Statistics & Vocabulary Construction**

      * Reads the aggregated JSON file (`merge_multiplicities.json`) from the previous step.
      * Sorts all tokens by frequency in descending order and selects the top $N$ tokens (e.g., 1,000) to construct the final vocabulary file `vocab.txt`.

5.  **`LPT_pretrain.py`**: **Model Pre-training**

      * A standard training script based on Hugging Face Transformers.
      * Uses the `vocab.txt` generated in Step 4 and pre-tokenized training data to pre-train a BERT model.
      * **Note**: This script requires the training data (`train_data_file`) to be a text file where each line consists of Token IDs (converted using `vocab.txt`) separated by spaces.

## Environment Setup

### 1\. Clone the Repository

```bash
git clone <your-repo-url>
cd <your-repo-name>
```

### 2\. Install Dependencies

It is recommended to use `conda` to create an isolated virtual environment.

```bash
conda create -n genomics_lm python=3.9
conda activate genomics_lm
```

Then, install the necessary libraries via `pip`.

```bash
pip install pandas==2.2.3 \
biopython==1.84 \
psutil==6.1.0 \
torch==2.6.0 \
transformers==4.47.1 \
tqdm==4.67.1
```

## Usage Instructions

### Step 1: Data Preprocessing

Ensure your reference genome file (e.g., `hg19.fa`) and gene annotation file (e.g., `hg19_zhengfu_sorted_output_file_end.csv`) are ready. Modify the file paths in `split_chrom.py`:

```python
# split_chrom.py
...
def main():
    fasta_file = r"./hg19.fa"  # <-- Update with your FASTA path
    csv_file = r"./hg19_zhengfu_sorted_output_file_end.csv"  # <-- Update with your annotation path
    output_file = r"./hg19_gene_annotation.csv" # <-- Output path
...
```

Run the script:

```bash
python split_chrom.py
```

Upon success, `hg19_gene_annotation.csv` will be generated.

### Step 2: Ladderpath Analysis

Modify `Ladderpath_multiprocess_DNASequence.py` to set the input file path, output directory, and the number of parallel workers.

```python
# Ladderpath_multiprocess_DNASequence.py
...
if __name__ == "__main__":
    max_length = 1000000  # Max total characters per batch
    output_folder = r'./hg19_json_file'  # <-- Output directory for POM JSON files
    file_path = r"./hg19_gene_annotation.csv"  # <-- Use the CSV generated in Step 1
    
    num_workers = 50 # <-- Adjust based on your CPU cores and memory
...
```

Run the script:

```bash
python Ladderpath_multiprocess_DNASequence.py
```

*Note: This process may take a significant amount of time. Upon completion, the `hg19_json_file` directory will contain numerous JSON files.*

### Step 3: Merge Parallel Results

Modify `Merge_the_multiplicities.py` to set the input directory and output file path.

```python
# Merge_the_multiplicities.py
...
# Define file paths
file_paths = [
    './hg19_json_file', # <-- Ensure this matches the output dir from Step 2
]
...
# Write results to JSON
output_file = './merge_multiplicities.json'  # <-- Output filename
...
```

Run the script:

```bash
python Merge_the_multiplicities.py
```

Upon success, `merge_multiplicities.json` will be generated containing global token frequencies.

### Step 4: Build Vocabulary

Modify `build_vocab.py` to **directly read the `merge_multiplicities.json` file** generated in the previous step (instead of re-scanning directories).

```python
# build_vocab.py (Recommended Modification)
import json

# Read the merged file directly
with open('./merge_multiplicities.json', 'r') as f:
    data_dict = json.load(f)

# Sort by frequency in descending order
sorted_items = sorted(data_dict.items(), key=lambda item: item[1], reverse=True)

# Extract Top K tokens (Modify K as needed)
top_k_keys = [key for key, value in sorted_items[:1000]]

# Write vocabulary
with open('./vocab.txt', 'w') as f:
    for key in top_k_keys:
        f.write(key + '\n')

print("Finished writing top keys to vocab.txt")
```

Run the modified script:

```bash
python build_vocab.py
```

Upon success, the final `vocab.txt` file will be generated.

### Step 5: Model Training

This is the final step. You must have the following files ready:

1.  **`vocab.txt`**: The vocabulary generated in Step 4.
2.  **`bert_config.json`**: Model configuration file (e.g., defining layers, hidden size, etc.).
3.  **Tokenized Training Data**: A `.txt` file where each line consists of space-separated Token IDs. *You will need to write a script to convert the sequences in `hg19_gene_annotation.csv` into this format using `vocab.txt`.*

Run the training script via command line. Example:

```bash
python LPT_pretrain.py \
    --output_dir ./lpt_model_output \
    --model_type bert \
    --do_train \
    --train_data_file ./path/to/your/tokenized_data.txt \
    --tokenizer_name ./ \
    --config_name ./bert_config.json \
    --mlm \
    --mlm_probability 0.15 \
    --per_gpu_train_batch_size 16 \
    --gradient_accumulation_steps 4 \
    --num_train_epochs 3 \
    --learning_rate 5e-5 \
    --save_steps 10000 \
    --logging_steps 500
```

Once training is complete, model checkpoints and related files will be saved in the directory specified by `--output_dir`.

## File Structure (Example)

```
.
├── split_chrom.py
├── Ladderpath_multiprocess_DNASequence.py
├── Merge_the_multiplicities.py
├── build_vocab.py
├── LPT_pretrain.py
├── README.md
│
├── data/
│   ├── hg19.fa
│   └── hg19_zhengfu_sorted_output_file_end.csv
│
├── hg19_gene_annotation.csv        # Output of Step 1
│
├── hg19_json_file/                   # Output of Step 2
│   ├── pom_data_batch_1.json
│   ├── pom_data_batch_2.json
│   └── ...
│
├── merge_multiplicities.json       # Output of Step 3
├── vocab.txt                       # Output of Step 4
│
├── tokenized_data.txt              # Input for Step 5 
├── bert_config.json                # Input for Step 5
│
└── lpt_model_output/                 # Output of Step 5
    ├── checkpoint-10000/
    ├── checkpoint-20000/
    └── ...
```
