from Bio import SeqIO
import pandas as pd


def parse_fasta(fasta_file):
    """ Parse the FASTA file and return a dictionary {chromosome: sequence} """
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = record.seq
    return sequences


def parse_csv(csv_file):
    """ Parse the CSV file and return a dictionary where keys are gene IDs and values are gene location information """
    genes = {}
    df = pd.read_csv(csv_file)

    for _, row in df.iterrows():
        gene_id = row['gene_id']
        seqid = row['seqid']
        start = row['start']
        end = row['end']
        strand = row['strand']
        # if strand == '+':
        # Process only records with a gene_id, and extract only the first row for the same gene_id
        genes[gene_id] = {
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand
            }

    return genes

def extract_gene_sequences(fasta_sequences, gff_genes):
    """ Extract DNA sequences of genes from the FASTA file based on location info, handling overlapping genes """
    gene_sequences = {}
    chromosome_boundaries = {}

    # Convert dictionary to list for indexing
    gene_list = list(gff_genes.items())

    for i, (gene_id, gene_info) in enumerate(gene_list):
        seqid = gene_info["seqid"]
        start = gene_info["start"]
        end = gene_info["end"]
        strand = gene_info["strand"]

        if seqid not in fasta_sequences:
            continue  # Skip genes if the chromosome is not in the FASTA file

        chromosome_seq = fasta_sequences[seqid]

        # Check if there is a next gene
        if i < len(gene_list) - 1:
            next_gene_id, next_gene_info = gene_list[i + 1]
            next_seqid = next_gene_info['seqid']

            # Check if on the same chromosome
            if seqid == next_seqid:
                # Check condition
                if start < next_gene_info['end']:
                    # Adjust start position
                    start = next_gene_info['end'] + 1

        # Ensure start and end positions are within valid range
        if start < 1:
            start = 1
        if end > len(chromosome_seq):
            end = len(chromosome_seq)

        # Extract gene sequence
        gene_seq = chromosome_seq[start - 1:end]  # GFF files are 1-based, Python is 0-based

        gene_sequences[f"{seqid}_{gene_id}"] = {
            "sequence": str(gene_seq),
            "strand": strand,
            "start": start,
            "end": end
        }

        # Record the start and end positions of the chromosome boundaries
        if seqid not in chromosome_boundaries:
            chromosome_boundaries[seqid] = {"start": start, "end": end}
        else:
            chromosome_boundaries[seqid]["start"] = min(chromosome_boundaries[seqid]["start"], start)
            chromosome_boundaries[seqid]["end"] = max(chromosome_boundaries[seqid]["end"], end)

    # Extract non-coding regions at the beginning and end of the chromosome
    for seqid, boundaries in chromosome_boundaries.items():
        chromosome_seq = fasta_sequences[seqid]
        # Beginning part
        if boundaries["start"] > 1:
            start_noncoding_seq = chromosome_seq[:boundaries["start"] - 1]
            gene_sequences[f"{seqid}_start_noncoding"] = {
                "sequence": str(start_noncoding_seq),
                "strand": "+",
                "start": 1,
                "end": boundaries["start"] - 1
            }
        # Ending part
        if boundaries["end"] < len(chromosome_seq):
            end_noncoding_seq = chromosome_seq[boundaries["end"]:]
            gene_sequences[f"{seqid}_end_noncoding"] = {
                "sequence": str(end_noncoding_seq),
                "strand": "+",
                "start": boundaries["end"] + 1,
                "end": len(chromosome_seq)
            }

    return gene_sequences

def extract_non_coding_regions_dict(fasta_sequences, gff_genes):
    """ Extract non-coding region sequences and store them in a dictionary """
    non_coding_dict = {}
    data = []
    for gene_id, info in gff_genes.items():
        data.append({
            "gene_id": gene_id,
            "seqid": info["seqid"],
            "start": info["start"],
            "end": info["end"],
            "strand": info["strand"]
        })

    # Create DataFrame and sort by chromosome (asc) and end position (desc)
    gene_df = pd.DataFrame(data)
    gene_df_sorted = gene_df.sort_values(['seqid', 'end'], ascending=[True, False])

    print("Sorted gene data:")
    print(gene_df_sorted.head())
    print(f"Total genes: {len(gene_df_sorted)}")

    for seqid, group in gene_df_sorted.groupby("seqid"):  # Group by chromosome
        if seqid not in fasta_sequences:
            print(f"Warning: Chromosome {seqid} not in reference sequence")
            continue

        chromosome_seq = fasta_sequences[seqid]
        previous_start = None  # Record the start position of the previous gene
        previous_gene_id = None  # Record the name of the previous gene

        # Reset counter
        non_coding_count = 0

        for _, row in group.iterrows():
            start, end = row["start"], row["end"]
            gene_id = row["gene_id"]
            strand = row["strand"]

            # Skip on the first iteration
            if previous_start is None:
                previous_start = start
                previous_gene_id = gene_id
                continue

            # Check for non-coding region
            if previous_start > end:
                non_coding_start = end + 1
                non_coding_end = previous_start - 1

                # Extra check for non-coding region length
                if non_coding_end >= non_coding_start:
                    try:
                        non_coding_seq = chromosome_seq[non_coding_start - 1:non_coding_end]

                        # Construct key name: Chromosome + Previous Gene Name + Non-coding Count + "_non_coding"
                        key = f"{seqid}_{previous_gene_id}_non_coding_{non_coding_count}"
                        non_coding_dict[key] = {
                            "sequence": str(non_coding_seq),
                            "strand": strand,
                            "start": non_coding_start,
                            "end": non_coding_end
                        }

                        # Print non-coding region information
                        print(f"Non-coding region: {key}")
                        print(f"  Length: {len(non_coding_seq)}")
                        print(f"  Start: {non_coding_start}")
                        print(f"  End: {non_coding_end}")

                        non_coding_count += 1
                    except IndexError:
                        print(f"Warning: Error extracting non-coding region {seqid}: {non_coding_start} - {non_coding_end}")

            # Update the start position and gene name for the next iteration
            previous_start = start
            previous_gene_id = gene_id

    print(f"Total non-coding regions: {len(non_coding_dict)}")
    return non_coding_dict



def save_to_csv(gene_sequences, output_file):
    """ Save gene IDs and sequences to a CSV file """
    data = []
    for gene_id, info in gene_sequences.items():
        data.append({
            "Gene_ID": gene_id,
            "DNA_Sequence": info["sequence"],
            "Strand": info["strand"],
            "Start": info["start"],
            "End": info["end"]
        })
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)


# Main program
def main():
    fasta_file = r"./hg19.fa"
    csv_file = r"./hg19_zhengfu_sorted_output_file_end.csv"
    output_file = r"./hg19_gene_annotation.csv"

    # Parse FASTA and CSV files
    fasta_sequences = parse_fasta(fasta_file)
    gff_genes = parse_csv(csv_file)
    print(gff_genes)

    # Extract gene sequences and non-coding regions
    gene_sequences = extract_gene_sequences(fasta_sequences, gff_genes)
    non_coding_regions = extract_non_coding_regions_dict(fasta_sequences, gff_genes)

    # Save results to a CSV file
    all_sequences = {**gene_sequences, **non_coding_regions}
    save_to_csv(all_sequences, output_file)


if __name__ == "__main__":
    main()