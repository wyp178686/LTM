import json
import os
import gc
import ladderpath as lp
import time
from multiprocessing import Process
import psutil
import pandas as pd

# Get current process
process = psutil.Process(os.getpid())

# Define set of valid bases
valid_bases = {"A", "T", "C", "G"}

def log_memory_usage(stage):
    """
    Print current memory usage to help identify memory spikes.
    """
    mem_info = process.memory_info()
    print(f"[{stage}] Current memory usage: {mem_info.rss / 1024 ** 2:.2f} MB")

def process_batch(batch_index, current_sequences, output_folder):
    """
    Process a single batch of data.
    """
    time.sleep(2)
    print(
        f"Processing batch {batch_index}, Total length: {sum(len(seq) for seq in current_sequences)}, List length: {len(current_sequences)}"
    )
    # Call ladderpath function
    lpjson = lp.get_ladderpath(
        current_sequences,
        info='V1.0.0.20240910_Alpha',
        estimate_eta=False,
        estimate_eta_para=[10, 'global'],
        save_file_name=f'./ladder_{batch_index}.json',
        show_version=True
    )

    # Check if the return value is valid
    if lpjson is None:
        print(f"Input for batch {batch_index} is invalid, skipping.")
    else:
        # Display the 3 Ladderpath indices
        index3 = lp.disp3index(lpjson)
        print(f"Ladderpath indices for batch {batch_index}: {index3}")

        # Obtain the Poset-Multiset (POM) representation of the Ladderpath
        pom, pom_str = lp.POM_from_JSON(lpjson, display_str=False)

        pom_file_path = f"{output_folder}/pom_data_batch_{batch_index}.json"
        with open(pom_file_path, 'w') as json_file:
            json.dump(pom, json_file, indent=4)
        print(f"POM dictionary for batch {batch_index} saved to {pom_file_path}")

    # Perform garbage collection after the subprocess completes the task
    print(f"Batch {batch_index} finished, starting garbage collection")
    gc.collect()  # Manually trigger garbage collection
    log_memory_usage(f"After processing batch {batch_index}")  # Memory checkpoint

# Define generator function
def sequence_generator(sequences, max_length):
    current_length = 0
    current_sequences = []
    for sequence in sequences:
        sequence = sequence.upper()  # Convert to uppercase
        # Split sequence: break at non-standard characters
        fragments = []
        fragment = []
        for char in sequence:
            if char in valid_bases:
                fragment.append(char)  # Build valid fragment
            else:
                if fragment:
                    fragments.append("".join(fragment))  # Save valid fragment
                    fragment = []  # Reset fragment
        if fragment:
            fragments.append("".join(fragment))  # Save the last valid fragment

        # Add split fragments to current_sequences
        for frag in fragments:
            frag_length = len(frag)

            # Skip if fragment length is <= 1
            if frag_length <= 1:
                continue

            # If fragment length exceeds max_length, split the fragment
            if frag_length > max_length:
                split_fragments = [frag[i:i + max_length] for i in range(0, len(frag), max_length)]
                for split_frag in split_fragments:
                    if current_length + len(split_frag) > max_length:
                        yield current_sequences
                        current_sequences = []
                        current_length = 0
                    current_sequences.append(split_frag)
                    current_length += len(split_frag)
            else:
                if current_length + frag_length > max_length:
                    yield current_sequences
                    current_sequences = []
                    current_length = 0
                current_sequences.append(frag)
                current_length += frag_length

    # Yield the last batch (if any remain)
    if current_sequences:
        yield current_sequences

# Main processing logic
if __name__ == "__main__":
    max_length = 1000000  # Max length limit
    output_folder = r'./hg19_json_file'
    file_path = r"./hg19_gene_annotation.csv"
    data = pd.read_csv(file_path)
    sequences = data['DNA_Sequence']

    num_workers = 50 # Control concurrency, adjust based on CPU and memory
    processes = []  # List to store processes
    batch_index = 1
    # Use generator to create arguments
    for current_sequences in sequence_generator(sequences, max_length):
        # Create and start each subprocess
        process = Process(target=process_batch, args=(batch_index, current_sequences, output_folder))
        process.start()
        processes.append(process)

        # Control concurrency, wait for some processes to finish
        if len(processes) >= num_workers:
            for p in processes:
                p.join()  # Wait for subprocess to finish
            processes = []  # Clear process list
        batch_index += 1

    # Wait for all processes to finish
    for p in processes:
        p.join()

    print("Main program execution completed")