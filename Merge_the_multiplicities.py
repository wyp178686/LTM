import json
import os
from collections import defaultdict

# Define file paths
file_paths = [
    './hg19_json_file_1000000',
]

# Collect all JSON file paths
files = []
for file_path in file_paths:
    files.extend([os.path.join(file_path, f) for f in os.listdir(file_path) if f.endswith('.json')])

# Use defaultdict to accumulate multiplicities
token_counts = defaultdict(int)

# Read files and accumulate multiplicities
for file in files:
    with open(file, 'r', encoding='utf-8') as f:
        data = json.load(f)
        for layer in data.values():  # Iterate through each layer
            for token, count in layer.items():
                # Only count if the key length is greater than or equal to 1
                if len(token) >= 1:
                    token_counts[token] += count  # Accumulate multiplicity

# Print dictionary information
print(f'Length of dictionary: {len(token_counts)}')

# Calculate the cumulative length of all keys
total_length = sum(len(key) for key in token_counts)
print("Cumulative length of all keys in dictionary:", total_length)

# Write results to a JSON file
output_file = './merge_multiplicities.json'  # Output filename
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(token_counts, f, ensure_ascii=False, indent=4)