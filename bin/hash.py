from collections import defaultdict
import sys

with open('data.txt', 'r') as file:
    data_dict=defaultdict(list)
    # Iterate through each line in the file
    for line in file:
        # Split the line into key and value
        key, value = line.strip().split()
        # Add the key-value pair to the dictionary
        data_dict[key].append(value)

# Print the resulting dictionary to file genes2Probes.txt
with open('genes2Probes.txt', 'w') as file:
    sys.stdout = file
    for key, value in data_dict.items():
        print(f"{key}: {value}")

# Reset stdout to default
    sys.stdout = sys.__stdout__

