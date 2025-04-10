import re
import numpy as np

def extract_energies_and_n_values(file_path):

    with open(file_path, 'r') as file:
        input_data = file.readlines()

    energies = []
    n_values = []

    pattern = r"E= ([\d\.]+)\s+\[([\d,\s]+)\]"

    for line in input_data:
        match = re.search(pattern, line)
        if match:
            energy = float(match.group(1))
            n_list = list(map(int, match.group(2).split(',')))
            energies.append(energy)
            n_values.append(n_list)

    energies = np.array(energies)
    n_values = np.array(n_values)
    return energies, n_values