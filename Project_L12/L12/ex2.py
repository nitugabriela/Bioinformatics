import numpy as np
import json
import random

def generate_dna_sequence(length=50):
    return "".join(random.choices("ACGT", k=length))

def calculate_transition_matrix(sequence, states):
    n = len(states)
    state_to_index = {state: i for i, state in enumerate(states)}

    counts = np.zeros((n, n))

    for i in range(len(sequence) - 1):
        current_s = sequence[i]
        next_s = sequence[i+1]
        
        row_idx = state_to_index[current_s]
        col_idx = state_to_index[next_s]
        
        counts[row_idx][col_idx] += 1

    transition_matrix = np.zeros((n, n))
    for i in range(n):
        row_sum = np.sum(counts[i])
        if row_sum > 0:
            transition_matrix[i] = counts[i] / row_sum
        else:
            pass

    return transition_matrix.tolist() 

def save_to_json(matrix, states, filename="dna_transition_matrix.json"):
    data = {
        "description": "Transition matrix derived from DNA sequence",
        "states": states,
        "matrix": matrix
    }
    
    with open(filename, "w") as f:
        json.dump(data, f, indent=4)
    print(f"Successfully saved matrix to {filename}")

dna_states = ['A', 'C', 'G', 'T']

sequence = generate_dna_sequence(50)
print(f"Analyzed Sequence (50 chars):\n{sequence}\n")

matrix = calculate_transition_matrix(sequence, dna_states)

save_to_json(matrix, dna_states)

print("Calculated Transition Matrix:")
print(f"      {'     '.join(dna_states)}")
for i, row in enumerate(matrix):
    print(f"{dna_states[i]}: {np.round(row, 2)}")