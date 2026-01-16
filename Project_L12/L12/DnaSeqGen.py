import json
import numpy as np

def generate_dna(json_file, length=50):
    with open(json_file, 'r') as f:
        data = json.load(f)

    states = data['states']
    matrix = np.array(data['matrix'])

    current_state = np.random.choice(states)
    sequence = [current_state]
    
    state_to_index = {state: i for i, state in enumerate(states)}
    
    for _ in range(length - 1):
        curr_idx = state_to_index[current_state]
        probs = matrix[curr_idx]

        if np.sum(probs) > 0:
            probs = probs / np.sum(probs)
            current_state = np.random.choice(states, p=probs)
            sequence.append(current_state)
        else:
            break
        
    return "".join(sequence)

if __name__ == "__main__":
    print("--- Generated DNA Sequence ---")
    print(generate_dna("dna_transition_matrix.json", length=50))