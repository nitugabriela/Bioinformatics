import json
import numpy as np
import random

def generate_text(json_file, word_count=50):
    with open(json_file, 'r') as f:
        data = json.load(f)

    word_to_sym = data['mapping_word_to_symbol']
    transitions = data['transitions']

    sym_to_word = {v: k for k, v in word_to_sym.items()}
    
    valid_starts = list(transitions.keys())
    current_symbol = random.choice(valid_starts)
    
    output = [sym_to_word[current_symbol]]
    
    for _ in range(word_count - 1):
        if current_symbol not in transitions:
            break
            
        next_options = transitions[current_symbol]
        candidates = list(next_options.keys())
        probs = list(next_options.values())

        total_prob = sum(probs)
        normalized_probs = [p / total_prob for p in probs]
        
        current_symbol = np.random.choice(candidates, p=normalized_probs)
        output.append(sym_to_word[current_symbol])
        
    return " ".join(output)

if __name__ == "__main__":
    print("--- Generated Text ---")
    print(generate_text("word_transitions.json", word_count=100))