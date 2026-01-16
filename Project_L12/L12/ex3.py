import numpy as np
import json
import re

def preprocess_text(text):
    words = re.findall(r'\b\w+\b', text.lower())
    return words

def build_mapping(words):
    unique_words = sorted(list(set(words)))
    word_to_symbol = {}
    symbol_list = []
    
    for i, word in enumerate(unique_words):
        symbol = f"S{i}" 
        word_to_symbol[word] = symbol
        symbol_list.append(symbol)
        
    return word_to_symbol, symbol_list

def calculate_word_transition_matrix(word_sequence, word_to_symbol, symbol_list):
    n = len(symbol_list)
    symbol_to_index = {sym: i for i, sym in enumerate(symbol_list)}

    counts = np.zeros((n, n))
 
    for i in range(len(word_sequence) - 1):
        current_word = word_sequence[i]
        next_word = word_sequence[i+1]
        
        current_sym = word_to_symbol[current_word]
        next_sym = word_to_symbol[next_word]
        
        row_idx = symbol_to_index[current_sym]
        col_idx = symbol_to_index[next_sym]
        
        counts[row_idx][col_idx] += 1

    transition_matrix = np.zeros((n, n))
    for i in range(n):
        row_sum = np.sum(counts[i])
        if row_sum > 0:
            transition_matrix[i] = counts[i] / row_sum
            
    return transition_matrix.tolist()

def save_word_data_to_json_sparse(matrix, word_to_symbol, symbol_list, filename="word_transitions.json"):
    sparse_transitions = {}
    
    for i, row in enumerate(matrix):
        current_symbol = symbol_list[i]
        transitions = {}
        
        for j, prob in enumerate(row):
            if prob > 0:
                target_symbol = symbol_list[j]
                transitions[target_symbol] = round(prob, 4)
        
        if transitions:
            sparse_transitions[current_symbol] = transitions

    data = {
        "description": "Word transition matrix (Sparse Format).",
        "mapping_word_to_symbol": word_to_symbol,
        "transitions": sparse_transitions
    }
    
    with open(filename, "w") as f:
        json.dump(data, f, indent=4)

text_input = """
Kindergarten - I was five years old when I started kindergarten, but with a September birthday, I turned six soon after. Mrs. Baker was round and wore dresses every day. The first day of school, a little girl in my class cried and cried. She cried every day for weeks. I was curious about her. I watched her come in with her mother. Outside the door, she didn’t shed a tear. But once her mother guided her inside our classroom door, the water works began. She cried all morning. When she wasn’t shriek-crying, she was sobbing, her shoulders raising and falling dramatically. I don’t remember tears, but I do remember that she had a hyper-productive snot gland. She stood there helplessly while Mrs. Baker tried to calm her. Every morning it was the same. Every morning, Mrs. Baker assured her mother that Miss Cries-a-lot just needed a little more time. Months passed. Suddenly one day the classroom was quiet. She was gone and I knew exactly where she was before Mrs. Baker even made the announcement. She was next door, in another teacher’s room. You see, that little crybaby didn’t cry because she missed her mommy. She cried because she wanted to be in the other kinder class – the one with more toys. That’s what she howled about every single day. But that’s not the story Mrs. Baker told. During circle time that morning, she announced that the little girl was now in the other classroom because there were too many kids in her class. I suspect that every single one of my classmates had the same thought…”move me…move me.”

After circle time, all the kids were playing when I approached Mrs. Baker. “She really got to move because she cried, right?” I asked. She didn’t answer. “Go play,” she said.
She would later tell my mother how smart I was. “She figured out exactly what happened.” She shared.

It didn’t take a rocket scientist.

Now it's your turn. Open up a page on your word procesor and begin with Kindergarten. Write just a paragraph or two of the first memory that comes to mind. Tie it up quickly and stop.
Then challenge yourself to add to it regularly. As a writer, I find it is a good "warm-up" to the writing I must get done on a particular day. Or you may want to use it if you are stuck for what to write. And in the mean time you are doing something even more amazing, you are taking something that is yours and yours alone- a memory- and recording it for generations to come. Good for you! Write Now, Karen
"""

word_sequence = preprocess_text(text_input)
word_map, symbols = build_mapping(word_sequence)
matrix = calculate_word_transition_matrix(word_sequence, word_map, symbols)
save_word_data_to_json_sparse(matrix, word_map, symbols)