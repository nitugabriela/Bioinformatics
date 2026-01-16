import numpy as np

def run_prediction(matrix, start_vector, steps=5):

    current_state = np.array(start_vector)
    transition_matrix = np.array(matrix)
    
    print(f"Initial State (Step 0): {current_state}")
    print("-" * 40)

    for step in range(1, steps + 1):

        current_state = np.dot(current_state, transition_matrix)

        formatted_state = np.round(current_state, 4)
        print(f"Step {step}: {formatted_state}")

M = [
    [0.3, 0.35, 0.35],
    [0.0, 0.0,  1.0],
    [0.9, 0.0,  0.1]
]

v = [1, 0, 0]

if __name__ == "__main__":
    run_prediction(M, v, steps=5)