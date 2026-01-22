import math
import pandas as pd
import matplotlib.pyplot as plt
import re

text_eminescu = """
A fost odată ca-n povești, a fost ca niciodată,
Din rude mari împărătești, o prea frumoasă fată.
Și era una la părinți și mândră-n toate cele,
Cum e Fecioara între sfinți și luna între stele.
Porni luceafărul de sus, mai stins cu o aripă,
Și-n locul lui a fost apus, în acea scurtă clipă.
Cobori în jos, luceafăr blând, alunecând pe o rază,
Pătrunde-n casă și în gând și viața-mi luminează!
"""

text_stanescu = """
Ea era frumoasă ca umbra unei idei,
a pielei de copil mirosea spinarea ei,
a piatră proaspăt spartă,
a strigăt într-o limbă moartă.
Nu avea greutate, ca respirarea.
Râzândă și plângândă cu lacrimi mari
era sărată ca sarea slăvită la ospețe de barbari.
Spune-mi, dacă te-aș prinde-ntr-o zi
și ți-aș săruta talpa piciorului,
nu-i așa că ai șchiopăta puțin,
de teamă să nu-mi strivești sărutul?
"""

suspect_text = """
A fost odată ca-n povești o piatră proaspăt spartă.
Din rude mari împărătești mirosea spinarea ei.
Cobori în jos, luceafăr blând, alunecând pe o rază,
dar nu avea greutate, ca respirarea.
Și era una la părinți și mândră-n toate cele,
dar dacă te-aș prinde-ntr-o zi și ți-aș săruta talpa piciorului.
Pătrunde-n casă și în gând și viața-mi luminează,
de teamă să nu-mi strivești sărutul.
"""

def tokenize(text):
    text = text.lower()
    words = re.findall(r'\w+', text)
    return words

def train_markov_model(text, vocabulary):
    words = tokenize(text)
    model = {w1: {w2: 1 for w2 in vocabulary} for w1 in vocabulary}
    
    for i in range(len(words) - 1):
        curr_w = words[i]
        next_w = words[i+1]
        if curr_w in model and next_w in model[curr_w]:
            model[curr_w][next_w] += 10
            
    prob_matrix = {}
    for w1, transitions in model.items():
        total = sum(transitions.values())
        prob_matrix[w1] = {w2: count/total for w2, count in transitions.items()}
        
    return prob_matrix

all_text = text_eminescu + text_stanescu + suspect_text
vocab = sorted(list(set(tokenize(all_text))))

prob_eminescu = train_markov_model(text_eminescu, vocab)
prob_stanescu = train_markov_model(text_stanescu, vocab)

def get_log_likelihood(w1, w2):
    p_e = prob_eminescu.get(w1, {}).get(w2, 1/len(vocab))
    p_s = prob_stanescu.get(w1, {}).get(w2, 1/len(vocab))
    return math.log2(p_e / p_s)

suspect_words = tokenize(suspect_text)
window_size = 5
results = []
plot_indices = []

print(f"{'WINDOW TEXT':<40} | {'SCORE':<10} | {'VERDICT'}")
print("-" * 70)

for i in range(len(suspect_words) - window_size + 1):
    window = suspect_words[i : i + window_size]
    window_score = 0
    
    for j in range(len(window) - 1):
        window_score += get_log_likelihood(window[j], window[j+1])
    
    if window_score > 0.5: verdict = "EMINESCU"
    elif window_score < -0.5: verdict = "STANESCU"
    else: verdict = "UNCERTAIN"
    
    results.append(window_score)
    plot_indices.append(i)
    
    snippet = " ".join(window)
    print(f"{snippet:<40} | {window_score:6.2f}     | {verdict}")

plt.figure(figsize=(12, 6))
colors = ['green' if s > 0 else 'red' for s in results]
plt.bar(plot_indices, results, color=colors, alpha=0.7)
plt.axhline(0, color='black', linewidth=1)
plt.title("Plagiarism Detection Analysis")
plt.xlabel("Window Position")
plt.ylabel("Log-Likelihood Score")
plt.grid(True, alpha=0.3)
plt.show()