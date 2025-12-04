import os
import matplotlib.pyplot as plt
from collections import Counter
import os
import matplotlib.pyplot as plt
from collections import Counter

def sliding_windows(seq, win):
    return [seq[i:i+win] for i in range(len(seq) - win + 1)]

def read_fasta(path):
    seq = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()

def compute_cgsw_values(seq, win):
    counts = Counter(seq)
    total = counts["A"] + counts["T"] + counts["C"] + counts["G"]
    cg_total = counts["C"] + counts["G"]
    CGTOT = 100 * cg_total / total
    xs = []
    for w in sliding_windows(seq, win):
        c = w.count("C")
        g = w.count("G")
        cg_sw = c + g
        xs.append((CGTOT * cg_sw) / len(w))
    return xs

def kappa_ic_window(window_seq):
    A = window_seq
    N = len(A) - 1
    if N <= 0:
        return 0.0
    T = 0
    for u in range(1, N+1):
        B = A[u:]
        C = 0
        for i in range(len(B)):
            if A[i] == B[i]:
                C += 1
        T += (C / len(B)) * 100
    return round(T / N, 2)

def compute_kappa_values(seq, win):
    return [kappa_ic_window(w) for w in sliding_windows(seq, win)]

def promoter_pattern(seq, win):
    xs = compute_cgsw_values(seq, win)
    ys = compute_kappa_values(seq, win)
    n = min(len(xs), len(ys))
    return xs[:n], ys[:n]

def center_of_weight(xs, ys):
    n = min(len(xs), len(ys))
    return sum(xs[:n]) / n, sum(ys[:n]) / n

if __name__ == "__main__":


    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "FastaFiles")
    window_size = 30
    influenza_files = [f"flu{i}.fasta" for i in range(1, 11)]
    covid_files = [f"covid{i}.fasta" for i in range(1, 11)]

    centers = []

    plt.figure(figsize=(9, 7))
    for fname in influenza_files:
        seq = read_fasta(os.path.join(folder, fname))
        xs, ys = promoter_pattern(seq, window_size)
        plt.scatter(xs, ys, s=3, alpha=0.4, color="blue")
        cx, cy = center_of_weight(xs, ys)
        centers.append((cx, cy, fname.replace(".fasta", ""), "flu"))

    for fname in covid_files:
        seq = read_fasta(os.path.join(folder, fname))
        xs, ys = promoter_pattern(seq, window_size)
        plt.scatter(xs, ys, s=3, alpha=0.4, color="red")
        cx, cy = center_of_weight(xs, ys)
        centers.append((cx, cy, fname.replace(".fasta", ""), "covid"))

    plt.xlabel("(C+G)% (CGSW)")
    plt.ylabel("Kappa IC")
    plt.title("Objective Digital Strains for 10 Influenza + 10 COVID-19 Genomes")
    plt.grid(True)

    plt.figure(figsize=(9, 7))
    for cx, cy, label, group in centers:
        color = "blue" if group == "flu" else "red"
        plt.scatter(cx, cy, s=60, facecolors='none', edgecolors=color)
        plt.text(cx + 0.1, cy + 0.1, label, fontsize=9)

    plt.xlabel("Center (C+G)%")
    plt.ylabel("Center Kappa IC")
    plt.title("Centers of Weight of ODS (20 genomes)")
    plt.grid(True)

    plt.show()
