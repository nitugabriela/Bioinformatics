import os
import time
import subprocess
from statistics import mean, median, pstdev

REPEATS = 3
EXCLUDE = {"sequence.fasta", "reads.fasta", "reconstructed_contigs.fasta"}

def read_fasta_single(path):
    seq = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    s = "".join(seq).upper()
    return "".join(c for c in s if c in "ACGT")

def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return 100.0 * (g + c) / len(seq) if seq else 0.0

def write_sequence_fasta(seq, name="seq"):
    with open("sequence.fasta", "w", encoding="utf-8") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

def time_ex1():
    t = []
    for _ in range(REPEATS):
        start = time.time()
        subprocess.run(["python3", "ex1.py"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        t.append(time.time() - start)
    return median(t)

def pearson_r(xs, ys):
    if len(xs) < 2 or len(ys) < 2:
        return float("nan")
    mx, my = mean(xs), mean(ys)
    sx, sy = pstdev(xs), pstdev(ys)
    if sx == 0 or sy == 0:
        return float("nan")
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / len(xs)
    return cov / (sx * sy)

# make_sequences.py
import random
def random_dna(length, gc_percent):
    gc = int(length * gc_percent / 100)
    at = length - gc
    seq = random.choices("GC", k=gc) + random.choices("AT", k=at)
    random.shuffle(seq)
    return "".join(seq)

for name, gc in [("seq1", 30), ("seq2", 50), ("seq3", 70)]:
    s = random_dna(2000, gc)
    with open(f"{name}.fasta", "w") as f:
        f.write(f">{name}\n{s}\n")

def main():
    fasta_files = [f for f in os.listdir(".") if f.endswith(".fasta") and f not in EXCLUDE]
    fasta_files.sort()
    if not fasta_files:
        print("No test FASTA files found (e.g., seq1.fasta, seq2.fasta). Place them next to ex1.py and re-run.")
        return

    rows = []
    print(f"{'Sequence':<20}{'Length':>8}  {'GC%':>7}  {'Time(s)':>8}")
    print("-" * 48)

    for fname in fasta_files:
        seq = read_fasta_single(fname)
        if not seq:
            continue
        gc = gc_content(seq)
        write_sequence_fasta(seq, name=os.path.splitext(fname)[0])
        t = time_ex1()
        rows.append((fname, len(seq), gc, t))
        print(f"{fname:<20}{len(seq):>8}  {gc:>7.2f}  {t:>8.3f}")

    if rows:
        gc_vals = [r[2] for r in rows]
        time_vals = [r[3] for r in rows]
        r = pearson_r(gc_vals, time_vals)
        print("\nPearson correlation r(GC%, Time) =", f"{r:.4f}")

        with open("gc_time_results.csv", "w", encoding="utf-8") as f:
            f.write("sequence,length,gc_percent,time_seconds\n")
            for fname, L, gc, t in rows:
                f.write(f"{fname},{L},{gc:.4f},{t:.6f}\n")
        print('Results saved to "gc_time_results.csv".')

if __name__ == "__main__":
    main()



"""
We measured assembly time for three 2,000 bp sequences with GC contents of 30%, 50%, and 70%. 
Measured times were 0.188 s, 0.176 s, and 0.160 s, respectively. 
The Pearson correlation between GC% and time was r = −0.9975; however, this analysis is based on only three samples and very small absolute timing differences (~0.03 s). 
Given the algorithm’s complexity (dominated by the number of k-mers for fixed-length, error-free reads), we conclude there is no meaningful effect of GC% on runtime for this assembler. 
To obtain a robust estimate, we recommend testing 10–20 sequences across a broad GC range, increasing timing repetitions, and fixing K across runs.
"""