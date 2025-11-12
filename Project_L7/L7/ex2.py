import os
from collections import Counter
import matplotlib.pyplot as plt

def read_fasta(filename):
    sequence = []
    with open(filename, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence.append(line.strip())
    return "".join(sequence)

def find_tandem_repeats(sequence, min_len=3, max_len=10):
    results = []
    seq_len = len(sequence)
    for length in range(min_len, max_len + 1):
        i = 0
        while i <= seq_len - length:
            motif = sequence[i:i+length]
            count = 1
            j = i + length
            while sequence[j:j+length] == motif:
                count += 1
                j += length
            if count > 1:
                results.append({
                    "motif": motif,
                    "count": count,
                    "start": i,
                    "end": j
                })
                i = j
            else:
                i += 1
    return results

def plot_tandem_repeats(repeats, genome_name):
    motif_counts = Counter([r["motif"] for r in repeats])
    if not motif_counts:
        print(f"No tandem repeats found in {genome_name}")
        return

    sorted_items = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)
    motifs, counts = zip(*sorted_items)

    plt.figure(figsize=(10, max(6, len(motifs) * 0.3)))
    plt.barh(motifs, counts, color="mediumseagreen")
    plt.xlabel("Frequency")
    plt.ylabel("Motif (3–10 bases)")
    plt.title(f"Tandem Repeats in {genome_name}")
    plt.tight_layout()
    plt.savefig(f"{genome_name}_repeats.png", dpi=300)
    plt.close()

def process_all_genomes():
    for i in range(1, 11):
        filename = f"flu{i}.fna"
        if not os.path.exists(filename):
            print(f"{filename} not found.")
            continue
        dna_seq = read_fasta(filename)
        repeats = find_tandem_repeats(dna_seq, 3, 10)
        plot_tandem_repeats(repeats, f"flu{i}")

if __name__ == "__main__":
    process_all_genomes()
