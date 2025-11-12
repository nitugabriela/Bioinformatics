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

def plot_all_tandem_repeats(repeats, save_as_png=False):
    motif_counts = Counter([r["motif"] for r in repeats])
    sorted_items = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)
    motifs, counts = zip(*sorted_items)

    plt.figure(figsize=(10, max(6, len(motifs) * 0.3)))
    plt.barh(motifs, counts, color="mediumseagreen")
    plt.xlabel("Frequency")
    plt.ylabel("Motif (3–10 bases)")
    plt.title("All Tandem Repeats in DNA Sequence")
    plt.tight_layout()
    if save_as_png:
        plt.savefig("tandem_repeats.png", dpi=300)
    plt.show()

def main():
    fasta_file = "sequence.fasta"
    dna_seq = read_fasta(fasta_file)
    print(f"Sequence length: {len(dna_seq)} nucleotides")

    repeats = find_tandem_repeats(dna_seq, 3, 10)
    if repeats:
        print("\nDetected Tandem Repeats (3–10 bases):")
        for r in repeats:
            print(f"{r['motif']} -> count: {r['count']}, start: {r['start']}, end: {r['end']}")
        plot_all_tandem_repeats(repeats, save_as_png=False)
    else:
        print("No tandem repeats of length 3–10 found.")

if __name__ == "__main__":
    main()
