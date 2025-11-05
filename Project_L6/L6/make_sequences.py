import random

def random_dna(length, gc_percent):
    gc = int(length * gc_percent / 100)
    at = length - gc
    seq = random.choices("GC", k=gc) + random.choices("AT", k=at)
    random.shuffle(seq)
    return "".join(seq)

def main():
    random.seed(42)
    configs = [
        ("seq1", 1800, 35),
        ("seq2", 2000, 45),
        ("seq3", 2200, 55),
        ("seq4", 2500, 60),
        ("seq5", 2700, 50),
    ]
    for name, length, gc in configs:
        s = random_dna(length, gc)
        with open(f"{name}.fasta", "w", encoding="utf-8") as f:
            f.write(f">{name}\n")
            for i in range(0, len(s), 80):
                f.write(s[i:i+80] + "\n")
        print(f"Wrote {name}.fasta ({length} bp, GC={gc}%)")

if __name__ == "__main__":
    main()
