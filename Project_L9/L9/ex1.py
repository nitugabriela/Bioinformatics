import random

MIN_LEN = 200
MAX_LEN = 400
NUM_TES = 4
WRITE_FASTA = False

def random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))

def insert_subsequence(seq, sub, pos):
    return seq[:pos] + sub + seq[pos:]

TEs = {
    "TE1": "TTAAGGCTACCTTAAGGCTACC",
    "TE2": "CGCGTATGATCGCGTATGATC",
    "TE3": "AAGCTTAGGCTAAGCTTAGGC",
    "TE4": "GTCAGTACCGATGTCAGTACCGA"
}

TE_NAMES = list(TEs.keys())[:NUM_TES]

target_len = random.randint(MIN_LEN, MAX_LEN)
background_len = target_len - sum(len(TEs[name]) for name in TE_NAMES)
background = random_dna(background_len)

insertion_positions = sorted(random.sample(range(len(background) + 1), k=len(TE_NAMES)))

sequence = background
te_locations = []
offset = 0

for pos, te_name in zip(insertion_positions, TE_NAMES):
    te_seq = TEs[te_name]
    real_pos = pos + offset
    sequence = insert_subsequence(sequence, te_seq, real_pos)
    start = real_pos + 1
    end = real_pos + len(te_seq)
    te_locations.append((te_name, start, end))
    offset += len(te_seq)

print(f"Final sequence length: {len(sequence)} bp")
print("\nTransposable elements locations (1-based):")
for name, start, end in te_locations:
    print(f"{name}: {start}-{end} ({end - start + 1} bp)")

print("\nFASTA format:\n")
header = f">artificial_seq_with_{len(TE_NAMES)}TEs_len{len(sequence)}"
print(header)
for i in range(0, len(sequence), 60):
    print(sequence[i:i+60])

if WRITE_FASTA:
    filename = "artificial_seq_with_TEs.fasta"
    with open(filename, "w") as f:
        f.write(header + "\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")
