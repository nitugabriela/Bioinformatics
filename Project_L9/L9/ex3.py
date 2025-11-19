def reverse_complement(seq):
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def read_fasta(path):
    header = None
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:]
            else:
                seq_parts.append(line)
    return "".join(seq_parts), (header if header is not None else path)


def find_inverted_repeats_sliding(seq, offset=0,
                                  min_len=4, max_len=6, max_dist=2000):
    seq = seq.upper()
    n = len(seq)
    results = []

    for i in range(n):
        for L in range(min_len, max_len + 1):
            if i + L > n:
                continue
            left = seq[i:i + L]
            rc = reverse_complement(left)

            search_start = i + L
            search_end = min(i + L + max_dist, n - L + 1)

            for j in range(search_start, search_end):
                if seq[j:j + L] == rc:
                    ls = offset + i
                    le = offset + i + L
                    rs = offset + j
                    re = offset + j + L
                    results.append((L, ls, le, rs, re))

    return results


genome_files = [
    ("Escherichia", "Escherichia.fasta"),
    ("Pseudomonas", "Pseudomonas.fasta"),
    ("Bacillus",     "Bacillus.fasta")  
]

MAX_BASES_PER_GENOME = 100000   # analizează doar primele 200k baze
MAX_TO_SHOW = 30
MAX_DIST = 1000

print("Script started.\n")

for label, path in genome_files:
    try:
        genome_seq, genome_id = read_fasta(path)
    except FileNotFoundError:
        print(f"===== Genome {label} ({path}) =====")
        print("File not found, skipping.\n")
        continue

    print(f"===== Genome {label} ({path}) =====")
    print(f"Header: {genome_id}")
    print(f"Full length: {len(genome_seq)} bp")

    subseq = genome_seq[:MAX_BASES_PER_GENOME]
    print(f"Using first {len(subseq)} bp for IR search.\n")

    irs = find_inverted_repeats_sliding(
        subseq,
        offset=0,
        min_len=4,
        max_len=6,
        max_dist=MAX_DIST
    )

    print(f"Total inverted repeat pairs "
          f"(len 4–6, distance ≤ {MAX_DIST} bp): {len(irs)}\n")

    for idx, (L, ls, le, rs, re) in enumerate(irs[:MAX_TO_SHOW], 1):
        left_start_1 = ls + 1
        left_end_1 = le
        right_start_1 = rs + 1
        right_end_1 = re
        candidate_start = left_start_1
        candidate_end = right_end_1
        print(
            f"IR {idx}: length = {L}, "
            f"left IR = {left_start_1}-{left_end_1}, "
            f"right IR = {right_start_1}-{right_end_1}, "
            f"candidate transposon = {candidate_start}-{candidate_end}"
        )

    if len(irs) > MAX_TO_SHOW:
        print(f"... ({len(irs) - MAX_TO_SHOW} more IR pairs not shown)\n")

    print("-" * 60 + "\n")

print("Done.")
