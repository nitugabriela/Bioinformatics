def find_exact(sequence, pattern):
    positions = []
    i = 0
    while True:
        i = sequence.find(pattern, i)
        if i == -1:
            break
        positions.append((i, i + len(pattern)))
        i += 1
    return positions

def reverse_complement(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]

def find_inverted_repeats(sequence, min_len=6, max_len=20):
    results = []
    n = len(sequence)
    for i in range(n):
        for j in range(i + min_len, min(i + max_len, n)):
            left = sequence[i:j]
            right = reverse_complement(left)
            k = sequence.find(right, j)
            if k != -1:
                results.append((i, j, k, k + len(right)))
    return results

def find_direct_repeats(sequence, min_len=6, max_len=20):
    results = []
    n = len(sequence)
    for i in range(n):
        for j in range(i + min_len, min(i + max_len, n)):
            motif = sequence[i:j]
            k = sequence.find(motif, j)
            if k != -1:
                results.append((i, j, k, k + len(motif)))
    return results

sequence = (
    "AACCCACCACAACTGTCCCACTACAAAGCTTGTTTAAGGCTACCTTAAGGCTACCTTAAG"
    "AGAGCCATATGCGCTGTGTATGGAAGGGCCTGGTGCCCGAGTGGAATTGGTAGCTGGCCT"
    "CTATCAGTCGGGAAGAGGCCGCCCGCCTTCCCGACCTAGGCGCACCACCGCAGGCAACTT"
    "CCGTCATAGGTTCCCTGATCGGCAAGCGCGTATGATCGCGTATGATCACTTCTCTTCGGC"
    "CGAAGCTTAGGCTAAGCTTAGGCACCTATCCTATGGGACGTGATTGGCGTTATCACAGAT"
    "CAAGGCCAGCATCTGACCCCCCGTCTTAGAGTGAAGTTCGGCACTCATTTTCACGTGGGG"
    "TGGATCGTCAGTACCGATGTCAGTACCGATG"
)

TEs = {
    "TE1": "TTAAGGCTACCTTAAGGCTACC",
    "TE2": "CGCGTATGATCGCGTATGATC",
    "TE3": "AAGCTTAGGCTAAGCTTAGGC",
    "TE4": "GTCAGTACCGATGTCAGTACCGA"
}

print("Transposable elements found:\n")
for name, te_seq in TEs.items():
    hits = find_exact(sequence, te_seq)
    for start, end in hits:
        start_1 = start + 1
        end_inc = end
        length = end - start
        print(f"{name}: start = {start_1}, end = {end_inc}, length = {length} bp")

irs = find_inverted_repeats(sequence)
print("\nInverted repeats (left IR → right IR):\n")
if not irs:
    print("None found")
else:
    for idx, (ls, le, rs, re) in enumerate(irs, 1):
        left_start = ls + 1
        left_end = le
        right_start = rs + 1
        right_end = re
        length = le - ls
        print(f"IR {idx}: left = {left_start}-{left_end}, right = {right_start}-{right_end}, length = {length} bp")

drs = find_direct_repeats(sequence)
print("\nDirect repeats (left DR → right DR):\n")
if not drs:
    print("None found")
else:
    for idx, (ls, le, rs, re) in enumerate(drs, 1):
        left_start = ls + 1
        left_end = le
        right_start = rs + 1
        right_end = re
        length = le - ls
        print(f"DR {idx}: left = {left_start}-{left_end}, right = {right_start}-{right_end}, length = {length} bp")
