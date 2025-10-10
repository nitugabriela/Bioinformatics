from itertools import product

S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

def generate_combinations(length):
    combinations = []
    for p in product('ACGT', repeat = length):
        combo = ""
        for char in p:
            combo += char
        combinations.append(combo)
    return combinations


def count_combinations(S, combinations, k):
    total_substr = len(S) - k + 1

    counts = {}
    for c in combinations:
        counts[c] = 0

    for i in range(total_substr):
        substring = S[i:i+k]
        if substring in counts:
            counts[substring] += 1  

    return counts, total_substr


def calculate_percentages(counts, total_substr):
    percentages = {}
    for c, count in counts.items():
        percentages[c] = round((count / total_substr) * 100, 2)
    return percentages


dinucleotides = generate_combinations(2)
di_counts, di_total = count_combinations(S, dinucleotides, 2)
di_percentages = calculate_percentages(di_counts, di_total)

print("Dinucleotide Percentages:")
for dinucleotide, percentage in di_percentages.items():
    print(f"{dinucleotide}: {percentage}%")

trinucleotides = generate_combinations(3)
tri_counts, tri_total = count_combinations(S, trinucleotides, 3)
tri_percentages = calculate_percentages(tri_counts, tri_total)

print("\nTrinucleotide Percentages:")
for trinucleotide, percentage in tri_percentages.items():
    print(f"{trinucleotide}: {percentage}%")

