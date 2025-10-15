import math


def calculate_tm_simple(dna_seq):
    dna_seq = dna_seq.upper()
    a_count = dna_seq.count('A')
    t_count = dna_seq.count('T')
    g_count = dna_seq.count('G')
    c_count = dna_seq.count('C')
    tm = 4 * (g_count + c_count) + 2 * (a_count + t_count)
    return tm

def calculate_tm_advanced(dna_seq, Na = 0.0001):
    dna_seq = dna_seq.upper()
    GC_percent = (dna_seq.count('G') + dna_seq.count('C')) / len(dna_seq) * 100
    tm = - (81.5 + 16.6 * math.log10(Na) + 0.41 * GC_percent - (600 / len(dna_seq)))
    return tm

while True:
    dna_seq = input("Enter a DNA sequence: ").strip()
    if all(c in 'ATGC' for c in dna_seq):
        break
    else:
        print("Invalid DNA sequence. Enter a sequence containing only A, T, G, and C.")

simple_tm = calculate_tm_simple(dna_seq)
advanced_tm = calculate_tm_advanced(dna_seq)

print(f"Simple Tm: {simple_tm} °C")
print(f"Advanced Tm: {advanced_tm:.2f} °C") 