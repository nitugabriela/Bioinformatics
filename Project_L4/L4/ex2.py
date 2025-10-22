from collections import Counter
import matplotlib.pyplot as plt

genetic_code = {
    'UUU': 'Phe', 'UUC': 'Phe',
    'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
    'AUG': 'Met',  
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'UAU': 'Tyr', 'UAC': 'Tyr',
    'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop',
    'CAU': 'His', 'CAC': 'His',
    'CAA': 'Gln', 'CAG': 'Gln',
    'AAU': 'Asn', 'AAC': 'Asn',
    'AAA': 'Lys', 'AAG': 'Lys',
    'GAU': 'Asp', 'GAC': 'Asp',
    'GAA': 'Glu', 'GAG': 'Glu',
    'UGU': 'Cys', 'UGC': 'Cys',
    'UGG': 'Trp',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AGU': 'Ser', 'AGC': 'Ser',
    'AGA': 'Arg', 'AGG': 'Arg',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

def read_fasta(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if not line.startswith('>')]
    seq = ''.join(lines).upper().replace('T', 'U')
    return seq

def count_codons_and_amino_acids(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]
    codon_counts = Counter(codons)
    amino_acid_counts = Counter(genetic_code.get(c, '') for c in codons if genetic_code.get(c, '') not in ('', 'Stop'))
    return codon_counts, amino_acid_counts

def plot_top10_codons(codon_counts, title):
    top10 = codon_counts.most_common(10)
    codons, counts = zip(*top10)
    plt.bar(codons, counts)
    plt.title(title)
    plt.xlabel("Codon")
    plt.ylabel("Frequency")
    plt.xticks(rotation=45)

if __name__ == "__main__":
    covid_seq = read_fasta("covid19.fasta")
    flu_seq = read_fasta("influenza.fasta")

    covid_codons, covid_amino = count_codons_and_amino_acids(covid_seq)
    flu_codons, flu_amino = count_codons_and_amino_acids(flu_seq)

    # a) and b) plot top 10
    plt.figure(figsize=(12,5))
    plt.subplot(1,2,1)
    plot_top10_codons(covid_codons, "Top 10 Codons - SARS-CoV-2 Genome")
    plt.subplot(1,2,2)
    plot_top10_codons(flu_codons, "Top 10 Codons - Influenza A Genome")
    plt.tight_layout()
    plt.show()

    # c) compare top codons
    covid_top10 = {codon for codon, _ in covid_codons.most_common(10)}
    flu_top10 = {codon for codon, _ in flu_codons.most_common(10)}
    common = covid_top10.intersection(flu_top10)
    print("\nCommon top codons between SARS-CoV-2 and Influenza A:", ', '.join(common) or "None")

    # d) top 3 amino acids
    print("\nTop 3 amino acids in SARS-CoV-2:")
    for aa, count in covid_amino.most_common(3):
        print(f"  {aa}: {count}")

    print("\nTop 3 amino acids in Influenza A:")
    for aa, count in flu_amino.most_common(3):
        print(f"  {aa}: {count}")
