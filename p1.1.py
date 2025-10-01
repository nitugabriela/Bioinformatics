def detectSeq(seq: str) -> str:
    dna = set("ACGT")
    rna = set("ACGU")
    protein = set("ACDEFGHIKLMNPQRSTVWY")

    letters = set(seq.upper())

    if letters.issubset(dna):
        return "DNA"
    elif letters.issubset(rna):
        return "RNA"
    elif letters.issubset(protein):
        return "Protein"
    else:
        return "Unknown"
    
seq = input("Introduce the sequence: ")
print("Type:", detectSeq(seq))
