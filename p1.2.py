S = "ACGGGCATATGCGC"

def percentages(seq:str) -> str:
    seq = seq.upper()

    freq = {}

    freq["A"] = 0
    freq["C"] = 0
    freq["G"] = 0
    freq["T"] = 0

    for i in seq:
        if i in freq:
            freq[i] += 1

    for key, value in freq.items():
        percentage = (value / len(seq)) * 100 if len(seq) else 0
        print(f"{key}: {percentage:.2f}%")

print("Sequence:", S)
percentages(S)
