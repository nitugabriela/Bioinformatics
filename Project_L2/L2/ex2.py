def existing_combinations(S, k):
    found = []
    for i in range(len(S) - k + 1):
        substring = S[i:i+k]
        if substring not in found:
            found.append(substring)
    return found

S = "ABBA"

dinucleotides = existing_combinations(S, 2)
trinucleotides = existing_combinations(S, 3)

print("Dinucleotides:", dinucleotides)
print("Trinucleotides:", trinucleotides)