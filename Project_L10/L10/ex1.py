import math
from typing import Dict, List, Tuple

MOTIFS = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT",
]

S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

BASES = "ACGT"


def make_count_matrix(motifs: List[str]) -> Dict[str, List[int]]:
    L = len(motifs[0])
    counts = {b: [0] * L for b in BASES}
    for m in motifs:
        if len(m) != L:
            raise ValueError("All motifs must have same length")
        for j, ch in enumerate(m):
            if ch not in BASES:
                raise ValueError(f"Invalid base {ch} in motif {m}")
            counts[ch][j] += 1
    return counts


def add_pseudocounts(counts: Dict[str, List[int]], alpha: int = 1) -> Dict[str, List[int]]:
    L = len(next(iter(counts.values())))
    weights = {b: [0] * L for b in BASES}
    for b in BASES:
        for j in range(L):
            weights[b][j] = counts[b][j] + alpha
    return weights


def make_frequency_matrix(weights: Dict[str, List[int]]) -> Dict[str, List[float]]:
    L = len(next(iter(weights.values())))
    freqs = {b: [0.0] * L for b in BASES}
    for j in range(L):
        col_sum = sum(weights[b][j] for b in BASES)
        for b in BASES:
            freqs[b][j] = weights[b][j] / col_sum
    return freqs


def make_log_likelihood_matrix(freqs: Dict[str, List[float]], null_prob: float = 0.25) -> Dict[str, List[float]]:
    L = len(next(iter(freqs.values())))
    llr = {b: [0.0] * L for b in BASES}
    for b in BASES:
        for j in range(L):
            p = freqs[b][j]
            llr[b][j] = math.log(p / null_prob)
    return llr


def score_window(window: str, llr: Dict[str, List[float]]) -> float:
    if len(window) != len(next(iter(llr.values()))):
        raise ValueError("Window length must match motif length")
    score = 0.0
    for j, ch in enumerate(window):
        if ch not in BASES:
            raise ValueError(f"Invalid base {ch} in window {window}")
        score += llr[ch][j]
    return score


def scan_sequence(seq: str, L: int, llr: Dict[str, List[float]]) -> List[Tuple[int, str, float]]:
    results = []
    for i in range(len(seq) - L + 1):
        w = seq[i:i+L]
        results.append((i, w, score_window(w, llr)))
    return results


def print_matrix(title: str, mat: Dict[str, List], fmt: str = "{:>8}") -> None:
    L = len(next(iter(mat.values())))
    print("\n" + title)
    header = "     " + "".join([fmt.format(j+1) for j in range(L)])
    print(header)
    for b in BASES:
        row = f"{b:>3}  " + "".join(fmt.format(v if not isinstance(v, float) else float(v)) for v in mat[b])
        print(row)


def main():
    L = len(MOTIFS[0])

    counts = make_count_matrix(MOTIFS)
    weights = add_pseudocounts(counts, alpha=1)
    freqs = make_frequency_matrix(weights)
    llr = make_log_likelihood_matrix(freqs, null_prob=0.25)

    print_matrix("1) Count matrix C", counts, fmt="{:>8}")
    print_matrix("2) Weight matrix W = C + 1", weights, fmt="{:>8}")

    print("\n3) Relative frequencies P (with pseudocounts)")
    print("     " + "".join([f"{j+1:>10}" for j in range(L)]))
    for b in BASES:
        print(f"{b:>3}  " + "".join([f"{freqs[b][j]:>10.4f}" for j in range(L)]))

    print("\n4) Log-likelihoods LLR = ln(P/0.25)")
    print("     " + "".join([f"{j+1:>10}" for j in range(L)]))
    for b in BASES:
        print(f"{b:>3}  " + "".join([f"{llr[b][j]:>10.4f}" for j in range(L)]))

    results = scan_sequence(S, L, llr)

    print("\n5) Sliding window scores over S")
    print(" idx  window       score")
    for i, w, sc in results:
        print(f"{i:>4}  {w}  {sc:>9.4f}")

    best = max(results, key=lambda x: x[2])
    top5 = sorted(results, key=lambda x: x[2], reverse=True)[:5]

    print("\nBest hit:")
    print(f" start={best[0]} window={best[1]} score={best[2]:.4f}")

    print("\nTop 5 hits:")
    for i, w, sc in top5:
        print(f" start={i:>2}  {w}  {sc:>9.4f}")

    scores = sorted(sc for _, _, sc in results)
    median = scores[len(scores)//2]
    print(f"\nMedian score: {median:.4f}")


if __name__ == "__main__":
    main()
