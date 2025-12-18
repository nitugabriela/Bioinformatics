import math
import os
from typing import Dict, List, Tuple, Optional

import matplotlib.pyplot as plt


# === Motifs (from your first exercise image) ===
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

BASES = "ACGT"


def read_fasta(path: str) -> str:
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line.upper())
    return "".join(seq_parts)


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
    return {b: [counts[b][j] + alpha for j in range(L)] for b in BASES}


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
            llr[b][j] = math.log(freqs[b][j] / null_prob)
    return llr


def score_window_llr(window: str, llr: Dict[str, List[float]]) -> Optional[float]:
    L = len(next(iter(llr.values())))
    if len(window) != L:
        raise ValueError("Window length must match motif length")
    s = 0.0
    for j, ch in enumerate(window):
        if ch not in BASES:
            return None
        s += llr[ch][j]
    return s


def scan_sequence_llr(seq: str, L: int, llr: Dict[str, List[float]]) -> List[float]:
    scores = []
    for i in range(len(seq) - L + 1):
        w = seq[i:i+L]
        sc = score_window_llr(w, llr)
        scores.append(float("nan") if sc is None else sc)
    return scores


def llr_to_lr(llr_scores: List[float]) -> List[float]:
    # LR = exp(LLR)  (clip extreme values to avoid overflow)
    out = []
    for x in llr_scores:
        if isinstance(x, float) and math.isnan(x):
            out.append(float("nan"))
        else:
            x_clip = max(min(x, 50.0), -50.0)
            out.append(math.exp(x_clip))
    return out


def normalize_0_1(values: List[float]) -> List[float]:
    valid = [v for v in values if not (isinstance(v, float) and math.isnan(v))]
    if not valid:
        return values
    mn, mx = min(valid), max(valid)
    if mx == mn:
        return [0.0 if not (isinstance(v, float) and math.isnan(v)) else float("nan") for v in values]
    return [((v - mn) / (mx - mn)) if not (isinstance(v, float) and math.isnan(v)) else float("nan") for v in values]


def percentile(values: List[float], p: float) -> float:
    valid = sorted(v for v in values if not (isinstance(v, float) and math.isnan(v)))
    if not valid:
        return float("nan")
    k = int(round((p / 100.0) * (len(valid) - 1)))
    return valid[k]


def top_k_hits(seq: str, llr_scores: List[float], L: int, k: int = 10) -> List[Tuple[int, str, float]]:
    hits = []
    for i, sc in enumerate(llr_scores):
        if not (isinstance(sc, float) and math.isnan(sc)):
            hits.append((i, seq[i:i+L], sc))
    hits.sort(key=lambda x: x[2], reverse=True)
    return hits[:k]


def plot_signal(genome_name: str, signal: List[float], hits: List[Tuple[int, str, float]], out_png: str) -> None:
    x = list(range(len(signal)))

    thr = percentile(signal, 95.0)  # threshold line like "background vs peak"
    plt.figure(figsize=(12, 4))
    plt.plot(x, signal)
    plt.title(f"{genome_name}: motif signal P along genome length L")
    plt.xlabel("L (window start position, 0-based)")
    plt.ylabel("P (normalized motif signal)")

    if not (isinstance(thr, float) and math.isnan(thr)):
        plt.axhline(thr, linestyle="--")
        plt.text(0, thr, "  95% threshold", va="bottom")

    # Mark top 5 most likely locations (peaks)
    for rank, (pos, _, _) in enumerate(hits[:5], start=1):
        plt.axvline(pos, linestyle="--")
        plt.text(pos, plt.ylim()[1], f"#{rank}", rotation=90, va="top", ha="right")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main():
    # Build motif model once
    counts = make_count_matrix(MOTIFS)
    weights = add_pseudocounts(counts, alpha=1)
    freqs = make_frequency_matrix(weights)
    llr = make_log_likelihood_matrix(freqs, null_prob=0.25)

    L = len(MOTIFS[0])

    script_dir = os.path.dirname(os.path.abspath(__file__))

    for idx in range(1, 11):
        fasta = os.path.join(script_dir, f"flu{idx}.fna")
        if not os.path.exists(fasta):
            print(f"[SKIP] Missing {fasta}")
            continue

        seq = read_fasta(fasta)
        if len(seq) < L:
            print(f"[SKIP] {fasta} too short")
            continue

        llr_scores = scan_sequence_llr(seq, L, llr)

        # Convert to probability-like signal:
        # LR = exp(LLR) then normalize to [0..1] so it matches the professor sketch
        lr_scores = llr_to_lr(llr_scores)
        P_signal = normalize_0_1(lr_scores)

        hits = top_k_hits(seq, llr_scores, L, k=10)

        print(f"\n=== {fasta} ===")
        print(f"Length: {len(seq)}")
        if hits:
            print(f"Best hit (by LLR): start={hits[0][0]} window={hits[0][1]} llr={hits[0][2]:.4f}")
            print("Top 10 hits:")
            for r, (pos, win, sc) in enumerate(hits, start=1):
                print(f"  #{r:02d} start={pos:<8} llr={sc:>9.4f}  {win}")
        else:
            print("No valid windows (maybe many Ns).")

        out_png = f"flu{idx}_signal.png"
        plot_signal(f"flu{idx}", P_signal, hits, out_png)
        print(f"Saved chart: {out_png}")

    print("\nDone.")


if __name__ == "__main__":
    main()
