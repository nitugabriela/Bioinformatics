import random, math, os
from collections import namedtuple
import matplotlib.pyplot as plt

SEQ_PATH = os.path.join(os.path.dirname(__file__), "sequence.fasta")
N_FRAGS = 10
FRAG_MIN = 100
FRAG_MAX = 3000
SEED = 7
OUT_FASTA = "fragments.fasta"
OUT_IMG = "gel.png"

Fragment = namedtuple("Fragment", ["name", "start", "length", "seq"])

def read_fasta_single(path):
    s = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                s.append(line.strip())
    s = "".join(s).upper()
    return "".join(c for c in s if c in "ACGT")

def write_fasta(frags, path):
    with open(path, "w", encoding="utf-8") as f:
        for fr in frags:
            f.write(f">{fr.name}|start={fr.start}|len={fr.length}\n")
            for i in range(0, len(fr.seq), 80):
                f.write(fr.seq[i:i+80] + "\n")

def sample_fragments(genome, n, min_len, max_len, seed=SEED):
    random.seed(seed)
    L = len(genome)
    max_len = min(max_len, L)
    frags = []
    for i in range(1, n+1):
        rlen = random.randint(min_len, max_len)
        start = random.randint(0, L - rlen)
        seq = genome[start:start+rlen]
        frags.append(Fragment(f"frag_{i}", start, rlen, seq))
    return frags

def size_to_position_bp(bp, bp_min=100, bp_max=3000):
    bp = max(bp_min, min(bp_max, bp))
    lo = math.log10(bp_min)
    hi = math.log10(bp_max)
    y_norm = (math.log10(bp) - lo) / (hi - lo)
    return 1.0 - y_norm

def plot_gel(frags, img_path):
    import math
    import matplotlib.pyplot as plt

    def size_to_position_bp(bp, bp_min=100, bp_max=3000):
        bp = max(bp_min, min(bp_max, bp))
        lo, hi = math.log10(bp_min), math.log10(bp_max)
        y = (math.log10(bp) - lo) / (hi - lo)
        return 1.0 - y

    fig = plt.figure(figsize=(6.5, 9), dpi=120)
    ax = plt.gca()
    ax.set_xlim(0, 1.45)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ladder_x, ladder_w = 0.10, 0.18
    sample_x, sample_w = 0.40, 0.50

    ax.add_patch(plt.Rectangle((ladder_x, 0.05), ladder_w, 0.90,
                               facecolor="#1f1f2e", edgecolor="#111111", linewidth=2))
    ax.add_patch(plt.Rectangle((ladder_x+0.01, 0.94), ladder_w-0.02, 0.02,
                               facecolor="none", edgecolor="#e6e6e6", linewidth=1.2))

    ladder_bps = [3000, 2000, 1500, 1200, 1000, 800, 700, 600, 500, 400, 300, 200]
    for bp in ladder_bps:
        y = 0.05 + 0.90 * size_to_position_bp(bp)
        th = 0.007
        ax.add_patch(plt.Rectangle((ladder_x+0.01, y - th/2),
                                   ladder_w-0.02, th,
                                   facecolor="#e8e8f0", edgecolor="#cfcfe8"))
        ax.text(ladder_x - 0.03, y, f"{bp} bp", va="center", ha="right",
                fontsize=10, color="#000000", fontweight="bold")

    ax.add_patch(plt.Rectangle((sample_x, 0.05), sample_w, 0.90,
                               facecolor="#1f1f2e", edgecolor="#111111", linewidth=2))
    ax.add_patch(plt.Rectangle((sample_x+0.02, 0.94), sample_w-0.04, 0.02,
                               facecolor="none", edgecolor="#e6e6e6", linewidth=1.5))

    for fr in sorted(frags, key=lambda f: f.length, reverse=True):
        y = 0.05 + 0.90 * size_to_position_bp(fr.length)
        thickness = max(0.004, 0.016 - 0.000003 * fr.length)
        ax.add_patch(plt.Rectangle((sample_x+0.01, y - thickness/2),
                                   sample_w-0.02, thickness,
                                   facecolor="#e8e8f0", alpha=0.95, edgecolor="#cfcfe8"))
        ax.text(sample_x + sample_w + 0.03, y, f"{fr.name} — {fr.length} bp",
                va="center", ha="left", fontsize=9, color="#000000", fontweight="bold")

    fig.tight_layout()
    fig.savefig(img_path)
    plt.close(fig)

def main():
    if not os.path.exists(SEQ_PATH):
        print(f"Missing {SEQ_PATH}. Put a 1–3 kb sequence there first.")
        return
    genome = read_fasta_single(SEQ_PATH)
    if not genome:
        print("Could not read a DNA sequence from sequence.fasta.")
        return
    print("Input sequence length:", len(genome))
    frags = sample_fragments(genome, N_FRAGS, FRAG_MIN, FRAG_MAX, seed=SEED)
    write_fasta(frags, OUT_FASTA)
    plot_gel(frags, OUT_IMG)
    print(f"Generated {len(frags)} fragments (100–3000 bp, clipped to genome length).")
    print(f"Wrote: {OUT_FASTA}")
    print(f"Wrote: {OUT_IMG}")

if __name__ == "__main__":
    main()
