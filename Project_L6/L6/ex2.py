import os, math, re
import matplotlib.pyplot as plt

SEQ_PATH = os.path.join(os.path.dirname(__file__), "sequence.fasta")
OUT_CSV  = "digest_fragments.csv"
OUT_IMG  = "gel_digest.png"

ENZYMES = [
    ("EcoRI",  "GAATTC",   1),
    ("BamHI",  "GGATCC",   1),
    ("HindIII","AAGCTT",   1),
    ("HaeIII", "GGCC",     2),
    ("NotI",   "GCGGCCGC", 2),
]

def read_fasta_single(path):
    s=[]
    with open(path,"r",encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                s.append(line.strip())
    s="".join(s).upper()
    return "".join(c for c in s if c in "ACGT")

def revcomp(s):
    return s.translate(str.maketrans("ACGT","TGCA"))[::-1]

def find_cuts(seq, site, cut_offset):
    pats = {site, revcomp(site)}
    cuts=set()
    for pat in pats:
        for m in re.finditer(re.escape(pat), seq):
            cuts.add(m.start()+cut_offset)
    cuts=[c for c in cuts if 0 < c < len(seq)]
    return sorted(cuts)

def digest(seq, site, cut_offset):
    cuts = find_cuts(seq, site, cut_offset)
    if not cuts:
        return [(0, len(seq), seq)]
    starts = [0] + cuts
    ends   = cuts + [len(seq)]
    return [(s, e-s, seq[s:e]) for s, e in zip(starts, ends)]

def size_to_position(bp, bp_min=100, bp_max=3000):
    bp = max(bp_min, min(bp_max, bp))
    lo, hi = math.log10(bp_min), math.log10(bp_max)
    y = (math.log10(bp)-lo)/(hi-lo)
    return 1.0 - y

def plot_gel(all_frag_sizes, img_path):
    fig = plt.figure(figsize=(8.5, 9), dpi=120)
    ax = plt.gca()
    ax.set_xlim(0, 1.55)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ladder_x, ladder_w = 0.05, 0.18
    lane_w = 0.18
    lane_gap = 0.12
    lane_xs=[]
    x = 0.32
    for _ in all_frag_sizes:
        lane_xs.append(x)
        x += lane_w + lane_gap
    ax.add_patch(plt.Rectangle((ladder_x, 0.05), ladder_w, 0.90,
                               facecolor="#1f1f2e", edgecolor="#111111", linewidth=2))
    ax.add_patch(plt.Rectangle((ladder_x+0.01, 0.94), ladder_w-0.02, 0.02,
                               facecolor="none", edgecolor="#e6e6e6", linewidth=1.2))
    ladder_bps = [3000, 2000, 1500, 1200, 1000, 800, 700, 600, 500, 400, 300, 200]
    for bp in ladder_bps:
        y = 0.05 + 0.90*size_to_position(bp)
        th = 0.007
        ax.add_patch(plt.Rectangle((ladder_x+0.01, y-th/2), ladder_w-0.02, th,
                                   facecolor="#e8e8f0", edgecolor="#cfcfe8"))
        ax.text(ladder_x-0.03, y, f"{bp} bp", va="center", ha="right",
                fontsize=10, color="#000000", fontweight="bold")
    ax.text(ladder_x+ladder_w/2, 0.965, "Ladder", ha="center", va="bottom",
            fontsize=10, color="#000000", fontweight="bold")
    for (label, sizes), lx in zip(all_frag_sizes, lane_xs):
        ax.add_patch(plt.Rectangle((lx, 0.05), lane_w, 0.90,
                                   facecolor="#1f1f2e", edgecolor="#111111", linewidth=2))
        ax.add_patch(plt.Rectangle((lx+0.02, 0.94), lane_w-0.04, 0.02,
                                   facecolor="none", edgecolor="#e6e6e6", linewidth=1.5))
        for s in sorted(sizes, reverse=True):
            y = 0.05 + 0.90*size_to_position(s)
            th = max(0.004, 0.016 - 0.000003*s)
            ax.add_patch(plt.Rectangle((lx+0.01, y-th/2), lane_w-0.02, th,
                                       facecolor="#e8e8f0", edgecolor="#cfcfe8"))
            ax.text(lx + lane_w + 0.02, y, f"{s} bp",va="center", ha="left", fontsize=8, color="#000000")
        ax.text(lx+lane_w/2, 0.965, label, ha="center", va="bottom",
                fontsize=10, color="#000000", fontweight="bold")
    fig.tight_layout()
    fig.savefig(img_path)
    plt.close(fig)

def main():
    if not os.path.exists(SEQ_PATH):
        print(f"Missing {SEQ_PATH}.")
        return
    seq = read_fasta_single(SEQ_PATH)
    if not seq:
        print("Empty or invalid FASTA.")
        return
    rows=[]
    lanes=[]
    for name, site, cut in ENZYMES:
        frags = digest(seq, site, cut)
        sizes = [L for _, L, _ in frags]
        lanes.append((name, sizes))
        for idx,(start,L,_) in enumerate(frags,1):
            rows.append((name, idx, start, L))
    rows.sort(key=lambda r:(r[0], r[1]))
    with open(OUT_CSV,"w",encoding="utf-8") as f:
        f.write("enzyme,fragment_index,start,length_bp\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]}\n")
    plot_gel(lanes, OUT_IMG)
    print(f"Wrote {OUT_CSV} and {OUT_IMG}.")
    print("Enzymes used:", ", ".join(e[0] for e in ENZYMES))

if __name__ == "__main__":
    main()
