import os, math, random
from collections import namedtuple
import matplotlib.pyplot as plt

N_FRAGS = 10
FRAG_MIN = 100
FRAG_MAX = 3000
SEED = 7
OUT_IMG = "gel_multi.png"
OUT_CSV = "fragments_multi.csv"
EXCLUDE = {"fragments.fasta", "sequence.fasta"} 

Fragment = namedtuple("Fragment", ["name", "start", "length", "seq"])

def read_fasta_single(path):
    s=[]
    with open(path,"r",encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                s.append(line.strip())
    s="".join(s).upper()
    return "".join(c for c in s if c in "ACGT")

def sample_fragments(genome, n, min_len, max_len, seed=SEED):
    random.seed(seed)
    L=len(genome)
    max_len=min(max_len, L)
    frags=[]
    for i in range(1, n+1):
        rlen=random.randint(min_len, max_len)
        start=random.randint(0, L-rlen)
        seq=genome[start:start+rlen]
        frags.append(Fragment(f"frag_{i}", start, rlen, seq))
    return frags

def size_to_position(bp, bp_min=100, bp_max=3000):
    bp=max(bp_min, min(bp_max, bp))
    lo,hi=math.log10(bp_min), math.log10(bp_max)
    y=(math.log10(bp)-lo)/(hi-hi if hi==hi else 1) 
    return 1.0 - ((math.log10(bp)-lo)/(hi-lo))

def plot_gel(lanes, img_path):
   
    fig = plt.figure(figsize=(2.5 + 1.8*len(lanes), 9), dpi=120)
    ax = plt.gca()
    ax.set_xlim(0, 0.28 + 0.30*len(lanes))
    ax.set_ylim(0, 1)
    ax.axis("off")

    ladder_x, ladder_w = 0.05, 0.18
    ax.add_patch(plt.Rectangle((ladder_x, 0.05), ladder_w, 0.90,
                               facecolor="#1f1f2e", edgecolor="#111111", linewidth=2))
    ax.add_patch(plt.Rectangle((ladder_x+0.01, 0.94), ladder_w-0.02, 0.02,
                               facecolor="none", edgecolor="#e6e6e6", linewidth=1.2))
    ladder_bps=[3000,2000,1500,1200,1000,800,700,600,500,400,300,200]
    for bp in ladder_bps:
        y=0.05+0.90*(1.0-(math.log10(max(100,min(3000,bp)))-math.log10(100))/(math.log10(3000)-math.log10(100)))
        th=0.007
        ax.add_patch(plt.Rectangle((ladder_x+0.01, y-th/2), ladder_w-0.02, th,
                                   facecolor="#e8e8f0", edgecolor="#cfcfe8"))
        ax.text(ladder_x-0.03, y, f"{bp} bp", va="center", ha="right", fontsize=10, color="#000000", fontweight="bold")
    ax.text(ladder_x+ladder_w/2, 0.965, "Ladder", ha="center", va="bottom", fontsize=10, color="#000000", fontweight="bold")

    x = 0.28
    lane_w = 0.18
    gap = 0.12
    for label, sizes in lanes:
        ax.add_patch(plt.Rectangle((x, 0.05), lane_w, 0.90,
                                   facecolor="#1f1f2e", edgecolor="#111111", linewidth=2))
        ax.add_patch(plt.Rectangle((x+0.02, 0.94), lane_w-0.04, 0.02,
                                   facecolor="none", edgecolor="#e6e6e6", linewidth=1.5))
        for s in sorted(sizes, reverse=True):
            y = 0.05 + 0.90*(1.0-(math.log10(max(100,min(3000,s)))-math.log10(100))/(math.log10(3000)-math.log10(100)))
            th = max(0.004, 0.016 - 0.000003*s)
            ax.add_patch(plt.Rectangle((x+0.01, y-th/2), lane_w-0.02, th,
                                       facecolor="#e8e8f0", edgecolor="#cfcfe8"))
        ax.text(x+lane_w/2, 0.965, label, ha="center", va="bottom", fontsize=10, color="#000000", fontweight="bold")
        x += lane_w + gap

    fig.tight_layout()
    fig.savefig(img_path)
    plt.close(fig)

def main():
    fasta_files=[f for f in os.listdir(".") if f.endswith(".fasta") and f not in EXCLUDE]
    if not fasta_files:
        print("No FASTA files found.")
        return

    rows=[]
    lanes=[]
    for fname in sorted(fasta_files):
        seq=read_fasta_single(fname)
        if not seq:
            continue
        frags=sample_fragments(seq, N_FRAGS, FRAG_MIN, FRAG_MAX, seed=SEED)
        sizes=[fr.length for fr in frags]
        lanes.append((os.path.splitext(fname)[0], sizes))
        for fr in frags:
            rows.append((fname, fr.name, fr.start, fr.length))

    with open(OUT_CSV,"w",encoding="utf-8") as f:
        f.write("sequence,fragment_name,start,length_bp\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]},{r[2]},{r[3]}\n")

    plot_gel(lanes, OUT_IMG)
    print(f"Processed {len(lanes)} sequences.")
    print(f"Wrote {OUT_IMG} and {OUT_CSV}.")

if __name__ == "__main__":
    main()
