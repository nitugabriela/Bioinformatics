import random
from collections import defaultdict, Counter

SEQ_PATH = "sequence.fasta"
READS_N = 2000
READ_MIN = 100
READ_MAX = 150
KS_TO_TRY = (41, 51, 61, 71, 81, 91)
MIN_KMER_COUNT = 2
SEED = 42

def read_fasta_single(path):
    s = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                s.append(line.strip())
    s = "".join(s).upper()
    return "".join(c for c in s if c in "ACGT")

def write_fasta(seqs, path, prefix):
    with open(path, "w", encoding="utf-8") as f:
        for i, s in enumerate(seqs, 1):
            f.write(f">{prefix}_{i}\n")
            for j in range(0, len(s), 80):
                f.write(s[j:j + 80] + "\n")

def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def canonical_kmer(km):
    rc = revcomp(km)
    return km if km <= rc else rc

def path_to_seq(nodes):
    if not nodes:
        return ""
    seq = nodes[0]
    for n in nodes[1:]:
        seq += n[-1]
    return seq

def sample_reads(genome, n, min_len, max_len, seed=SEED):
    random.seed(seed)
    L = len(genome)
    reads = []
    for _ in range(n):
        rlen = random.randint(min_len, max_len)
        i = random.randint(0, L - rlen)
        r = genome[i:i + rlen]
        reads.append(r)
    return reads

def build_dbg(reads, k, canonical=False, min_kmer_count=2):
    kc = Counter()
    for r in reads:
        if len(r) < k: 
            continue
        for i in range(len(r) - k + 1):
            km = r[i:i+k]
            kc[km] += 1
    edges = defaultdict(list)
    indeg = Counter()
    outdeg = Counter()
    for km, c in kc.items():
        if c < min_kmer_count:
            continue
        u, v = km[:-1], km[1:]
        edges[u].append(v)
        outdeg[u] += 1
        indeg[v] += 1
        _ = indeg[u]
        _ = outdeg[v]
    return edges, indeg, outdeg

def _pick_start(nodes, E, indeg, outdeg):
    for n in nodes:
        if outdeg[n] - indeg[n] == 1 and E.get(n):
            return n
    for n in nodes:
        if E.get(n):
            return n
    return None

def euler_all_contigs(edges):
    E = {u: list(vs) for u, vs in edges.items()}
    indeg = Counter()
    outdeg = Counter()
    for u, vs in E.items():
        outdeg[u] += len(vs)
        for v in vs:
            indeg[v] += 1
            _ = indeg[u]
            _ = outdeg[v]
    nodes = set(indeg) | set(outdeg)
    contigs = []
    while True:
        s = _pick_start(nodes, E, indeg, outdeg)
        if s is None:
            break
        st = [s]
        circuit = []
        while st:
            v = st[-1]
            if E.get(v) and E[v]:
                st.append(E[v].pop())
            else:
                circuit.append(st.pop())
                if not E.get(v):
                    E.pop(v, None)
        circuit.reverse()
        if len(circuit) > 1:
            contigs.append(path_to_seq(circuit))
    contigs.sort(key=len, reverse=True)
    return contigs

def assemble_best(reads, ks=KS_TO_TRY, min_kmer_count=MIN_KMER_COUNT):
    best_contigs = []
    best_k = None
    for k in ks:
        if k < 2 or k >= READ_MIN:
            continue
        edges, indeg, outdeg = build_dbg(reads, k=k, canonical=False, min_kmer_count=min_kmer_count)
        if not edges:
            continue
        contigs = euler_all_contigs(edges)
        if contigs and (not best_contigs or len(contigs[0]) > len(best_contigs[0])):
            best_contigs = contigs
            best_k = k
    return best_k, best_contigs

def main():
    genome = read_fasta_single(SEQ_PATH)
    if not genome:
        print(f"Could not read a DNA sequence from {SEQ_PATH}.")
        return
    print("Input sequence length:", len(genome))
    reads = sample_reads(genome, READS_N, READ_MIN, READ_MAX, seed=SEED)
    write_fasta(reads, "reads.fasta", "read")
    print(f"Sampled {len(reads)} reads (len {READ_MIN}-{READ_MAX}).")
    best_k, contigs = assemble_best(reads, ks=KS_TO_TRY, min_kmer_count=MIN_KMER_COUNT)
    if not contigs:
        print("No contigs assembled. Try a smaller K, lower MIN_KMER_COUNT, or longer reads.")
        return
    write_fasta(contigs, "reconstructed_contigs.fasta", "contig")
    print(f"Best K: {best_k}")
    print("Assembled contigs:", len(contigs))
    print("Longest contig length:", len(contigs[0]))
    genome_rc = revcomp(genome)
    contig = contigs[0]
    print("Forward rotation:", contig in (genome + genome))
    print("Reverse rotation:", contig in (genome_rc + genome_rc))
    print("Δ length:", len(genome) - len(contig))
    approx_input_len = len(genome)
    longest = len(contigs[0])
    if abs(longest - approx_input_len) <= 20:
        print("Longest contig matches the input length (±20 bp). Likely a perfect rebuild (or rotation).")
    else:
        print("Longest contig differs from input length. Likely repeats caused fragmentation or branching.")

if __name__ == "__main__":
    main()

