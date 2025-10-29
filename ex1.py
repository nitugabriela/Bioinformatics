import random
from collections import defaultdict, Counter

SEQ_PATH = "sequence.fasta"
READS_N = 2000
READ_MIN = 100
READ_MAX = 150
K = 61
MIN_KMER_COUNT = 2
SEED = 42

def read_fasta_single(path):
    s = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(">"):
                s.append(line.strip())
    s = "".join(s).upper()
    s = "".join(c for c in s if c in "ACGT")
    return s

def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def canonical_kmer(km):
    rc = revcomp(km)
    return km if km <= rc else rc

def sample_reads(genome, n, min_len, max_len):
    random.seed(SEED)
    L = len(genome)
    reads = []
    for _ in range(n):
        rlen = random.randint(min_len, max_len)
        i = random.randint(0, L - rlen)
        r = genome[i:i + rlen]
        if random.random() < 0.5:
            r = revcomp(r)
        reads.append(r)
    return reads

def build_dbg(reads, k=K, canonical=True, min_kmer_count=MIN_KMER_COUNT):
    kc = Counter()
    for r in reads:
        if len(r) < k:
            continue
        for i in range(len(r) - k + 1):
            km = r[i:i + k]
            if canonical:
                km = canonical_kmer(km)
            kc[km] += 1
    keep = {km for km, c in kc.items() if c >= min_kmer_count}
    edges = defaultdict(list)
    indeg = Counter()
    outdeg = Counter()
    for km in keep:
        u, v = km[:-1], km[1:]
        edges[u].append(v)
        outdeg[u] += 1
        indeg[v] += 1
        _ = indeg[u]
        _ = outdeg[v]
    return edges, indeg, outdeg

def euler_paths(edges, indeg, outdeg):
    E = {u: list(vs) for u, vs in edges.items()}
    nodes = set(indeg) | set(outdeg)
    starts = [n for n in nodes if outdeg[n] - indeg[n] == 1]
    if not starts:
        starts = [n for n, d in outdeg.items() if d > 0]
    paths = []
    for s in starts:
        if not E.get(s):
            continue
        st = [s]
        circuit = []
        while st:
            v = st[-1]
            if E.get(v):
                st.append(E[v].pop())
            else:
                circuit.append(st.pop())
        circuit.reverse()
        if len(circuit) > 1:
            paths.append(circuit)
    return paths

def path_to_seq(nodes):
    if not nodes:
        return ""
    seq = nodes[0]
    for n in nodes[1:]:
        seq += n[-1]
    return seq

def assemble_contigs(reads, k=K, canonical=True, min_kmer_count=MIN_KMER_COUNT):
    edges, indeg, outdeg = build_dbg(reads, k, canonical, min_kmer_count)
    paths = euler_paths(edges, indeg, outdeg)
    contigs = [path_to_seq(p) for p in paths]
    contigs.sort(key=len, reverse=True)
    return contigs

def write_fasta(seqs, path, prefix):
    with open(path, "w") as f:
        for i, s in enumerate(seqs, 1):
            f.write(f">{prefix}_{i}\n")
            for j in range(0, len(s), 80):
                f.write(s[j:j + 80] + "\n")

def main():
    genome = read_fasta_single(SEQ_PATH)
    print("Sequence length:", len(genome))
    reads = sample_reads(genome, READS_N, READ_MIN, READ_MAX)
    write_fasta(reads, "reads.fasta", "read")
    contigs = assemble_contigs(reads, K, True, MIN_KMER_COUNT)
    if not contigs:
        print("No contigs assembled. Try a smaller K or lower MIN_KMER_COUNT.")
    else:
        write_fasta(contigs, "reconstructed_contigs.fasta", "contig")
        print("Contigs:", len(contigs), "| Longest:", len(contigs[0]))

if __name__ == "__main__":
    main()
