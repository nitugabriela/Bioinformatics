import matplotlib.pyplot as plt
from collections import Counter
from typing import List, Tuple

# 1-2

def sliding_windows(seq: str, win: int) -> List[str]:
    seq = seq.upper()
    return [seq[i:i+win] for i in range(len(seq) - win + 1)]

# 3

def compute_cgtot(seq: str) -> float:
    seq = seq.upper()
    counts = Counter(seq)
    total = counts['A'] + counts['T'] + counts['C'] + counts['G']
    cg = counts['C'] + counts['G']
    if total == 0:
        return 0.0
    cgtot = 100.0 * cg / total
    return round(cgtot, 2)

def compute_cgsw_values(seq: str, win: int) -> List[float]:
    seq = seq.upper()
    counts = Counter(seq)
    total_tot = counts['A'] + counts['T'] + counts['C'] + counts['G']
    cg_tot = counts['C'] + counts['G']
    if total_tot == 0 or cg_tot == 0:
        return []
    cgtot = 100.0 * cg_tot / total_tot 

    windows = sliding_windows(seq, win)
    cgsw_values: List[float] = []
    for w in windows:
        c = w.count('C')
        g = w.count('G')
        total_w = len(w)
        cg_sw = c + g
        if total_w == 0 or cg_sw == 0:
            cgsw_values.append(0.0)
        else:
            cgsw = cgtot * cg_sw / total_w
            cgsw_values.append(cgsw)
    return cgsw_values

# 4
def kappa_ic_window(window_seq: str) -> float:
    A = window_seq.upper()
    N = len(A) - 1
    if N <= 0:
        return 0.0

    T = 0.0
    for u in range(1, N + 1):  
        B = A[u:]              
        C = 0
        for i in range(len(B)): 
            if A[i] == B[i]:
                C += 1
        T += (C / len(B)) * 100.0

    IC = round(T / N, 2)
    return IC

def compute_kappa_values(seq: str, win: int) -> List[float]:
    windows = sliding_windows(seq, win)
    return [kappa_ic_window(w) for w in windows]

# 5-7

def promoter_pattern(seq: str, win: int) -> Tuple[List[float], List[float]]:
    x_vals = compute_cgsw_values(seq, win)
    y_vals = compute_kappa_values(seq, win)
    n = min(len(x_vals), len(y_vals))
    return x_vals[:n], y_vals[:n]

def center_of_weight(xs: List[float], ys: List[float]) -> Tuple[float, float]:
    if not xs or not ys:
        return (0.0, 0.0)
    n = min(len(xs), len(ys))
    xs = xs[:n]
    ys = ys[:n]
    cx = sum(xs) / n
    cy = sum(ys) / n
    return (cx, cy)

def plot_pattern(xs: List[float], ys: List[float],
                 center: Tuple[float, float] | None = None,
                 title: str = "Promoter pattern") -> None:
    plt.figure()
    plt.scatter(xs, ys, s=10)
    if center is not None:
        cx, cy = center
        plt.scatter([cx], [cy], marker="x", s=80)
    plt.xlabel("(C+G)% (CGSW)")
    plt.ylabel("Kappa IC")
    plt.title(title)
    plt.grid(True)

def plot_centers(centers: List[Tuple[float, float]],
                 title: str = "Centers of promoter patterns") -> None:
    if not centers:
        return
    xs = [c[0] for c in centers]
    ys = [c[1] for c in centers]
    plt.figure()
    plt.scatter(xs, ys)
    plt.xlabel("Center (C+G)%")
    plt.ylabel("Center Kappa IC")
    plt.title(title)
    plt.grid(True)

# 8

if __name__ == "__main__":
    S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"

    window_size = 30

    CG = compute_cgtot(S)
    print(f"CGTOT (C+G)% for S = {CG}") 

    IC = kappa_ic_window(S)
    print(f"Kappa IC for S = {IC}")  

    xs, ys = promoter_pattern(S, window_size)

    center_S = center_of_weight(xs, ys)
    print(f"Center of weight for S pattern: {center_S}")

    plot_pattern(xs, ys, center=center_S,
                 title="Pattern for test sequence S")
    
    centers = [center_S]
    plot_centers(centers, title="Centers of patterns (example with S only)")

    plt.show()
