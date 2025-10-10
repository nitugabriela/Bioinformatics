import tkinter as tk
from tkinter import filedialog, messagebox
from collections import Counter
import matplotlib.pyplot as plt
import os

ALPHABET = {"A", "C", "G", "T"}

# -------- FASTA loading --------
def load_fasta():
    path = filedialog.askopenfilename(
        title="Choose FASTA file",
        filetypes=[("FASTA files", "*.fa *.fasta *.fna *.ffn *.frn"), ("All files", "*.*")]
    )
    if not path:
        return "", ""
    try:
        parts = []
        with open(path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    continue
                parts.append(line.strip())
        seq = "".join(parts).upper()
        if not seq:
            raise ValueError("No sequence content found.")
        return seq, path
    except Exception as e:
        messagebox.showerror("Read error", f"Could not read FASTA:\n{e}")
        return "", ""

# -------- Sliding-window frequencies (A/C/G/T only) --------
def sliding_window_freqs(sequence, window_size=30):
    n = len(sequence)
    if window_size <= 0:
        raise ValueError("Window size must be a positive integer.")
    if n < window_size:
        raise ValueError(f"Sequence too short ({n}) for window size {window_size}.")

    freqs = {b: [] for b in ALPHABET}
    centers = []

    for i in range(n - window_size + 1):
        chunk = sequence[i:i+window_size]
        counts = Counter(ch for ch in chunk if ch in ALPHABET)
        denom = sum(counts[b] for b in ALPHABET)

        if denom == 0:
            for b in ALPHABET:
                freqs[b].append(0.0)
        else:
            for b in ALPHABET:
                freqs[b].append(counts[b] / denom)

        # center position of window (1-based for readability)
        centers.append(i + window_size // 2 + 1)

    return centers, freqs

# -------- Simple moving-average smoothing --------
def smooth(values, window=5):
    """Centered moving average (handles edges by shrinking the window)."""
    if window is None or window < 2:
        return values
    m = []
    half = window // 2
    n = len(values)
    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        m.append(sum(values[start:end]) / (end - start))
    return m

# -------- Plotting --------
def plot_freqs(x, freqs, title, smooth_window=5):
    plt.figure(figsize=(12, 6))
    for base in sorted(freqs.keys()):  # A, C, G, T
        y = smooth(freqs[base], smooth_window) if smooth_window and smooth_window > 1 else freqs[base]
        plt.plot(x, y, label=base, linewidth=1.8, alpha=0.95)

    plt.title(title)
    plt.xlabel("Position (bp, window center)")
    plt.ylabel("Relative Frequency")
    plt.ylim(0, 1)
    plt.grid(True, linestyle="--", linewidth=0.6, alpha=0.7)
    plt.legend(title="Base")
    plt.tight_layout()
    plt.show()

# -------- GUI wiring --------
def run_analysis(window_entry, smooth_entry):
    # Parse window size
    try:
        w = int(window_entry.get())
        if w <= 0:
            raise ValueError
    except ValueError:
        messagebox.showerror("Invalid input", "Window size must be a positive integer.")
        return

    # Parse smoothing window (allow 1 = no smoothing)
    try:
        sw = int(smooth_entry.get())
        if sw < 1:
            raise ValueError
    except ValueError:
        messagebox.showerror("Invalid input", "Smoothing window must be an integer ≥ 1.")
        return

    sequence, path = load_fasta()
    if not sequence:
        return

    try:
        x, freqs = sliding_window_freqs(sequence, w)
    except Exception as e:
        messagebox.showerror("Error", str(e))
        return

    title = f"{os.path.basename(path)} — window={w}, smooth={sw}"
    plot_freqs(x, freqs, title=title, smooth_window=sw)

# -------- Main window --------
root = tk.Tk()
root.title("DNA Sliding-Window Frequency Analyzer")

frm = tk.Frame(root)
frm.pack(padx=12, pady=12)

tk.Label(frm, text="Window size:").grid(row=0, column=0, sticky="e", padx=(0,6))
window_entry = tk.Entry(frm, width=6)
window_entry.insert(0, "30")
window_entry.grid(row=0, column=1, padx=(0,12))

tk.Label(frm, text="Smoothing window:").grid(row=0, column=2, sticky="e", padx=(0,6))
smooth_entry = tk.Entry(frm, width=6)
smooth_entry.insert(0, "5")  # set to "1" for no smoothing
smooth_entry.grid(row=0, column=3, padx=(0,12))

btn = tk.Button(frm, text="Choose FASTA & Analyze", command=lambda: run_analysis(window_entry, smooth_entry))
btn.grid(row=0, column=4)

root.mainloop()
