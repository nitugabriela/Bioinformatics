import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import math

def calculate_tm_simple(dna_seq):
    dna_seq = dna_seq.upper()
    a_count = dna_seq.count('A')
    t_count = dna_seq.count('T')
    g_count = dna_seq.count('G')
    c_count = dna_seq.count('C')
    tm = 4 * (g_count + c_count) + 2 * (a_count + t_count)
    return tm

def calculate_tm_advanced(dna_seq, Na = 0.0001):
    dna_seq = dna_seq.upper()
    GC_percent = (dna_seq.count('G') + dna_seq.count('C')) / len(dna_seq) * 100
    tm = - (81.5 + 16.6 * math.log10(Na) + 0.41 * GC_percent - (600 / len(dna_seq)))
    return tm

def sliding_window(dna_seq, window = 8):
    tm_result = []
    for start in range(len(dna_seq) - window + 1):
        window_seq = dna_seq[start:start + window]
        if all(c in 'ATGC' for c in window_seq):
            tm = calculate_tm_simple(window_seq)
            tm_adv = calculate_tm_advanced(window_seq)
            tm_result.append((start + 1, window_seq, tm, tm_adv))
    return tm_result

def read_FASTA(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return seq

def open_and_process():
    path = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files", "*.fa *.fasta"), ("All files", "*.*")]
    )
    if not path:
        return
    try:
        seq = read_FASTA(path)
    except Exception as e:
        messagebox.showerror("Error", f"Failed to read file:\n{e}")
        return
    results = sliding_window(seq, window=8)
    text_area.config(state="normal")
    text_area.delete("1.0", "end")
    text_area.insert("end", f"File: {path}\n")
    text_area.insert("end", f"Sequence length: {len(seq)}\n")
    text_area.insert("end", f"Windows (window=8): {len(results)}\n\n")
    text_area.insert("end", "Start\tWindow\tTm_simple (°C)\tTm_advanced (°C)\n")
    text_area.insert("end", "-" * 80 + "\n")
    for start, wseq, tm, tm_adv in results:
        adv_str = f"{tm_adv:.2f}°C" if tm_adv is not None else "N/A"
        text_area.insert("end", f"{start}\t{wseq}\t{tm}°C\t{adv_str}\n")
    if not results:
        text_area.insert("end", "\nNo valid 8-nt windows with only A/T/G/C were found.\n")
    text_area.config(state="disabled")

def create_gui():
    root = tk.Tk()
    root.title("Sliding-window Tm (window=8)")
    frm = tk.Frame(root, padx=8, pady=8)
    frm.pack(fill="both", expand=True)
    btn = tk.Button(frm, text="Select FASTA file...", command=open_and_process)
    btn.pack(anchor="w")
    global text_area
    text_area = scrolledtext.ScrolledText(frm, width=100, height=24, state="disabled")
    text_area.pack(fill="both", expand=True, pady=(8,0))
    root.mainloop()

if __name__ == "__main__":
    create_gui()
