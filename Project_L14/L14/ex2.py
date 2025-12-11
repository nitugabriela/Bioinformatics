import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import numpy as np

class GenomeUtils:
    @staticmethod
    def read_fasta(filename):
        sequence = ""
        try:
            with open(filename, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line.startswith(">"):
                        sequence += line
            return sequence.upper()
        except FileNotFoundError:
            return None

class AlignmentLayer:
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.k = 6

    def find_potential_matches(self):
        n = len(self.seq1)
        m = len(self.seq2)
        match_coords = []
        
        kmer_index = {}
        for i in range(0, m - self.k + 1):
            kmer = self.seq2[i:i+self.k]
            if kmer not in kmer_index:
                kmer_index[kmer] = []
            kmer_index[kmer].append(i)
            
        step = 10
        for j in range(0, n - self.k + 1, step):
            kmer = self.seq1[j:j+self.k]
            if kmer in kmer_index:
                for pos_in_seq2 in kmer_index[kmer]:
                    match_coords.append((j, pos_in_seq2))
                    
        return match_coords

class GenomeApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Genome Aligner - Density Simulation")
        self.root.geometry("1100x700")
        
        self.control_frame = ttk.Frame(root, padding="10")
        self.control_frame.pack(side=tk.LEFT, fill=tk.Y)
        
        self.viz_frame = ttk.Frame(root, padding="10")
        self.viz_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        
        ttk.Label(self.control_frame, text="Genome 1 (X-axis):").pack(anchor=tk.W)
        self.lbl_g1 = ttk.Label(self.control_frame, text="Not Loaded", foreground="red")
        self.lbl_g1.pack(anchor=tk.W, pady=5)

        ttk.Label(self.control_frame, text="Genome 2 (Y-axis):").pack(anchor=tk.W)
        self.lbl_g2 = ttk.Label(self.control_frame, text="Not Loaded", foreground="red")
        self.lbl_g2.pack(anchor=tk.W, pady=5)
        
        self.btn_load = ttk.Button(self.control_frame, text="1. Load Genomes", command=self.load_data)
        self.btn_load.pack(fill=tk.X, pady=10)
        
        self.btn_sim = ttk.Button(self.control_frame, text="2. Run Simulation (Heatmap)", command=self.run_simulation, state=tk.DISABLED)
        self.btn_sim.pack(fill=tk.X, pady=10)
        
        self.log_text = tk.Text(self.control_frame, height=20, width=40, font=("Courier", 8))
        self.log_text.pack(pady=10)

        self.fig, self.ax = plt.subplots(figsize=(6, 6))
        self.ax.set_title("Waiting for Data...")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.viz_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        self.seq1 = ""
        self.seq2 = ""

    def log(self, msg):
        self.log_text.insert(tk.END, msg + "\n")
        self.log_text.see(tk.END)
        self.root.update()

    def load_data(self):
        self.log("Checking for files in script directory...")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        
        file1_name = "covid1.fasta"
        file2_name = "flu1.fasta"
        
        path1 = os.path.join(script_dir, file1_name)
        path2 = os.path.join(script_dir, file2_name)
        
        self.seq1 = GenomeUtils.read_fasta(path1)
        self.seq2 = GenomeUtils.read_fasta(path2)
        
        if self.seq1 and self.seq2:
            self.lbl_g1.config(text=f"{file1_name} ({len(self.seq1)} bp)", foreground="green")
            self.lbl_g2.config(text=f"{file2_name} ({len(self.seq2)} bp)", foreground="green")
            self.btn_sim.config(state=tk.NORMAL)
            self.log(f"Loaded successfully.")
        else:
            self.log(f"Error: Could not find files.")

    def run_simulation(self):
        if not self.seq1 or not self.seq2:
            return
            
        self.log("--- Starting Simulation ---")
        layer = AlignmentLayer(self.seq1, self.seq2)
        
        self.log("Step 1: Sampling matches (Seeds)...")
        matches = layer.find_potential_matches()
        
        self.log(f"Step 2: Found {len(matches)} seeds.")
        self.log("Step 3: Generating Density Heatmap...")
        
        x_val = [m[0] for m in matches]
        y_val = [m[1] for m in matches]
        
        self.ax.clear()
        
        h = self.ax.hist2d(x_val, y_val, bins=[50, 50], cmap='plasma', cmin=1)
        
        if hasattr(self, 'cbar') and self.cbar:
            self.cbar.remove()
        self.cbar = self.fig.colorbar(h[3], ax=self.ax)
        self.cbar.set_label('Match Density')
        
        self.ax.set_title("Genome Homology Simulation\n(Heatmap of Matching Regions)")
        self.ax.set_xlabel(f"Genome 1 Position")
        self.ax.set_ylabel(f"Genome 2 Position")
        self.ax.set_xlim(0, len(self.seq1))
        self.ax.set_ylim(0, len(self.seq2))
        self.ax.grid(True, linestyle='--', alpha=0.3)
        
        self.log("Simulation Complete.")
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = GenomeApp(root)
    root.mainloop()