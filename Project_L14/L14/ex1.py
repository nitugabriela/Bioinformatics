import tkinter as tk
from tkinter import ttk, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.colors as mcolors

def run_needleman_wunsch_core(s1, s2, gap_pen, match_sc, mismatch_sc):
    n = len(s1)
    m = len(s2)

    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    tb_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(m + 1):
        score_matrix[i][0] = i * gap_pen
        tb_matrix[i][0] = 2 
    for j in range(n + 1):
        score_matrix[0][j] = j * gap_pen
        tb_matrix[0][j] = 3 
    tb_matrix[0][0] = 0 

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s2[i - 1] == s1[j - 1]:
                diagonal_score = score_matrix[i - 1][j - 1] + match_sc
            else:
                diagonal_score = score_matrix[i - 1][j - 1] + mismatch_sc
            
            up_score = score_matrix[i - 1][j] + gap_pen
            left_score = score_matrix[i][j - 1] + gap_pen
            
            max_score = max(diagonal_score, up_score, left_score)
            score_matrix[i][j] = max_score

            if max_score == diagonal_score:
                tb_matrix[i][j] = 1
            elif max_score == up_score:
                tb_matrix[i][j] = 2
            else:
                tb_matrix[i][j] = 3

    align_s1 = ""
    align_s2 = ""
    i, j = m, n
    match_count = 0
    path_coords = set() 
    path_coords.add((i, j))

    while i > 0 and j > 0:
        current_score = score_matrix[i][j]
        diagonal_val = score_matrix[i - 1][j - 1]
        up_val = score_matrix[i - 1][j]
        left_val = score_matrix[i][j - 1]
        
        is_match = (s2[i - 1] == s1[j - 1])
        step_score = match_sc if is_match else mismatch_sc
        
        if current_score == diagonal_val + step_score:
            align_s1 += s1[j - 1]
            align_s2 += s2[i - 1]
            if is_match: match_count += 1
            i -= 1; j -= 1
        elif current_score == up_val + gap_pen:
            align_s1 += "-"
            align_s2 += s2[i - 1]
            i -= 1
        elif current_score == left_val + gap_pen:
            align_s1 += s1[j - 1]
            align_s2 += "-"
            j -= 1
        path_coords.add((i, j))

    while i > 0:
        align_s1 += "-"
        align_s2 += s2[i - 1]
        i -= 1
        path_coords.add((i, j))
    while j > 0:
        align_s1 += s1[j - 1]
        align_s2 += "-"
        j -= 1
        path_coords.add((i, j))

    align_s1 = align_s1[::-1]
    align_s2 = align_s2[::-1]
    
    results = {
        's1_aligned': align_s1, 's2_aligned': align_s2,
        'matches': match_count, 'matrix': score_matrix,
        'path': path_coords, 'dims': (m, n)
    }
    return results

class AlignmentApp:
    def __init__(self, root):
        self.root = root
        self.root.title("NW Alignment App (Python/Tkinter)")
        self.root.geometry("1200x700")
        
        style = ttk.Style()
        style.theme_use('clam')

        left_pane = ttk.Frame(root, padding="10")
        left_pane.pack(side=tk.LEFT, fill=tk.Y)
        
        right_pane = ttk.Frame(root, padding="10")
        right_pane.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        seq_frame = ttk.LabelFrame(left_pane, text="Sequences", padding="10")
        seq_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(seq_frame, text="Sq 1").grid(row=0, column=0)
        self.s1_entry = ttk.Entry(seq_frame, width=25, font=('Courier', 10))
        self.s1_entry.insert(0, "ACCGTGAAGCCAATAC")
        self.s1_entry.grid(row=0, column=1, padx=5)
        
        ttk.Label(seq_frame, text="Sq 2").grid(row=1, column=0)
        self.s2_entry = ttk.Entry(seq_frame, width=25, font=('Courier', 10))
        self.s2_entry.insert(0, "AGCGTGCAGCCAATAC")
        self.s2_entry.grid(row=1, column=1, padx=5)

        param_frame = ttk.LabelFrame(left_pane, text="Parameters", padding="10")
        param_frame.pack(fill=tk.X, pady=10)

        ttk.Label(param_frame, text="Gap =").grid(row=0, column=0, sticky=tk.E)
        self.gap_entry = ttk.Entry(param_frame, width=5, justify=tk.CENTER)
        self.gap_entry.insert(0, "0")
        self.gap_entry.grid(row=0, column=1)

        ttk.Label(param_frame, text="Mach =").grid(row=1, column=0, sticky=tk.E)
        self.match_entry = ttk.Entry(param_frame, width=5, justify=tk.CENTER)
        self.match_entry.insert(0, "1")
        self.match_entry.grid(row=1, column=1)

        ttk.Label(param_frame, text="MMach =").grid(row=2, column=0, sticky=tk.E)
        self.mismatch_entry = ttk.Entry(param_frame, width=5, justify=tk.CENTER)
        self.mismatch_entry.insert(0, "-1")
        self.mismatch_entry.grid(row=2, column=1)

        self.align_btn = ttk.Button(left_pane, text="Align", command=self.perform_alignment)
        self.align_btn.pack(fill=tk.X, pady=20)

        presets_frame = ttk.LabelFrame(left_pane, text="Presets (Visual placeholders)", padding="10")
        presets_frame.pack(fill=tk.X, side=tk.BOTTOM)
        for i in range(4):
            ttk.Button(presets_frame, text=f"Setting {i+1}").pack(fill=tk.X, pady=2)

        plots_container = ttk.Frame(right_pane)
        plots_container.pack(fill=tk.BOTH, expand=True)
        
        self.heatmap_frame = ttk.LabelFrame(plots_container, text="Graphic representation of alignment matrix")
        self.heatmap_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)
        self.heatmap_canvas = None

        self.traceback_frame = ttk.LabelFrame(plots_container, text="Traceback path")
        self.traceback_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5)
        self.traceback_canvas = None

        output_frame = ttk.LabelFrame(right_pane, text="Results Output")
        output_frame.pack(fill=tk.X, pady=10)
        self.result_text = scrolledtext.ScrolledText(output_frame, height=12, font=('Courier', 10))
        self.result_text.pack(fill=tk.BOTH)

    def perform_alignment(self):
        s1 = self.s1_entry.get()
        s2 = self.s2_entry.get()
        try:
            gap = int(self.gap_entry.get())
            match = int(self.match_entry.get())
            mismatch = int(self.mismatch_entry.get())
        except ValueError:
            self.result_text.delete('1.0', tk.END)
            self.result_text.insert(tk.END, "Error: Parameters must be integers.")
            return

        res = run_needleman_wunsch_core(s1, s2, gap, match, mismatch)

        visual_line = ""
        for k in range(len(res['s1_aligned'])):
            if res['s1_aligned'][k] == res['s2_aligned'][k]:
                visual_line += "|"
            else:
                visual_line += " "
        
        final_len = len(res['s1_aligned'])
        similarity = int((res['matches'] / final_len) * 100) if final_len > 0 else 0
        
        output_str = f"Show Alignment:\n\n"
        output_str += f"{res['s1_aligned']}\n"
        output_str += f"{visual_line}\n"
        output_str += f"{res['s2_aligned']}\n\n"
        output_str += f"Matches = {res['matches']}\n"
        output_str += f"Length = {final_len}\n"
        output_str += f"Similarity = {similarity} %\n"
        output_str += f"Tracing back: M[{res['dims'][0]},{res['dims'][1]}]"

        self.result_text.delete('1.0', tk.END)
        self.result_text.insert(tk.END, output_str)

        self.plot_heatmap(res['matrix'])
        self.plot_traceback_grid(res['dims'][0], res['dims'][1], res['path'])

    def plot_heatmap(self, matrix):
        if self.heatmap_canvas:
            self.heatmap_canvas.get_tk_widget().destroy()

        fig, ax = plt.subplots(figsize=(5, 4))
        cmap = mcolors.LinearSegmentedColormap.from_list("", ["#4a004a", "#ff0055"])
        cax = ax.imshow(matrix, cmap=cmap, interpolation='nearest', aspect='auto')
        ax.set_title("Scoring Matrix Heatmap", fontsize=10)
        ax.axis('off')
        plt.tight_layout()

        self.heatmap_canvas = FigureCanvasTkAgg(fig, master=self.heatmap_frame)
        self.heatmap_canvas.draw()
        self.heatmap_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def plot_traceback_grid(self, rows, cols, path_set):
        if self.traceback_canvas:
            self.traceback_canvas.get_tk_widget().destroy()

        fig, ax = plt.subplots(figsize=(5, 4))
        
        grid_data = [[0 for _ in range(cols + 1)] for _ in range(rows + 1)]
        for r in range(rows + 1):
            for c in range(cols + 1):
                if (r, c) in path_set:
                    grid_data[r][c] = 1

        cmap = mcolors.ListedColormap(['#ffffcc', '#cc0000'])
        
        ax.imshow(grid_data, cmap=cmap, interpolation='nearest', aspect='auto')
        
        ax.set_xticks([x - 0.5 for x in range(1, cols + 2)], minor=True)
        ax.set_yticks([y - 0.5 for y in range(1, rows + 2)], minor=True)
        ax.grid(which="minor", color="black", linestyle='-', linewidth=1)
        ax.tick_params(which="minor", bottom=False, left=False)
        ax.set_xticks([]); ax.set_yticks([])

        ax.set_title("Optimal Path (Red)", fontsize=10)
        plt.tight_layout()

        self.traceback_canvas = FigureCanvasTkAgg(fig, master=self.traceback_frame)
        self.traceback_canvas.draw()
        self.traceback_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

if __name__ == "__main__":
    root = tk.Tk()
    app = AlignmentApp(root)
    root.mainloop()
