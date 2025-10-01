import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from collections import Counter

def read_fasta(path):
    header = None
    chunks = []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)

def clean_sequence(seq: str) -> str:
    s = seq.upper()
    return "".join(ch for ch in s if ch.isalpha())

def percentages(seq: str):
    s = clean_sequence(seq)
    if not s:
        return [], 0
    cnt = Counter(s)
    total = sum(cnt.values())
    rows = [(sym, cnt[sym], (cnt[sym]/total)*100) for sym in sorted(cnt)]
    return rows, total

class FastaGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Percentages")
        self.geometry("720x480")

        top = ttk.Frame(self, padding=8)
        top.pack(side=tk.TOP, fill=tk.X)
        ttk.Button(top, text="Open FASTA…", command=self.open_fasta).pack(side=tk.LEFT)
        self.path_var = tk.StringVar(value="No file selected.")
        ttk.Label(top, textvariable=self.path_var).pack(side=tk.LEFT, padx=8)

        main = ttk.Frame(self, padding=8)
        main.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(main)
        left.pack(side=tk.LEFT, fill=tk.Y)
        ttk.Label(left, text="Sequences").pack(anchor="w")
        self.listbox = tk.Listbox(left, height=20, exportselection=False)
        self.listbox.pack(fill=tk.Y, expand=True)
        self.listbox.bind("<<ListboxSelect>>", self.on_select)

        right = ttk.Frame(main)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(12, 0))

        self.info_var = tk.StringVar(value="Open a FASTA file to begin.")
        ttk.Label(right, textvariable=self.info_var, font=("TkDefaultFont", 10, "bold")).pack(anchor="w", pady=(0,6))

        self.tree = ttk.Treeview(right, columns=("symbol", "count", "pct"), show="headings", height=16)
        self.tree.heading("symbol", text="Symbol")
        self.tree.heading("count", text="Count")
        self.tree.heading("pct", text="Percent")
        self.tree.column("symbol", width=100, anchor="center")
        self.tree.column("count", width=100, anchor="e")
        self.tree.column("pct", width=120, anchor="e")
        self.tree.pack(fill=tk.BOTH, expand=True)

        self.records = []  
    def open_fasta(self):
        path = filedialog.askopenfilename(
            title="Open FASTA file",
            filetypes=[("FASTA files", "*.fa *.fasta *.fna *.faa *.txt"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            recs = list(read_fasta(path))
        except Exception as e:
            messagebox.showerror("Error", f"Failed to read file:\n{e}")
            return
        if not recs:
            messagebox.showwarning("Warning", "No FASTA records found in file.")
            return

        self.path_var.set(path)
        self.records.clear()
        self.listbox.delete(0, tk.END)
        self.tree.delete(*self.tree.get_children())

        for idx, (header, seq) in enumerate(recs, start=1):
            rows, total = percentages(seq)
            item = {
                "header": header or f"record_{idx}",
                "rows": rows,
                "total": total
            }
            self.records.append(item)
            self.listbox.insert(tk.END, item["header"])

        self.info_var.set("Select a sequence on the left to view composition.")
        if self.records:
            self.listbox.selection_set(0)
            self.on_select(None)

    def on_select(self, _event):
        sel = self.listbox.curselection()
        if not sel:
            return
        i = sel[0]
        rec = self.records[i]

        self.info_var.set(f"{rec['header']}  |  Length: {rec['total']}")
        self.tree.delete(*self.tree.get_children())
        if rec["total"] == 0:
            self.tree.insert("", tk.END, values=("-", 0, "0.00%"))
            return
        for sym, cnt, pct in rec["rows"]:
            self.tree.insert("", tk.END, values=(sym, cnt, f"{pct:.2f}%"))

if __name__ == "__main__":
    app = FastaGUI()
    app.mainloop()
