import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from ttkbootstrap.dialogs import Messagebox
from tkinter import filedialog
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import threading


def calculate_tm(sequence, Na=50, Mg=1.5):
    return mt.Tm_NN(sequence, Na=Na, Mg=Mg)


def extend_primer(sequence, start, tm_target, direction='forward', Na=50, Mg=1.5, min_length=20):
    if direction == 'forward':
        for end in range(start + min_length, len(sequence)):
            primer = sequence[start:end]
            if calculate_tm(primer, Na, Mg) >= tm_target:
                return primer
    elif direction == 'reverse':
        for end in range(start, min_length, -1):
            primer = sequence[end:start]
            if len(primer) >= min_length and calculate_tm(primer, Na, Mg) >= tm_target:
                return str(Seq(primer).reverse_complement())
    return None


def design_primers(sequence, tile_size, tm_target, overlap=35):
    first_fwd_primer = extend_primer(sequence, 0, tm_target, 'forward')
    first_fwd_primer_pos = 0
    if not first_fwd_primer:
        Messagebox.show_error("Unable to find a valid first forward primer with the desired Tm.", "Error")
        return [], []

    last_rev_primer = extend_primer(sequence, len(sequence), tm_target, 'reverse')
    last_rev_primer_pos = len(sequence) - len(last_rev_primer)
    if not last_rev_primer:
        Messagebox.show_error("Unable to find a valid last reverse primer with the desired Tm.", "Error")
        return [], []

    slice_positions = [i for i in range(tile_size, len(sequence), tile_size)]

    fwd_primers = [(first_fwd_primer, first_fwd_primer_pos)]
    rev_primers = []

    for pos in slice_positions:
        rev_start = pos + 35
        rev_primer = extend_primer(sequence, rev_start, tm_target, 'reverse')
        if not rev_primer:
            Messagebox.show_error(f"Unable to find a valid reverse primer at position {rev_start} with the desired Tm.",
                                  "Error")
            continue
        rev_primers.append((rev_primer, rev_start - len(rev_primer)))

        fwd_start = pos - 35
        fwd_primer = extend_primer(sequence, fwd_start, tm_target, 'forward')
        if not fwd_primer:
            Messagebox.show_error(f"Unable to find a valid forward primer at position {fwd_start} with the desired Tm.",
                                  "Error")
            continue
        fwd_primers.append((fwd_primer, fwd_start))

    rev_primers.append((last_rev_primer, last_rev_primer_pos))

    for i, (fwd, rev) in enumerate(zip(fwd_primers, rev_primers)):
        fwd_primer, fwd_pos = fwd
        rev_primer, rev_pos = rev
        if fwd_primer and rev_primer:
            if fwd_pos < rev_pos < (fwd_pos + len(fwd_primer)):
                Messagebox.show_warning(f"Forward primer {fwd_primer} and reverse primer {rev_primer} overlap.",
                                        "Warning")
            print(f"Tile {i + 1}:")
            print(f"  Forward primer: {fwd_primer} (start: {fwd_pos})")
            print(f"  Reverse primer: {rev_primer} (start: {rev_pos})")

    return fwd_primers, rev_primers


def main(input_fasta, output_fasta, tile_size, tm_target, log_text_widget):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    output_records = []

    for record in records:
        sequence = str(record.seq)
        fwd_primers, rev_primers = design_primers(sequence, tile_size, tm_target)

        for i, (fwd, rev) in enumerate(zip(fwd_primers, rev_primers)):
            fwd_primer, fwd_pos = fwd
            rev_primer, rev_pos = rev
            if fwd_primer and rev_primer:
                fwd_record = SeqRecord(Seq(fwd_primer), id=f"{record.id}_fwd_primer_{i + 1}",
                                       description="Forward Primer")
                rev_record = SeqRecord(Seq(rev_primer), id=f"{record.id}_rev_primer_{i + 1}",
                                       description="Reverse Primer")
                output_records.append(fwd_record)
                output_records.append(rev_record)
                log_text_widget.insert('end', f"Tile {i + 1}:\n")
                log_text_widget.insert('end', f"  Forward primer: {fwd_primer} (start: {fwd_pos})\n")
                log_text_widget.insert('end', f"  Reverse primer: {rev_primer} (start: {rev_pos})\n")

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(output_records, output_handle, "fasta")

    log_text_widget.insert('end', "Primer design completed and saved to output file.\n")
    log_text_widget.see('end')


class App:
    def __init__(self, root):
        self.root = root
        self.root.title("TileCreator - Primer Designer")
        self.style = ttk.Style()
        self.style.theme_use('lumen')

        self.frame = ttk.Frame(self.root, padding="10")
        self.frame.grid(row=0, column=0, sticky=(W, E, N, S))

        self.lbl_input = ttk.Label(self.frame, text="Input FASTA File:")
        self.lbl_input.grid(row=0, column=0, padx=10, pady=10, sticky=W)
        self.entry_input = ttk.Entry(self.frame, width=50)
        self.entry_input.grid(row=0, column=1, padx=10, pady=10)
        self.btn_input = ttk.Button(self.frame, text="Browse", command=self.browse_input)
        self.btn_input.grid(row=0, column=2, padx=10, pady=10)

        self.lbl_output = ttk.Label(self.frame, text="Output File:")
        self.lbl_output.grid(row=1, column=0, padx=10, pady=10, sticky=W)
        self.entry_output = ttk.Entry(self.frame, width=50)
        self.entry_output.grid(row=1, column=1, padx=10, pady=10)
        self.btn_output = ttk.Button(self.frame, text="Browse", command=self.browse_output)
        self.btn_output.grid(row=1, column=2, padx=10, pady=10)

        self.lbl_tile_size = ttk.Label(self.frame, text="Tile Size:")
        self.lbl_tile_size.grid(row=2, column=0, padx=10, pady=10, sticky=W)
        self.entry_tile_size = ttk.Entry(self.frame, width=20)
        self.entry_tile_size.grid(row=2, column=1, padx=10, pady=10, sticky=W)

        self.lbl_tm_target = ttk.Label(self.frame, text="Primer Tm Target:")
        self.lbl_tm_target.grid(row=3, column=0, padx=10, pady=10, sticky=W)
        self.entry_tm_target = ttk.Entry(self.frame, width=20)
        self.entry_tm_target.grid(row=3, column=1, padx=10, pady=10, sticky=W)

        self.btn_run = ttk.Button(self.frame, text="Run", command=self.run_olaygo)
        self.btn_run.grid(row=4, column=0, columnspan=3, padx=10, pady=20)

        self.log_text = ttk.Text(self.frame, width=70, height=10)
        self.log_text.grid(row=5, column=0, columnspan=3, padx=10, pady=10)

    def browse_input(self):
        self.input_file = filedialog.askopenfilename(
            filetypes=[("FASTA files", "*.fasta"), ("FA files", "*.fa"), ("FNA files", "*.fna"), ("All files", "*.*")])
        self.entry_input.delete(0, 'end')
        self.entry_input.insert(0, self.input_file)

    def browse_output(self):
        self.output_file = filedialog.asksaveasfilename(defaultextension=".fasta",
                                                        filetypes=[("FASTA files", "*.fasta"), ("FA files", "*.fa"),
                                                                   ("FNA files", "*.fna"), ("All files", "*.*")])
        self.entry_output.delete(0, 'end')
        self.entry_output.insert(0, self.output_file)

    def run_olaygo(self):
        input_file = self.entry_input.get()
        output_file = self.entry_output.get()
        tile_size = int(self.entry_tile_size.get())
        tm_target = float(self.entry_tm_target.get())

        threading.Thread(target=self.run_main, args=(input_file, output_file, tile_size, tm_target)).start()

    def run_main(self, input_file, output_file, tile_size, tm_target):
        main(input_file, output_file, tile_size, tm_target, self.log_text)
        Messagebox.show_info("Primer design completed and saved to output file.", "Success")


if __name__ == "__main__":
    root = ttk.Window(themename="lumen")
    app = App(root)
    root.mainloop()
