import sys
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
    # Force the first forward primer
    first_fwd_primer = extend_primer(sequence, 0, tm_target, 'forward')
    first_fwd_primer_pos = 0
    if not first_fwd_primer:
        print("Error: Unable to find a valid first forward primer with the desired Tm.")
        return [], []

    # Force the last reverse primer
    last_rev_primer = extend_primer(sequence, len(sequence), tm_target, 'reverse')
    last_rev_primer_pos = len(sequence) - len(last_rev_primer)
    if not last_rev_primer:
        print("Error: Unable to find a valid last reverse primer with the desired Tm.")
        return [], []

    slice_positions = [i for i in range(tile_size, len(sequence), tile_size)]

    fwd_primers = [(first_fwd_primer, first_fwd_primer_pos)]
    rev_primers = []

    for pos in slice_positions:
        # Create reverse primer 35 bases downstream of the slice point
        rev_start = pos + 35
        rev_primer = extend_primer(sequence, rev_start, tm_target, 'reverse')
        if not rev_primer:
            print(f"Error: Unable to find a valid reverse primer at position {rev_start} with the desired Tm.")
            continue
        rev_primers.append((rev_primer, rev_start - len(rev_primer)))

        # Create forward primer 35 bases upstream of the slice point
        fwd_start = pos - 35
        fwd_primer = extend_primer(sequence, fwd_start, tm_target, 'forward')
        if not fwd_primer:
            print(f"Error: Unable to find a valid forward primer at position {fwd_start} with the desired Tm.")
            continue
        fwd_primers.append((fwd_primer, fwd_start))

    rev_primers.append((last_rev_primer, last_rev_primer_pos))

    # Check for overlaps and print primer details
    for i, (fwd, rev) in enumerate(zip(fwd_primers, rev_primers)):
        fwd_primer, fwd_pos = fwd
        rev_primer, rev_pos = rev
        if fwd_primer and rev_primer:
            if fwd_pos < rev_pos < (fwd_pos + len(fwd_primer)):
                print(f"Warning: Forward primer {fwd_primer} and reverse primer {rev_primer} overlap.")
            print(f"Tile {i + 1}:")
            print(f"  Forward primer: {fwd_primer} (start: {fwd_pos})")
            print(f"  Reverse primer: {rev_primer} (start: {rev_pos})")

    return fwd_primers, rev_primers


def main(input_fasta, output_fasta, tile_size, tm_target):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    primer_records = []

    for record in records:
        sequence = str(record.seq)
        fwd_primers, rev_primers = design_primers(sequence, tile_size, tm_target)

        for i, (fwd, rev) in enumerate(zip(fwd_primers, rev_primers)):
            fwd_primer, fwd_pos = fwd
            rev_primer, rev_pos = rev
            if fwd_primer and rev_primer:
                fwd_record = SeqRecord(Seq(fwd_primer), id=f"{record.id}_tile_{i + 1}_fwd",
                                       description="Forward primer")
                rev_record = SeqRecord(Seq(rev_primer), id=f"{record.id}_tile_{i + 1}_rev",
                                       description="Reverse primer")
                primer_records.append(fwd_record)
                primer_records.append(rev_record)

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(primer_records, output_handle, "fasta")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python primer_design.py <input_fasta> <output_fasta> <tile_size> <tm_target>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    tile_size = int(sys.argv[3])
    tm_target = float(sys.argv[4])

    main(input_fasta, output_fasta, tile_size, tm_target)
