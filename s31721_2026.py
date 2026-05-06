# ============================================================
# Album number : s12345
# Date         : 2026-04-30
# Description  : Random DNA sequence generator in FASTA format.
#                Includes validation, statistics, name embedding,
#                and additional bioinformatics features.
# ============================================================

import random
import csv
import os


# ------------------------------------------------------------
# CORE GENERATION
# ------------------------------------------------------------

def generate_sequence(length: int) -> str:
    """Returns a random DNA sequence of the specified length.

    Uses random.choices for uniform nucleotide distribution.
    Each position is independently drawn from {A, C, G, T}.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choices(nucleotides, k=length))


def generate_sequence_weighted(length: int, weights: dict) -> str:
    """Returns a random DNA sequence with user-defined nucleotide distribution.

    Args:
        length  : number of nucleotides to generate
        weights : dict with keys A, C, G, T and float percentage values (sum=100)

    Uses random.choices with custom weights so each nucleotide appears
    roughly at the requested frequency over long sequences.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    w = [weights[n] for n in nucleotides]
    return ''.join(random.choices(nucleotides, weights=w, k=length))


# ------------------------------------------------------------
# STATISTICS
# ------------------------------------------------------------

def calculate_stats(sequence: str) -> dict:
    """Returns a dictionary of sequence statistics.

    Keys: 'A', 'C', 'G', 'T'  (float, percentage of each nucleotide),
          'GC'                  (float, GC-content percentage).

    The input sequence must contain only uppercase nucleotide characters.
    Lowercase characters (e.g. embedded names) should be stripped before
    calling this function.
    """
    n = len(sequence)
    if n == 0:
        return {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0, 'GC': 0.0}

    counts = {nuc: sequence.count(nuc) for nuc in 'ACGT'}
    stats = {nuc: round(counts[nuc] / n * 100, 2) for nuc in 'ACGT'}
    gc = round((counts['G'] + counts['C']) / n * 100, 2)
    stats['GC'] = gc
    return stats


def print_stats(stats: dict, length: int) -> None:
    """Prints sequence statistics to stdout in the expected format."""
    print(f"\nSequence statistics (n={length}):")
    for nuc in ['A', 'C', 'G', 'T']:
        print(f"  {nuc}: {stats[nuc]:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")


# ------------------------------------------------------------
# NAME EMBEDDING
# ------------------------------------------------------------

def insert_name(sequence: str, name: str) -> str:
    """Inserts a name at a random position in the sequence.

    The name is written in lowercase so it is visually distinguishable
    from the uppercase nucleotides. The inserted letters do NOT count
    as part of the biological sequence — statistics must be calculated
    on the original sequence before calling this function.

    Args:
        sequence : uppercase DNA string
        name     : user's name (any characters accepted)

    Returns a new string with the lowercase name spliced in at a random index.
    """
    if not name:
        return sequence

    pos = random.randint(0, len(sequence))
    return sequence[:pos] + name.lower() + sequence[pos:]


# ------------------------------------------------------------
# FASTA FORMATTING
# ------------------------------------------------------------

def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """Returns a formatted FASTA record as a string.

    The header line starts with '>' followed by seq_id. If description
    is non-empty, a space and the description are appended.
    The sequence is then wrapped at line_width characters per line.

    Args:
        seq_id      : identifier without whitespace
        description : optional free-text description
        sequence    : nucleotide string (may contain embedded name chars)
        line_width  : characters per sequence line (default 80)
    """
    if description.strip():
        header = f">{seq_id} {description}"
    else:
        header = f">{seq_id}"

    # Break sequence into fixed-width lines
    lines = [sequence[i:i + line_width] for i in range(0, len(sequence), line_width)]
    return header + '\n' + '\n'.join(lines)


def save_fasta(filepath: str, fasta_content: str) -> None:
    """Writes a FASTA record to a file, appending a trailing newline.

    Appends rather than overwrites so multi-FASTA batch files can be
    built incrementally by calling this function multiple times.
    """
    with open(filepath, 'a', encoding='utf-8') as fh:
        fh.write(fasta_content + '\n')


# ------------------------------------------------------------
# INPUT VALIDATION
# ------------------------------------------------------------

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """Gets an integer from the user in the range [min_val, max_val].

    Repeats the prompt on invalid input instead of raising an exception.
    Handles non-numeric input and out-of-range values separately so the
    error message is always accurate.
    """
    while True:
        raw = input(prompt).strip()
        try:
            value = int(raw)
        except ValueError:
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")
            continue

        if not (min_val <= value <= max_val):
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")
            continue

        return value


def validate_seq_id(prompt: str) -> str:
    """Gets a sequence ID from the user, rejecting any input containing whitespace.

    FASTA identifiers must not contain spaces or tab characters because
    the space character is the separator between ID and description in
    the header line.
    """
    while True:
        seq_id = input(prompt).strip()
        if not seq_id:
            print("Error: ID cannot be empty.")
            continue
        if any(ch.isspace() for ch in seq_id):
            print("Error: ID cannot contain whitespace characters.")
            continue
        return seq_id


# ------------------------------------------------------------
# ADDITIONAL FEATURE 1 — BATCH MODE
# ------------------------------------------------------------

def batch_mode() -> None:
    """Generates multiple sequences and saves them all to a single multi-FASTA file.

    The user specifies the number of sequences and a common base ID.
    Each sequence receives a unique numeric suffix (e.g. Seq_001, Seq_002).
    All records are written to <base_id>_batch.fasta.
    """
    count = validate_positive_int("How many sequences to generate? ", min_val=1, max_val=1000)
    base_id = validate_seq_id("Enter base ID (e.g. Seq): ")
    length = validate_positive_int("Enter sequence length: ")
    description = input("Enter description (optional): ").strip()
    name = input("Enter your name (optional, leave blank to skip): ").strip()

    filepath = f"{base_id}_batch.fasta"

    # Remove any existing file so we start fresh
    if os.path.exists(filepath):
        os.remove(filepath)

    for i in range(1, count + 1):
        seq_id = f"{base_id}_{i:03d}"
        seq = generate_sequence(length)

        # Statistics on pure nucleotide sequence (before name insertion)
        stats = calculate_stats(seq)

        if name:
            seq_with_name = insert_name(seq, name)
        else:
            seq_with_name = seq

        fasta_str = format_fasta(seq_id, description, seq_with_name)
        save_fasta(filepath, fasta_str)

        print(f"  [{i}/{count}] {seq_id} saved — GC: {stats['GC']:.2f}%")

    print(f"\nBatch complete. All sequences saved to: {filepath}")


# ------------------------------------------------------------
# ADDITIONAL FEATURE 2 — CONFIGURABLE NUCLEOTIDE DISTRIBUTION
# ------------------------------------------------------------

def get_custom_weights() -> dict:
    """Asks the user for custom nucleotide percentages that sum to 100.

    Repeats the prompt if the percentages are invalid (non-numeric,
    negative values, or sum not equal to 100).

    Returns a dict with keys A, C, G, T and float values.
    """
    while True:
        print("Enter nucleotide percentages (must sum to 100):")
        try:
            a = float(input("  A (%): "))
            c = float(input("  C (%): "))
            g = float(input("  G (%): "))
            t = float(input("  T (%): "))
        except ValueError:
            print("Error: please enter numeric values.")
            continue

        if any(v < 0 for v in [a, c, g, t]):
            print("Error: percentages cannot be negative.")
            continue

        total = a + c + g + t
        if abs(total - 100.0) > 0.01:
            print(f"Error: percentages must sum to 100 (got {total:.2f}).")
            continue

        return {'A': a, 'C': c, 'G': g, 'T': t}


# ------------------------------------------------------------
# ADDITIONAL FEATURE 3 — MOTIF SEARCH
# ------------------------------------------------------------

def find_motif(sequence: str, motif: str) -> list:
    """Searches for all occurrences of a motif in a DNA sequence.

    Positions are reported in 1-based biological indexing.
    Overlapping occurrences are detected (the search advances by 1 each time).

    Args:
        sequence : uppercase DNA string
        motif    : pattern to search for (uppercase)

    Returns a list of 1-based start positions.
    """
    positions = []
    start = 0
    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        positions.append(pos + 1)  # convert to 1-based
        start = pos + 1
    return positions


def motif_search_interactive(sequence: str) -> None:
    """Prompts the user for a motif and prints all found positions."""
    motif = input("Enter motif to search for (e.g. ATG): ").strip().upper()
    if not motif:
        print("No motif entered, skipping search.")
        return

    positions = find_motif(sequence, motif)
    if positions:
        print(f"Motif '{motif}' found {len(positions)} time(s) at position(s): {positions}")
    else:
        print(f"Motif '{motif}' not found in the sequence.")


# ------------------------------------------------------------
# ADDITIONAL FEATURE 4 — COMPLEMENT AND REVERSE COMPLEMENT
# ------------------------------------------------------------

def complement(sequence: str) -> str:
    """Returns the complementary strand of a DNA sequence (5'->3' direction).

    Uses standard Watson-Crick base pairing: A<->T, C<->G.
    Non-standard characters are left unchanged.
    """
    table = str.maketrans('ACGTacgt', 'TGCAtgca')
    return sequence.translate(table)


def reverse_complement(sequence: str) -> str:
    """Returns the reverse complement of a DNA sequence.

    This represents the antiparallel complementary strand read 5'->3'.
    Biologically, this is the sequence of the template strand oriented
    in the conventional left-to-right direction.
    """
    return complement(sequence)[::-1]


def add_complement_records(filepath: str, seq_id: str, sequence: str) -> None:
    """Appends complementary and reverse complementary records to a FASTA file.

    The complementary strand is tagged with suffix '_comp' and the
    reverse complement with '_revcomp'.
    """
    comp_seq = complement(sequence)
    revcomp_seq = reverse_complement(sequence)

    save_fasta(filepath, format_fasta(f"{seq_id}_comp", "Complementary strand", comp_seq))
    save_fasta(filepath, format_fasta(f"{seq_id}_revcomp", "Reverse complement strand", revcomp_seq))
    print(f"Complementary and reverse complement records added to {filepath}.")


# ------------------------------------------------------------
# ADDITIONAL FEATURE 5 — IN SILICO TRANSCRIPTION
# ------------------------------------------------------------

def transcribe_to_mrna(sequence: str) -> str:
    """Returns the mRNA sequence produced by in silico transcription.

    Transcription replaces each thymine (T) with uracil (U).
    The rest of the sequence is unchanged.
    """
    return sequence.replace('T', 'U').replace('t', 'u')


def add_mrna_record(filepath: str, seq_id: str, sequence: str) -> None:
    """Transcribes the sequence and appends the mRNA record to the FASTA file."""
    mrna = transcribe_to_mrna(sequence)
    save_fasta(filepath, format_fasta(f"{seq_id}_mRNA", "mRNA transcript (T->U)", mrna))
    print(f"mRNA transcript record added to {filepath}.")


# ------------------------------------------------------------
# ADDITIONAL FEATURE 6 — SLIDING WINDOW GC ANALYSIS
# ------------------------------------------------------------

def sliding_window_gc(sequence: str, window_size: int, step: int = 1) -> list:
    """Calculates GC-content in a sliding window across the sequence.

    Args:
        sequence    : uppercase DNA string
        window_size : width of the window in nucleotides
        step        : advance per iteration (default 1)

    Returns a list of (start_position, gc_content) tuples where
    start_position is 1-based and gc_content is a float percentage.
    """
    results = []
    n = len(sequence)

    for start in range(0, n - window_size + 1, step):
        window = sequence[start:start + window_size]
        gc_count = window.count('G') + window.count('C')
        gc_pct = round(gc_count / window_size * 100, 2)
        results.append((start + 1, gc_pct))  # 1-based position

    return results


def save_gc_csv(results: list, filepath: str) -> None:
    """Saves sliding window GC results to a CSV file.

    Columns: start_position, gc_content
    """
    with open(filepath, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(['start_position', 'gc_content'])
        writer.writerows(results)
    print(f"Sliding window GC data saved to: {filepath}")


def sliding_window_interactive(sequence: str, seq_id: str) -> None:
    """Prompts user for window size and step, runs analysis, saves CSV."""
    window_size = validate_positive_int(
        "Window size for GC analysis: ", min_val=1, max_val=len(sequence)
    )
    step = validate_positive_int("Step size: ", min_val=1, max_val=len(sequence))

    results = sliding_window_gc(sequence, window_size, step)
    csv_path = f"{seq_id}_gc_sliding.csv"
    save_gc_csv(results, csv_path)


# ------------------------------------------------------------
# MAIN FLOW
# ------------------------------------------------------------

def main():
    """Main program flow.

    The primary interaction matches the expected flow from the assignment:
      1. Enter sequence length
      2. Enter sequence ID
      3. Enter description
      4. Enter name
      5. Sequence is saved and statistics are printed
      6. Additional features are offered one by one

    Batch mode and custom nucleotide distribution are offered only AFTER
    the core sequence has been generated and saved, so the opening prompts
    stay identical to the expected interaction example.
    """
    # ---- Core flow (matches assignment interaction example) ----

    length = validate_positive_int("Enter sequence length: ")
    seq_id = validate_seq_id("Enter sequence ID: ")
    description = input("Enter a description of the sequence: ").strip()
    name = input("Enter your name: ").strip()

    # Generate uniform random sequence
    sequence = generate_sequence(length)

    # Statistics on the pure nucleotide sequence — BEFORE name insertion
    stats = calculate_stats(sequence)

    # Embed name at a random position (lowercase, visual only)
    sequence_with_name = insert_name(sequence, name)

    # Write FASTA file (remove stale file from a previous run first)
    filepath = f"{seq_id}.fasta"
    if os.path.exists(filepath):
        os.remove(filepath)

    fasta_record = format_fasta(seq_id, description, sequence_with_name)
    save_fasta(filepath, fasta_record)

    print(f"\nSequence saved to file: {filepath}")
    print_stats(stats, length)

    # ---- Additional features (offered after core save) --------

    print("\nAdditional features:")

    # Feature: Motif search
    do_motif = input("Search for a motif? [y/N]: ").strip().lower()
    if do_motif == 'y':
        motif_search_interactive(sequence)

    # Feature: Complement and reverse complement
    do_comp = input("Add complement/reverse complement records? [y/N]: ").strip().lower()
    if do_comp == 'y':
        add_complement_records(filepath, seq_id, sequence)

    # Feature: mRNA transcription
    do_mrna = input("Add mRNA transcript record? [y/N]: ").strip().lower()
    if do_mrna == 'y':
        add_mrna_record(filepath, seq_id, sequence)

    # Feature: Sliding window GC analysis (only meaningful for sequences >=2 nt)
    if length >= 2:
        do_gc = input("Run sliding window GC analysis? [y/N]: ").strip().lower()
        if do_gc == 'y':
            sliding_window_interactive(sequence, seq_id)

    # Feature: Custom nucleotide distribution — generate an extra sequence
    do_weighted = input("Generate an additional sequence with custom nucleotide distribution? [y/N]: ").strip().lower()
    if do_weighted == 'y':
        weights = get_custom_weights()
        w_length = validate_positive_int("Enter length for weighted sequence: ")
        w_seq_id = validate_seq_id("Enter ID for weighted sequence: ")
        w_desc = input("Enter description (optional): ").strip()
        w_seq = generate_sequence_weighted(w_length, weights)
        w_stats = calculate_stats(w_seq)
        w_seq_named = insert_name(w_seq, name)
        w_filepath = f"{w_seq_id}.fasta"
        if os.path.exists(w_filepath):
            os.remove(w_filepath)
        save_fasta(w_filepath, format_fasta(w_seq_id, w_desc, w_seq_named))
        print(f"Weighted sequence saved to: {w_filepath}")
        print_stats(w_stats, w_length)

    # Feature: Batch mode — generate many sequences into one multi-FASTA file
    do_batch = input("Run batch mode (generate multiple sequences)? [y/N]: ").strip().lower()
    if do_batch == 'y':
        batch_mode()

    print("\nDone.")


if __name__ == "__main__":
    main()
