#!/usr/bin/env python3
import argparse
import sys

def parse_fasta(filename):
    """
    Parses a FASTA file and returns a dictionary mapping the record ID (first word after '>') 
    to the concatenated sequence (ignoring newline breaks).
    """
    sequences = {}
    header = None
    seq_lines = []
    try:
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        # Save the previous record
                        sequences[header] = "".join(seq_lines)
                    # Use only the first designator (ignoring additional text)
                    header = line[1:].split()[0]
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if header is not None:
                sequences[header] = "".join(seq_lines)
    except Exception as e:
        print(f"Error reading '{filename}': {e}")
        sys.exit(1)
    return sequences

def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                  'N': 'N', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def align_oligo_to_ref(oligo, ref_seq):
    """
    Attempts to align the given oligo to the reference sequence.
    First checks for an exact match using str.find.
    If not found, scans the reference for any substring (of the same length as oligo)
    that has exactly one mismatch.
    Returns a tuple (position, mismatch_count) if found, otherwise (None, None).
    """
    # Try exact match
    pos = ref_seq.find(oligo)
    if pos != -1:
        return pos, 0
    # Try approximate match: allow exactly one mismatch
    oligo_len = len(oligo)
    for i in range(len(ref_seq) - oligo_len + 1):
        segment = ref_seq[i:i+oligo_len]
        mismatches = sum(1 for a, b in zip(oligo, segment) if a != b)
        if mismatches == 1:
            return i, 1
    return None, None

def create_alignment_string(oligo, ref_seq, pos):
    """
    Creates an alignment string of the same length as the reference.
    It consists entirely of '-' (gaps) except where the oligo sequence is inserted at the given position.
    """
    alignment = list("-" * len(ref_seq))
    for j, base in enumerate(oligo):
        alignment[pos+j] = base
    return "".join(alignment)

def main():
    parser = argparse.ArgumentParser(description="Align oligos to a reference sequence.")
    parser.add_argument("-ref", required=True, help="Reference sequence file in FASTA format")
    parser.add_argument("-oligo", required=True, help="Oligo sequences file in FASTA format")
    parser.add_argument("-output", required=True, help="Output file for aligned sequences in FASTA format")
    args = parser.parse_args()

    # Load reference; expect one record
    ref_records = parse_fasta(args.ref)
    if not ref_records:
        print("Error: Reference file is empty or improperly formatted.")
        sys.exit(1)
    # Use the first record as the reference sequence
    ref_id = list(ref_records.keys())[0]
    ref_seq = ref_records[ref_id]

    # Load oligo sequences
    oligo_records = parse_fasta(args.oligo)
    if not oligo_records:
        print("Error: Oligo file is empty or improperly formatted.")
        sys.exit(1)
    
    total_oligos = len(oligo_records)
    success_count = 0
    rc_count = 0
    mismatch_count = 0
    
    # Dictionaries to store alignment results
    aligned_results = {}   # header -> aligned string
    mismatch_results = {}  # for oligos aligned with exactly one mismatch

    # Process each oligo record
    for oligo_id, oligo_seq in oligo_records.items():
        candidate = oligo_seq
        used_rc = False
        pos, diff = align_oligo_to_ref(candidate, ref_seq)
        
        # If not found in original orientation, try reverse complement
        if pos is None:
            candidate_rc = reverse_complement(oligo_seq)
            pos, diff = align_oligo_to_ref(candidate_rc, ref_seq)
            if pos is not None:
                used_rc = True
                candidate = candidate_rc
        
        # If an alignment was found (either exact or with one mismatch)
        if pos is not None:
            # Create alignment string (with dashes as gaps)
            aligned_str = create_alignment_string(candidate, ref_seq, pos)
            header = oligo_id + ("-RC" if used_rc else "")
            
            aligned_results[header] = aligned_str
            success_count += 1
            if used_rc:
                rc_count += 1
            if diff == 1:
                mismatch_count += 1
                mismatch_results[header] = aligned_str
        # If no alignment was found in either orientation, omit this oligo from results

    # Write FASTA formatted output
    try:
        with open(args.output, "w") as outf:
            # Write reference record (unaltered)
            outf.write(f">reference\n{ref_seq}\n")
            # Write each aligned oligo record
            for header, alignment in aligned_results.items():
                outf.write(f">{header}\n{alignment}\n")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)
    
    # Write mismatch error file (only oligos that aligned with one mismatch)
    if mismatch_count > 0:
        try:
            with open("mismatch.err", "w") as errf:
                for header, alignment in mismatch_results.items():
                    errf.write(f">{header}\n{alignment}\n")
        except Exception as e:
            print(f"Error writing mismatch error file: {e}")
            sys.exit(1)

    # Print summary messages
    print(f"{success_count} out of {total_oligos} oligos successfully aligned.")
    print(f"{rc_count} out of {total_oligos} oligos were reverse-complemented and successfully aligned.")
    if mismatch_count > 0:
        print(f"{mismatch_count} out of {total_oligos} oligos is/are a mismatch. See 'mismatch.err' for details.")

if __name__ == "__main__":
    main()