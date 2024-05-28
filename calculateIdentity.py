import sys
from Bio import SeqIO

def calculate_percentage_identity(seq1, seq2, include_gaps=True):
    matches = 0
    length = len(seq1)
    for a, b in zip(seq1, seq2):
        if include_gaps:
            if a == b:
                matches += 1
        else:
            if a != '-' and b != '-':
                if a == b:
                    matches += 1
    if include_gaps:
        return (matches / length) * 100
    else:
        aligned_positions = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
        return (matches / aligned_positions) * 100 if aligned_positions > 0 else 0

def read_fasta_and_calculate_identity(fasta_file):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    if len(sequences) != 2:
        raise ValueError("FASTA file should contain exactly two sequences")
    seq1 = str(sequences[0].seq)
    seq2 = str(sequences[1].seq)
    
    identity_with_gaps = calculate_percentage_identity(seq1, seq2, include_gaps=True)
    identity_without_gaps = calculate_percentage_identity(seq1, seq2, include_gaps=False)
    
    return identity_with_gaps, identity_without_gaps

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python calculate_identity.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    
    try:
        identity_with_gaps, identity_without_gaps = read_fasta_and_calculate_identity(fasta_file)
        print(f"Percentage Identity (Including Gaps): {identity_with_gaps:.2f}%")
        print(f"Percentage Identity (Excluding Gaps): {identity_without_gaps:.2f}%")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

