# Optimized CKSAAP (Composition of K-Spaced Amino Acid Pairs) feature script
from itertools import product
from collections import Counter

AMINO_ACIDS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']


DIPEPTIDES = [a + b for a, b in product(AMINO_ACIDS, repeat=2)]

def calculate_cksaap(sequence, k=0):
    """Calculate CKSAAP feature vector with specified spacing k."""
    pairs = [sequence[i] + sequence[i + k + 1] for i in range(len(sequence) - k - 1)]
    counts = Counter(pairs)
    return [counts[p] for p in DIPEPTIDES]

def run_cksaap_demo():
    examples = [
        ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQKSTVWY"), 
        ("MKTIIALSYIFCLVFA", "MKTIIALSYVFCLVFA")           
    ]

    for wt_seq, mut_seq in examples:
        print(f"\nWild-type sequence: {wt_seq}")
        print(f"Mutant sequence    : {mut_seq}")

        for k in [0, 1, 2]:  # CKSAAP with k = 0, 1, 2
            features = calculate_cksaap(mut_seq, k=k)
            print(f"  CKSAAP-k={k} feature vector length: {len(features)}")
            print(f"  First 10 values: {features[:10]}")  # Print first 10 as example


run_cksaap_demo()
