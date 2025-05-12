

AA_PROPERTIES = {
    'accessible_area': {'A': 115, 'R': 225, 'N': 160, 'D': 150, 'C': 135, 'Q': 180, 'E': 190,
                        'G': 75, 'H': 195, 'I': 175, 'L': 170, 'K': 200, 'M': 185, 'F': 210,
                        'P': 145, 'S': 115, 'T': 140, 'W': 255, 'Y': 230, 'V': 155, '*': 0},
    'flexibility': {'A': 0.357, 'R': 0.529, 'N': 0.463, 'D': 0.511, 'C': 0.346, 'Q': 0.493,
                    'E': 0.497, 'G': 0.544, 'H': 0.323, 'I': 0.462, 'L': 0.365, 'K': 0.466,
                    'M': 0.295, 'F': 0.314, 'P': 0.509, 'S': 0.507, 'T': 0.444, 'W': 0.305,
                    'Y': 0.420, 'V': 0.386, '*': 0},
    'frequency': {'A': 0.077, 'R': 0.051, 'N': 0.043, 'D': 0.052, 'C': 0.020, 'Q': 0.041,
                  'E': 0.062, 'G': 0.074, 'H': 0.023, 'I': 0.053, 'L': 0.091, 'K': 0.059,
                  'M': 0.024, 'F': 0.040, 'P': 0.051, 'S': 0.069, 'T': 0.059, 'W': 0.014,
                  'Y': 0.032, 'V': 0.066, '*': 0},
    'hydrophobicity': {'A': 0.02, 'R': -0.42, 'N': -0.77, 'D': -1.04, 'C': 0.77, 'Q': -1.10,
                       'E': -1.14, 'G': -0.80, 'H': 0.26, 'I': 1.81, 'L': 1.14, 'K': -0.41,
                       'M': 1.00, 'F': 1.35, 'P': -0.09, 'S': -0.97, 'T': -0.77, 'W': 1.71,
                       'Y': 1.11, 'V': 1.13, '*': 0},
    'polarizability': {'A': 0.046, 'R': 0.291, 'N': 0.134, 'D': 0.105, 'C': 0.128, 'Q': 0.180,
                       'E': 0.151, 'G': 0.000, 'H': 0.230, 'I': 0.186, 'L': 0.186, 'K': 0.219,
                       'M': 0.221, 'F': 0.290, 'P': 0.131, 'S': 0.062, 'T': 0.108, 'W': 0.409,
                       'Y': 0.298, 'V': 0.140, '*': 0},
    'volume': {'A': 91.5, 'R': 196.1, 'N': 138.3, 'D': 135.2, 'C': 114.4, 'Q': 156.4,
               'E': 154.6, 'G': 67.5, 'H': 163.2, 'I': 162.6, 'L': 163.4, 'K': 162.5,
               'M': 165.9, 'F': 198.8, 'P': 123.4, 'S': 102.0, 'T': 126.0, 'W': 209.8,
               'Y': 237.2, 'V': 138.4, '*': 0},
    'side_chain': {'A': 2.87, 'R': 7.82, 'N': 4.58, 'D': 4.74, 'C': 4.47, 'Q': 6.11,
                   'E': 5.97, 'G': 2.06, 'H': 5.23, 'I': 4.92, 'L': 4.92, 'K': 6.89,
                   'M': 6.36, 'F': 4.62, 'P': 4.11, 'S': 3.97, 'T': 4.11, 'W': 7.68,
                   'Y': 4.73, 'V': 4.11, '*': 0}
}

def compute_feature(sequence, prop):
    values = [AA_PROPERTIES[prop].get(aa, 0) for aa in sequence]
    score = sum(values[i] * values[j] for i in range(len(values)) for j in range(i + 1, len(values)))
    return round(score / 10, 4)

def compute_all_sequence_features(sequence):
    return {name: compute_feature(sequence, name) for name in AA_PROPERTIES}


def net_volume(w, m): return AA_PROPERTIES['volume'].get(m, 0) - AA_PROPERTIES['volume'].get(w, 0)
def net_hydrophobicity(w, m): return AA_PROPERTIES['hydrophobicity'].get(m, 0) - AA_PROPERTIES['hydrophobicity'].get(w, 0)
def net_flexibility(w, m): return AA_PROPERTIES['flexibility'].get(m, 0) - AA_PROPERTIES['flexibility'].get(w, 0)

def mutation_hydrophobicity(w, m):
    def group(a): return 0 if a in 'ACFILMVW' else 1 if a in 'GHPSTY' else 2
    return group(w) * 3 + group(m)

def mutation_polarity(w, m):
    def group(a): return 0 if a in 'HKR' else 1 if a in 'ACFGILMPVW' else 2 if a in 'NQSTY' else 3
    return group(w) * 4 + group(m)

def size(w, m):
    def group(a): return 0 if a in 'CDNPT' else 1 if a in 'EHQV' else 2 if a in 'IKLMR' else 3 if a in 'FWY' else 4
    return group(w) * 5 + group(m)

def hydrogen_bond(w, m):
    def group(a): return 0 if a in 'ACFGILMPV' else 1 if a in 'KRW' else 2 if a in 'HNQSTY' else 3
    return group(w) * 4 + group(m)

def chemical_property(w, m):
    def group(a): return 0 if a in 'HKR' else 1 if a in 'NQ' else 2 if a in 'DE' else 3 if a in 'CM' else 4 if a in 'ST' else 5 if a in 'FWY' else 6
    return group(w) * 7 + group(m)

def mutation_type(w, m):
    aa_list = list(AA_PROPERTIES['volume'].keys())
    return aa_list.index(w) * 21 + aa_list.index(m)


examples = [
    ("MKTIIALSYIFCLVFA", "MKTIIALSYVFCLVFA"),  
    ("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQKSTVWY") 
]

results = []

for wt_seq, mut_seq in examples:
    for i, (w, m) in enumerate(zip(wt_seq, mut_seq)):
        if w != m:
            mutation_features = [
                net_volume(w, m), net_hydrophobicity(w, m), net_flexibility(w, m),
                mutation_hydrophobicity(w, m), mutation_polarity(w, m), mutation_type(w, m),
                size(w, m), hydrogen_bond(w, m), chemical_property(w, m)
            ]
            sequence_features = compute_all_sequence_features(mut_seq)
            combined = {**sequence_features, **{
                "net_volume": mutation_features[0],
                "net_hydrophobicity": mutation_features[1],
                "net_flexibility": mutation_features[2],
                "mutation_hydrophobicity": mutation_features[3],
                "mutation_polarity": mutation_features[4],
                "mutation_type": mutation_features[5],
                "size": mutation_features[6],
                "hydrogen_bond": mutation_features[7],
                "chemical_property": mutation_features[8]
            }}
            results.append((mut_seq, combined))
            break  # one mutation per sample

for seq, feats in results:
    print(f"\nMutant Sequence: {seq}")
    for k, v in feats.items():
        print(f"  {k}: {v}")

