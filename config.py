import numpy as np

# ==============================================================================
# GLOBAL VARIABLES for VOCABULARY
# ==============================================================================


# Special vocabulary symbols - we always put them at the start.
_PAD = "_PAD"
_GO = "_GO"
_EOS = "_EOS"
_START_VOCAB = [_PAD, _GO, _EOS]

PAD_ID = 0
GO_ID = 1
EOS_ID = 2
assert PAD_ID == 0




vocab_reverse = ['A',
                 'R',
                 'N',
#                  'N(Deamidation)',
                 'D',
#                  'C',
                 'C(Carbamidomethylation)',
                 'E',
                 'Q',
#                  'Q(Deamidation)',
                 'G',
                 'H',
                 'I',
                 'L',
                 'K',
                 'M',
                 'M(Oxidation)',
                 'F',
                 'P',
                 'S',
#                  'S(Phosphorylation)',
                 'T',
#                  'T(Phosphorylation)',
                 'W',
                 'Y',
#                  'Y(Phosphorylation)',
                 'V',
                 ]

vocab_reverse = _START_VOCAB + vocab_reverse
print("Training vocab_reverse ", vocab_reverse)

vocab = dict([(x, y) for (y, x) in enumerate(vocab_reverse)])
print("Training vocab ", vocab)

vocab_size = len(vocab_reverse)
print("Training vocab_size ", vocab_size)

# mass value
mass_H = 1.0078
mass_H2O = 18.0106
mass_NH3 = 17.0265
mass_N_terminus = 1.0078
mass_C_terminus = 17.0027
mass_CO = 27.9949
mass_Phosphorylation = 79.96633

# mass_AA should be comprehensive, including the mass for all common ptm
mass_AA = {'_PAD': 0.0,
           '_GO': mass_N_terminus - mass_H,
           '_EOS': mass_C_terminus + mass_H,
           'A': 71.03711,  # 0
           'R': 156.10111,  # 1
           'N': 114.04293,  # 2
           'N(Deamidation)': 115.02695,
           'D': 115.02694,  # 3
           'C': 103.00919,  # 4
           'C(Carbamidomethylation)': 160.03065,  # C(+57.02)
           # ~ 'C(Carbamidomethylation)': 161.01919, # C(+58.01) # orbi
           'E': 129.04259,  # 5
           'Q': 128.05858,  # 6
           'Q(Deamidation)': 129.0426,
           'G': 57.02146,  # 7
           'H': 137.05891,  # 8
           'I': 113.08406,  # 9
           'L': 113.08406,  # 10
           'K': 128.09496,  # 11
           'M': 131.04049,  # 12
           'M(Oxidation)': 147.0354,
           'F': 147.06841,  # 13
           'P': 97.05276,  # 14
           'S': 87.03203,  # 15
           'S(Phosphorylation)': 87.03203 + mass_Phosphorylation,
           'T': 101.04768,  # 16
           'T(Phosphorylation)': 101.04768 + mass_Phosphorylation,
           'W': 186.07931,  # 17
           'Y': 163.06333,  # 18
           'Y(Phosphorylation)': 163.06333 + mass_Phosphorylation,
           'V': 99.06841,  # 19
           }

mass_ID = [mass_AA[vocab_reverse[x]] for x in range(vocab_size)]
mass_ID_np = np.array(mass_ID, dtype=np.float32)

mass_AA_min = mass_AA["G"]  # 57.02146

MZ_MAX = 8000.0

