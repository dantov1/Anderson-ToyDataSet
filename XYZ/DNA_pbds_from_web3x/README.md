# DNA Structure Generation for PBAE-DNA Complex Dataset

## Overview

This document describes the process of generating DNA structures in XYZ format for use in the Anderson-ToyDataSet project. These DNA structures will be combined with PBAE (poly(beta-amino ester)) polymers to create PBAE-DNA complexes for training CGCNN models.

---

## Table of Contents

1. Final Solution: w3DNA
2. Output Files
3. Conversion Script
4. Failed Approaches (For Reference)
5. DNA Structure Background
6. Future Usage

---

## Final Solution: w3DNA

### What is w3DNA?

w3DNA (wDSSR) is a web interface to X3DNA-DSSR for dissecting and modeling 3D nucleic acid structures. It provides pre-computed fiber models based on crystallographic data.

Website: http://web.x3dna.org/

### Step-by-Step Process

#### 1. Navigate to w3DNA
- Go to http://web.x3dna.org/
- Click on "MODEL" in the top navigation menu

#### 2. Select Fiber Model
- You will see a table of 56 available fiber models
- Select Row 4: "B-DNA (calf thymus; generic sequence: A, C, G and T)"
  - Twist: 36.0 degrees
  - Rise: 3.375 Angstroms
- Click the "Action" button for that row

#### 3. Enter Sequence
Paste your desired sequence. We generated structures for:

20 bp:
ATCGATCGATCGATCGATCG

100 bp:
TCCCGAGTGGGTATACAATTAGGCAACGGGATTTAAAGGGCGCATAGCGTGCTAGTCAGTACGTACGATCGATCGTAGCTAGCTAGCTGATCGATCGTAG

150 bp:
CCTTTGATCCCAGTCTGAGTGGTTGGTAGTGTACATCAAAGCTAGCTAGCTGATCGATCGTAGCTAGCTGCATGCATGCTAGCTAGCTGATCGATCGTAGCTAGCTGCATGCATGCTAGCTAGCTGATCGATCGTAGCTAGCTGCAT

200 bp:
CTACATATTGTCGAGTACAAGCTATCAATGGCAGTGTCCTGCTAGCTAGCTGATCGATCGTAGCTAGCTGCATGCATGCTAGCTAGCTGATCGATCGTAGCTAGCTGCATGCATGCTAGCTAGCTGATCGATCGTAGCTAGCTGCATGCATGCTAGCTAGCTGATCGATCGTAGCTAGCTGCATGCATGCTAG

#### 4. Generate and Download
- Click "Generate" or "Submit"
- Download the resulting PDB file
- Rename appropriately (e.g., dna_100bp.pdb)

#### 5. Save Location
Save downloaded PDB files to:
Anderson-ToyDataSet/XYZ/DNA_pbds_from_web3x/

---

## Output Files

### Directory Structure

Anderson-ToyDataSet/
└── XYZ/
    └── DNA_pbds_from_web3x/
        ├── dna_20bp.pdb
        ├── dna_20bp.xyz
        ├── dna_100bp.pdb
        ├── dna_100bp.xyz
        ├── dna_150bp.pdb
        ├── dna_150bp.xyz
        ├── dna_200bp.pdb
        └── dna_200bp.xyz

### File Sizes (Approximate)

| Structure | PDB Size | XYZ Size | Atoms |
|-----------|----------|----------|-------|
| dna_20bp  | ~65 KB   | ~37 KB   | ~800  |
| dna_100bp | 324 KB   | 184 KB   | ~4,000|
| dna_150bp | 477 KB   | 271 KB   | ~6,000|
| dna_200bp | 626 KB   | 356 KB   | ~8,000|

### Atom Count Formula
- Approximately 40 atoms per base pair (both strands)
- 100 bp = ~4,000 atoms
- 150 bp = ~6,000 atoms
- 200 bp = ~8,000 atoms

---

## Conversion Script

### pdb_to_xyz.py

This script converts PDB files downloaded from w3DNA to XYZ format.

Location: scripts/pdb_to_xyz.py

Usage:
    cd scripts
    python pdb_to_xyz.py

The script:
1. Finds all .pdb files in DNA_pbds_from_web3x folder
2. Extracts atom coordinates and element symbols
3. Writes .xyz files in same folder

### PDB Parsing Details

The script reads PDB ATOM records with this column layout:
- Columns 31-38: X coordinate
- Columns 39-46: Y coordinate
- Columns 47-54: Z coordinate
- Columns 77-78: Element symbol

If element is not in columns 77-78, it derives from atom name (column 13-16).

---

## Failed Approaches (For Reference)

### Attempt 1: Simple Tiling

Approach: Download a short DNA structure (PDB 4C64, 12 bp) and tile/repeat it to create longer structures.

Problem: The tiles were not covalently connected. At each junction, the O3'-P bond distance was much larger than the ideal 1.6 Angstroms (typically 5-30 Angstrom gaps).

Visual Result: In molecular viewers, the structure appeared as separate 12 bp chunks floating near each other, not a continuous helix.

### Attempt 2: Build Nucleotide-by-Nucleotide (Single Chain Templates)

Approach: Extract nucleotide templates from PDB 4C64, then build DNA one nucleotide at a time, positioning each to maintain proper O3'-P bonds (~1.6 Angstroms).

Problem: Used the same templates for both sense and antisense strands, applying rotations to flip the antisense strand. However, the backbone geometry is different for each strand direction.

Result: 
- Chain A (sense): All O3'-P bonds correct (~1.6 Angstroms) - SUCCESS
- Chain B (antisense): O3'-P distances 5-40 Angstroms - FAILED

### Attempt 3: Chain-Specific Templates

Approach: Extract separate templates from Chain A (sense) and Chain B (antisense) of PDB 4C64, then use the appropriate template for each strand.

Problem: While backbone connectivity was perfect (all bonds 1.6 Angstroms), the two strands were not intertwined - they appeared as two separate parallel helices instead of the characteristic double helix structure.

Result:
- Backbone bonds: All correct - SUCCESS
- Double helix structure: Not formed - FAILED

### Why w3DNA Works

w3DNA uses pre-computed fiber models based on experimental crystallographic data. These models include:

1. Correct helical parameters (twist, rise, inclination, etc.)
2. Proper strand-strand relationships (base pairing geometry)
3. Accurate backbone torsion angles (alpha, beta, gamma, delta, epsilon, zeta)
4. Validated atomic coordinates from decades of structural biology research

Building these from scratch requires solving a complex geometric optimization problem that goes far beyond simple template positioning.

---

## DNA Structure Background

### B-DNA Parameters

| Parameter            | Value        |
|----------------------|--------------|
| Helix Type           | Right-handed |
| Base pairs per turn  | 10.5         |
| Twist per bp         | 36.0 degrees |
| Rise per bp          | 3.375 A      |
| Helix diameter       | ~20 A        |
| Major groove width   | ~22 A        |
| Minor groove width   | ~12 A        |

### XYZ Format

<number of atoms>
<comment line>
<element>  <x>  <y>  <z>
<element>  <x>  <y>  <z>
...

### DNA Nucleotide Atoms (Approximate)

| Nucleotide | Atoms |
|------------|-------|
| dA (deoxyadenosine)  | 21 |
| dT (deoxythymidine)  | 20 |
| dG (deoxyguanosine)  | 22 |
| dC (deoxycytidine)   | 16 |

Average: ~20 atoms per nucleotide x 2 strands = ~40 atoms per base pair

---

## Future Usage

### Adding New DNA Lengths

1. Go to http://web.x3dna.org/
2. Click MODEL then Row 4 (B-DNA)
3. Enter sequence of desired length
4. Download PDB
5. Save to Anderson-ToyDataSet/XYZ/DNA_pbds_from_web3x/
6. Run: python pdb_to_xyz.py

### Sequence Generation

For random sequences with ~50% GC content, use Python:

    import random
    
    def generate_sequence(length, gc_content=0.5, seed=42):
        random.seed(seed)
        seq = []
        for _ in range(length):
            if random.random() < gc_content:
                seq.append(random.choice(['G', 'C']))
            else:
                seq.append(random.choice(['A', 'T']))
        return ''.join(seq)
    
    print(generate_sequence(50))   # 50 bp
    print(generate_sequence(100))  # 100 bp

### Verification

Always verify generated structures in a molecular viewer:
- Mol* (Molstar): https://molstar.org/viewer/
- NGL Viewer: https://nglviewer.org/ngl/

Check for:
- Double helix shape
- Two intertwined strands
- Continuous backbone (no gaps)
- Proper base pairing

### Alternative Web Tools

If w3DNA is unavailable:

1. make-na: http://structure.usc.edu/make-na/server.html
2. 3D-DART: https://haddock.science.uu.nl/3ddart/

---

## References

- w3DNA/DSSR: http://web.x3dna.org/
- X3DNA: https://x3dna.org/
- PDB Format: https://www.wwpdb.org/documentation/file-format
- B-DNA Structure: Dickerson, R.E. et al. (1982) "The anatomy of A-, B-, and Z-DNA" Science

---

## Authors & Date

- Generated: 2024
- Project: Anderson-ToyDataSet (PBAE-DNA complexes for CGCNN)
- Institution: Johns Hopkins University

---

## Appendix: Reference PDB Structures

If you need to examine real DNA crystal structures:

| PDB ID | Description              | Base Pairs |
|--------|--------------------------|------------|
| 1BNA   | Dickerson-Drew dodecamer | 12         |
| 4C64   | B-DNA dodecamer          | 12         |
| 1D49   | B-DNA decamer            | 10         |

Download from: https://www.rcsb.org/