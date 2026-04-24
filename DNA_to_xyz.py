"""
dna_from_pdb.py

Downloads PDB 4C64 and converts to XYZ format.
Tiles the structure to create 100, 150, 200 bp DNA.

Usage:
    python dna_from_pdb.py
"""

import urllib.request
import numpy as np
from pathlib import Path
from collections import defaultdict

# =============================================================================
# CONFIGURATION
# =============================================================================

PDB_ID = "4C64"
OUTPUT_SUBDIR = Path("Anderson-ToyDataSet") / "XYZ" / "DNA"
TARGET_LENGTHS = [100, 150, 200]
RISE_PER_BP = 3.38
TWIST_PER_BP = 36.0


# =============================================================================
# PATH UTILITIES
# =============================================================================

def find_project_root():
    """Find project root by looking for Anderson-ToyDataSet folder."""
    current = Path.cwd()
    
    for path in [current, current.parent, current.parent.parent]:
        if (path / "Anderson-ToyDataSet").exists():
            return path
    
    print(f"WARNING: Could not find Anderson-ToyDataSet folder.")
    return current


# =============================================================================
# DOWNLOAD PDB
# =============================================================================

def download_pdb(pdb_id, output_dir):
    """Download PDB file from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = Path(output_dir) / f"{pdb_id}.pdb"
    
    if output_path.exists():
        print(f"  [OK] {pdb_id}.pdb already exists")
        return output_path
    
    print(f"  Downloading {pdb_id} from RCSB...")
    urllib.request.urlretrieve(url, output_path)
    print(f"  [OK] Saved to {output_path.name}")
    
    return output_path


# =============================================================================
# PARSE PDB
# =============================================================================

def parse_pdb(pdb_file):
    """Parse PDB file and extract atoms."""
    atoms = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            
            try:
                atom = {
                    'serial': int(line[6:11]),
                    'name': line[12:16].strip(),
                    'resname': line[17:20].strip(),
                    'chain': line[21].strip(),
                    'resnum': int(line[22:26]),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                }
                
                if len(line) >= 78:
                    atom['element'] = line[76:78].strip()
                
                if not atom.get('element'):
                    name = atom['name']
                    if name[0].isdigit():
                        atom['element'] = name[1] if len(name) > 1 else 'X'
                    else:
                        atom['element'] = name[0]
                
                atoms.append(atom)
                
            except (ValueError, IndexError):
                continue
    
    return atoms


def is_dna_residue(resname):
    """Check if residue is DNA."""
    dna_names = {
        'DA', 'DT', 'DG', 'DC', 'DU', 'DI',
        'DA3', 'DA5', 'DT3', 'DT5', 'DG3', 'DG5', 'DC3', 'DC5',
        'A', 'T', 'G', 'C', 'U',
    }
    return resname.upper() in dna_names


def filter_dna_atoms(atoms):
    """Keep only DNA atoms."""
    return [a for a in atoms if is_dna_residue(a['resname'])]


def analyze_structure(atoms):
    """Analyze and print structure info."""
    residues = defaultdict(list)
    for atom in atoms:
        key = (atom['chain'], atom['resnum'], atom['resname'])
        residues[key].append(atom)
    
    chains = defaultdict(list)
    for (chain, resnum, resname) in residues.keys():
        chains[chain].append((resnum, resname))
    
    for chain in chains:
        chains[chain].sort(key=lambda x: x[0])
    
    print(f"\n  Structure analysis:")
    print(f"    Total atoms: {len(atoms)}")
    print(f"    Chains: {list(chains.keys())}")
    
    total_bp = 0
    for chain, res_list in sorted(chains.items()):
        seq = ''.join(r[1][-1] if len(r[1]) <= 2 else r[1][1] for r in res_list)
        print(f"    Chain {chain}: {len(res_list)} residues")
        print(f"      Sequence: {seq[:50]}{'...' if len(seq) > 50 else ''}")
        total_bp += len(res_list)
    
    bp_count = total_bp // 2 if len(chains) >= 2 else total_bp
    print(f"    Estimated base pairs: {bp_count}")
    
    return {
        'atoms': len(atoms),
        'chains': dict(chains),
        'residues': len(residues),
        'bp_estimate': bp_count,
    }


# =============================================================================
# COORDINATE MANIPULATION
# =============================================================================

def get_coordinates(atoms):
    """Extract coordinates as numpy array."""
    return np.array([[a['x'], a['y'], a['z']] for a in atoms])


def center_structure(atoms):
    """Center structure at origin."""
    coords = get_coordinates(atoms)
    center = coords.mean(axis=0)
    
    for atom in atoms:
        atom['x'] -= center[0]
        atom['y'] -= center[1]
        atom['z'] -= center[2]
    
    return atoms


def get_helix_axis(atoms):
    """Estimate helix axis using PCA."""
    coords = get_coordinates(atoms)
    centered = coords - coords.mean(axis=0)
    cov = np.cov(centered.T)
    eigenvalues, eigenvectors = np.linalg.eig(cov)
    axis_idx = np.argmax(eigenvalues)
    axis = eigenvectors[:, axis_idx].real
    return axis, eigenvalues[axis_idx].real


def rotation_matrix_axis_angle(axis, angle_deg):
    """Create rotation matrix for rotation around arbitrary axis."""
    angle = np.radians(angle_deg)
    axis = axis / np.linalg.norm(axis)
    
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    x, y, z = axis
    
    return np.array([
        [t*x*x + c,    t*x*y - s*z,  t*x*z + s*y],
        [t*x*y + s*z,  t*y*y + c,    t*y*z - s*x],
        [t*x*z - s*y,  t*y*z + s*x,  t*z*z + c]
    ])


def tile_dna(atoms, original_bp, target_bp):
    """
    Tile DNA structure to reach target length.
    Updates residue numbers for each copy so viewers display correctly.
    """
    if target_bp <= original_bp:
        return [a.copy() for a in atoms]
    
    n_copies = int(np.ceil(target_bp / original_bp))
    
    coords = get_coordinates(atoms)
    axis, _ = get_helix_axis(atoms)
    
    if axis[2] < 0:
        axis = -axis
    
    translation = axis * (original_bp * RISE_PER_BP)
    rotation_angle = original_bp * TWIST_PER_BP
    
    # Find max residue number per chain
    max_resnum = {}
    for atom in atoms:
        chain = atom['chain']
        if chain not in max_resnum:
            max_resnum[chain] = 0
        max_resnum[chain] = max(max_resnum[chain], atom['resnum'])
    
    tiled_atoms = []
    
    for i in range(n_copies):
        rot = rotation_matrix_axis_angle(axis, i * rotation_angle)
        
        for atom in atoms:
            new_atom = atom.copy()
            
            pos = np.array([atom['x'], atom['y'], atom['z']])
            pos = rot @ pos
            pos = pos + i * translation
            
            new_atom['x'] = pos[0]
            new_atom['y'] = pos[1]
            new_atom['z'] = pos[2]
            new_atom['serial'] = len(tiled_atoms) + 1
            
            # FIX: Update residue number for each tile
            chain = atom['chain']
            new_atom['resnum'] = atom['resnum'] + (i * max_resnum.get(chain, 12))
            
            tiled_atoms.append(new_atom)
    
    return tiled_atoms


# =============================================================================
# OUTPUT FILES
# =============================================================================

def write_xyz(atoms, filename, title="DNA Structure"):
    """Write XYZ file."""
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{title}\n")
        
        for atom in atoms:
            elem = atom['element']
            x, y, z = atom['x'], atom['y'], atom['z']
            f.write(f"{elem:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
    
    return filename


def write_pdb(atoms, filename, title="DNA Structure"):
    """Write PDB file with proper residue numbering."""
    with open(filename, 'w') as f:
        f.write(f"TITLE     {title}\n")
        f.write(f"REMARK    Generated from {PDB_ID}\n")
        
        for i, atom in enumerate(atoms, 1):
            name = atom['name']
            resname = atom['resname']
            chain = atom['chain'] if atom['chain'] else 'A'
            resnum = atom['resnum'] % 10000  # PDB limit
            x, y, z = atom['x'], atom['y'], atom['z']
            elem = atom['element']
            
            f.write(
                f"ATOM  {i:5d} {name:<4s} {resname:3s} {chain:1s}{resnum:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}\n"
            )
        
        f.write("END\n")
    
    return filename


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print(f"DNA STRUCTURE GENERATOR FROM PDB {PDB_ID}")
    print("=" * 70)
    
    project_root = find_project_root()
    output_dir = project_root / OUTPUT_SUBDIR
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\nProject root: {project_root}")
    print(f"Output directory: {output_dir}")
    
    # Download
    print(f"\n[1/4] Downloading {PDB_ID}...")
    pdb_file = download_pdb(PDB_ID, output_dir)
    
    # Parse
    print(f"\n[2/4] Parsing structure...")
    all_atoms = parse_pdb(pdb_file)
    print(f"  Total atoms in file: {len(all_atoms)}")
    
    dna_atoms = filter_dna_atoms(all_atoms)
    print(f"  DNA atoms: {len(dna_atoms)}")
    
    if len(dna_atoms) == 0:
        print("\n[ERROR] No DNA atoms found!")
        return
    
    info = analyze_structure(dna_atoms)
    original_bp = info['bp_estimate']
    
    # Center
    print(f"\n[3/4] Processing structure...")
    dna_atoms = center_structure(dna_atoms)
    print(f"  Centered at origin")
    
    # Generate files
    print(f"\n[4/4] Generating output files...")
    
    # Original
    orig_xyz = output_dir / f"{PDB_ID}_dna.xyz"
    orig_pdb = output_dir / f"{PDB_ID}_dna.pdb"
    write_xyz(dna_atoms, orig_xyz, f"{PDB_ID} DNA - {original_bp} bp")
    write_pdb(dna_atoms, orig_pdb, f"{PDB_ID} DNA - {original_bp} bp")
    print(f"  [OK] {orig_xyz.name} ({len(dna_atoms)} atoms)")
    print(f"  [OK] {orig_pdb.name}")
    
    # Tiled structures
    for target_bp in TARGET_LENGTHS:
        print(f"\n  Creating {target_bp} bp structure...")
        
        if target_bp <= original_bp:
            tiled = [a.copy() for a in dna_atoms]
            print(f"    Using original ({original_bp} bp >= {target_bp} bp)")
        else:
            tiled = tile_dna(dna_atoms, original_bp, target_bp)
            n_copies = int(np.ceil(target_bp / original_bp))
            print(f"    Tiled {n_copies}x ({original_bp} bp x {n_copies} = ~{original_bp * n_copies} bp)")
        
        xyz_file = output_dir / f"dna_{target_bp}bp.xyz"
        pdb_file = output_dir / f"dna_{target_bp}bp.pdb"
        
        write_xyz(tiled, xyz_file, f"B-DNA {target_bp} bp (from {PDB_ID})")
        write_pdb(tiled, pdb_file, f"B-DNA {target_bp} bp (from {PDB_ID})")
        
        print(f"    [OK] {xyz_file.name} ({len(tiled)} atoms)")
        print(f"    [OK] {pdb_file.name}")
    
    # Summary
    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)
    print(f"\nFiles in: {output_dir}\n")
    
    for f in sorted(output_dir.glob("*")):
        size = f.stat().st_size
        size_str = f"{size/1024:.1f} KB" if size > 1024 else f"{size} B"
        print(f"  {f.name:<30} {size_str:>10}")


if __name__ == "__main__":
    main()