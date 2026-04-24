"""
pdb_to_xyz.py

Converts downloaded DNA PDB files from w3DNA to XYZ format.

Usage:
    python pdb_to_xyz.py
"""

from pathlib import Path


def pdb_to_xyz(pdb_file, xyz_file):
    """
    Convert PDB file to XYZ format.
    """
    atoms = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                element = ''
                
                if len(line) >= 78:
                    element = line[76:78].strip()
                
                if not element:
                    atom_name = line[12:16].strip()
                    if atom_name[0].isdigit():
                        element = atom_name[1] if len(atom_name) > 1 else 'H'
                    else:
                        element = atom_name[0]
                
                atoms.append((element, x, y, z))
                
            except (ValueError, IndexError):
                continue
    
    with open(xyz_file, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"Converted from {pdb_file.name}\n")
        
        for elem, x, y, z in atoms:
            f.write(f"{elem:2s}  {x:12.6f}  {y:12.6f}  {z:12.6f}\n")
    
    return len(atoms)


def main():
    print("=" * 60)
    print("PDB TO XYZ CONVERTER")
    print("=" * 60)
    
    # Input and output in same folder
    dna_dir = Path(r"C:\Users\danie\OneDrive - Johns Hopkins\Desktop\Research\5. AI and ML HEHEHEHEH\bitchass toy dataset from 2005\Anderson-ToyDataSet\XYZ\DNA_pbds_from_web3x")
    
    print(f"\nDirectory: {dna_dir}")
    
    # Find PDB files
    pdb_files = list(dna_dir.glob("*.pdb"))
    
    if not pdb_files:
        print(f"\n[ERROR] No PDB files found!")
        return
    
    print(f"\nFound {len(pdb_files)} PDB file(s):")
    for p in pdb_files:
        print(f"  - {p.name}")
    
    # Convert each PDB to XYZ
    print(f"\nConverting...")
    print("-" * 60)
    
    for pdb_file in sorted(pdb_files):
        xyz_file = pdb_file.with_suffix('.xyz')
        
        try:
            atom_count = pdb_to_xyz(pdb_file, xyz_file)
            print(f"  [OK] {pdb_file.name:30s} -> {xyz_file.name} ({atom_count} atoms)")
        except Exception as e:
            print(f"  [ERROR] {pdb_file.name}: {e}")
    
    # Summary
    print("\n" + "=" * 60)
    print("OUTPUT FILES")
    print("=" * 60)
    print(f"\nFiles in: {dna_dir}\n")
    
    for f in sorted(dna_dir.glob("*")):
        size = f.stat().st_size
        size_str = f"{size/1024:.1f} KB" if size > 1024 else f"{size} B"
        print(f"  {f.name:30s} {size_str:>10}")
    
    print("\n" + "=" * 60)
    print("DONE!")
    print("=" * 60)


if __name__ == "__main__":
    main()