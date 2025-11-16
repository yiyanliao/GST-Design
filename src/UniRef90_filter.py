#!/usr/bin/env python3
"""
Filter a FASTA file based on alignment identity from an M8 file,
and copy corresponding PDB files for the sequences that are kept.
"""

import argparse
import os
import shutil  # Import the shutil module for file copying

def main(args):
        
    print(f"Reading M8 file: {args.m8_file}")
    seq_ids_to_remove = set()
    try:
        with open(args.m8_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                
                if len(parts) == 12:
                    try:
                        query_id = parts[0]
                        identity = float(parts[2])
                        
                        # If identity is *greater* than threshold, mark for removal
                        if identity > args.identity_threshold:
                            seq_ids_to_remove.add(query_id)
                            
                    except ValueError:
                        print(f"Warning: Skipping malformed line in M8 file: {line}")
                else:
                     print(f"Warning: Skipping malformed line (not 12 columns) in M8 file: {line}")

    except FileNotFoundError:
        print(f"Error: M8 file not found at {args.m8_file}")
        return
    
    print(f"Found {len(seq_ids_to_remove)} unique sequences to remove based on identity > {args.identity_threshold}")
        
    print(f"Reading source FASTA file: {args.fasta_file}")
    sequences = {}
    try:
        with open(args.fasta_file, 'r') as f:
            current_seq_id = None
            current_seq = []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq_id is not None:
                        sequences[current_seq_id] = ''.join(current_seq)
                    
                    current_seq_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_seq_id is not None:
                sequences[current_seq_id] = ''.join(current_seq)
                
    except FileNotFoundError:
        print(f"Error: Source FASTA file not found at {args.fasta_file}")
        return

    print(f"Loaded {len(sequences)} total sequences from FASTA.")

    if os.path.exists(args.output_fasta_file):
        os.remove(args.output_fasta_file)

    os.makedirs(args.pdb_output_dir, exist_ok=True)
    print(f"Will copy filtered PDBs to: {args.pdb_output_dir}")

    kept_count = 0
    removed_count = 0
    pdb_copied_count = 0
    
    with open(args.output_fasta_file, 'a') as f_out:
        for seq_id, sequence in sequences.items():
            
            if seq_id not in seq_ids_to_remove:
                f_out.write(f">{seq_id}\n")
                f_out.write(f"{sequence}\n")
                kept_count += 1
                
                source_pdb_path = os.path.join(args.pdb_input_dir, f"{seq_id}.pdb")
                dest_pdb_path = os.path.join(args.pdb_output_dir, f"{seq_id}.pdb")
                
                try:
                    shutil.copy(source_pdb_path, dest_pdb_path)
                    pdb_copied_count += 1
                except FileNotFoundError:
                    print(f"Warning: PDB file not found at {source_pdb_path}. Cannot copy.")
                except Exception as e:
                    print(f"Warning: Failed to copy {source_pdb_path} due to: {e}")
                    
            else:
                removed_count += 1
    
    print(f"\n--- Filtering Complete ---")
    print(f"Total sequences processed: {len(sequences)}")
    print(f"Sequences kept: {kept_count}")
    print(f"Sequences removed: {removed_count}")
    print(f"Filtered sequences saved to: {args.output_fasta_file}")
    print(f"Corresponding PDBs copied: {pdb_copied_count}")


def get_args():
    parser = argparse.ArgumentParser(
        description="Filter FASTA and copy PDBs based on sequence identity hits from an M8 file."
    )
    
    parser.add_argument("--fasta_file", type=str, required=True, help="Path to the input FASTA file to be filtered.")
    parser.add_argument("--m8_file", type=str, required=True, help="Path to the MMseqs2/BLAST alignment file in M8 format.")
    parser.add_argument("--pdb_input_dir", type=str, default="./output/esmfold", help="Directory containing the original PDB files.")
    
    parser.add_argument("--output_fasta_file", type=str, default="./output/identity_filtered_sequences.fasta", help="Path to save the new, filtered FASTA file.")
    parser.add_argument("--pdb_output_dir", type=str, default="./output/identity_filtered_esmfold", help="Directory to save the PDBs for filtered sequences.")

    parser.add_argument("--identity_threshold", type=float, default=0.3, help="Identity threshold. Sequences with *any* hit > this value will be removed.")

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = get_args()
    main(args)