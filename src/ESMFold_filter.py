#!/usr/bin/env python3
"""
Folding proteins with ESMFold to check foldability of generated sequences.
"""

import argparse
import torch
import esm
import os
import numpy as np
from torch.utils._pytree import tree_map

def main(args):
    model = esm.pretrained.esmfold_v1()
    model = model.eval().cuda()

    os.makedirs(args.structure_output_path, exist_ok=True)

    if os.path.exists(args.filtered_sequences_file):
        os.remove(args.filtered_sequences_file)

    with open(args.sequences_file, 'r') as f:
        sequences = {}
        current_seq_id = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq_id is not None:
                    sequences[current_seq_id] = ''.join(current_seq)
                current_seq_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        

    for seq_id, sequence in sequences.items():
        with torch.no_grad():
        # output = model.infer_pdb(sequence)
            output = model.infer(sequence)
        
        pdb_str = model.output_to_pdb(output)[0]

        output = tree_map(lambda x: x.cpu().numpy(), output)
        ptm = output["ptm"][0]
        plddt = output["plddt"][0,...,1].mean()
        pae = (output["aligned_confidence_probs"][0] * np.arange(64)).mean(-1) * 31
        mask = output["atom37_atom_exists"][0,:,1] == 1
        pae = pae[mask,:][:,mask]
        pae = pae.mean()

        pass_flag = False

        if (plddt.mean().item() >= args.plddt_threshold and
            ptm.item() >= args.ptm_threshold and
            pae.mean().item() <= args.pae_threshold):
            pass_flag = True

            pdb_path = os.path.join(args.structure_output_path, f"{seq_id}.pdb")
            with open(pdb_path, 'w') as pdb_file:
                pdb_file.write(pdb_str)

            with open(args.filtered_sequences_file, 'a') as f:
                f.write(f">{seq_id}\n")
                f.write(f"{sequence}\n")

        pass_flag = "✅" if pass_flag else "❌"

        print(f"Sequence ID: {seq_id}, pLDDT: {plddt.mean().item():.2f}, pTM: {ptm.item():.2f}, PAE: {pae.mean().item():.2f}, Foldable: {pass_flag}")

    

def get_args():
    parser = argparse.ArgumentParser(
        description="Fold generated sequences with ESMFold to check foldability."
    )
    parser.add_argument("--sequences_file", type=str, required=True, help="Path to the file containing generated sequences in FASTA format.")
    parser.add_argument("--structure_output_path", type=str, default="./output/esmfold", help="Directory to save the folding results.")
    parser.add_argument("--filtered_sequences_file", type=str, default="./output/esmfold_filtered_sequences.fasta", help="Path to save the filtered foldable sequences.")

    
    parser.add_argument("--plddt_threshold", type=float, default=80.0, help="Minimum average pLDDT score to consider a sequence as foldable.")
    parser.add_argument("--pae_threshold", type=float, default=7.0, help="Maximum average PAE score to consider a sequence as foldable.")
    parser.add_argument("--ptm_threshold", type=float, default=0.8, help="Minimum pTM score to consider a sequence as foldable.")

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = get_args()

    main(args)