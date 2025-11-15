#!/usr/bin/env python3
"""
ESM3 Protein Sequence Generator.
"""

import os
import argparse
import sys
import torch
from esm.sdk.api import ESMProtein, GenerationConfig
from esm.utils.structure.protein_chain import ProteinChain
from esm.models.esm3 import ESM3
import time

def main(args):
    model = ESM3.from_pretrained(
        args.model,   
        device=torch.device(args.device)
    )

    template_gst = ESMProtein.from_protein_chain(
        ProteinChain.from_rcsb(args.pdb_id, chain_id=args.chain_id)
    )
    template_gst_tokens = model.encode(template_gst)

    prompt_sequence = ["_"] * len(template_gst.sequence)
    prompt_sequence[5] = "Y"
    prompt_sequence[6] = "W"
    prompt_sequence[43] = "K"
    prompt_sequence[52] = "N"
    prompt_sequence[53] = "L"
    prompt_sequence[65] = "Q"
    prompt_sequence[66] = "S"
    prompt_sequence[99] = "D"
    prompt_sequence[109] = "Y"
    prompt_sequence = "".join(prompt_sequence)

    print(template_gst.sequence)
    print(prompt_sequence)

    prompt = model.encode(ESMProtein(sequence=prompt_sequence))

    prompt.structure = torch.full_like(prompt.sequence, 4096)
    prompt.structure[0] = 4098
    prompt.structure[-1] = 4097

    prompt.structure[5] = template_gst_tokens.structure[5]
    prompt.structure[6] = template_gst_tokens.structure[6]
    prompt.structure[43] = template_gst_tokens.structure[43]
    prompt.structure[52] = template_gst_tokens.structure[52]
    prompt.structure[53] = template_gst_tokens.structure[53]
    prompt.structure[65] = template_gst_tokens.structure[65]
    prompt.structure[66] = template_gst_tokens.structure[66]
    prompt.structure[99] = template_gst_tokens.structure[99]
    prompt.structure[109] = template_gst_tokens.structure[109]

    print("".join(["✔" if st < 4096 else "_" for st in prompt.structure]))

    num_tokens_to_decode_structure = min((prompt.structure == 4096).sum().item(), 20)
    num_tokens_to_decode_sequence = min((prompt.sequence == 32).sum().item(), 20)

    num_sequence_generated = 0
    sequences = []

    template_chain = template_gst.to_protein_chain()
    constrained_site_positions = [5, 6, 43, 52, 53, 65, 66, 99, 109]

    config_structure = GenerationConfig(
        track="structure",
        num_steps=num_tokens_to_decode_structure,
        temperature=args.temperature,
    )
    config_sequence = GenerationConfig(
        track="sequence", 
        num_steps=num_tokens_to_decode_sequence, 
        temperature=args.temperature,
    )
    config_refold = GenerationConfig(
        track="structure", 
        num_steps=1, 
        temperature=0.0
    )

    while num_sequence_generated < args.num_sequences:
        prompts_list = [prompt] * args.batch_size
        configs_structure_list = [config_structure] * args.batch_size
        configs_sequence_list = [config_sequence] * args.batch_size
        configs_refold_list = [config_refold] * args.batch_size

        with torch.no_grad():
            structure_tensors = model.batch_generate(
                prompts_list,
                configs_structure_list
            )

            sequence_tensors = model.batch_generate(
                structure_tensors,
                configs_sequence_list
            )

            del structure_tensors

            for tensor in sequence_tensors:
                tensor.structure = None

            refolded_tensors = model.batch_generate(
                sequence_tensors,
                configs_refold_list
            )

            del sequence_tensors
        

        for i, final_tensor in enumerate(refolded_tensors):
            if num_sequence_generated >= args.num_sequences:
                break 

            protein = model.decode(final_tensor)
            generated_sequence = protein.sequence
            generation_chain = protein.to_protein_chain()

            constrained_site_rmsd = template_chain[constrained_site_positions].rmsd(
                generation_chain[constrained_site_positions]
            )
            if constrained_site_rmsd >= args.max_constraint_rmsd:
                del protein, generation_chain, final_tensor 
                continue

            backbone_rmsd = template_chain.rmsd(generation_chain)
            if backbone_rmsd <= args.min_backbone_rmsd:
                del protein, generation_chain, final_tensor 
                continue

            num_sequence_generated += 1
            print(f"Collected {num_sequence_generated}/{args.num_sequences}: Constrained Site RMSD={constrained_site_rmsd:.3f}, Backbone RMSD={backbone_rmsd:.3f}")

            sequences.append(generated_sequence)

            del protein, generation_chain, final_tensor
        
        del refolded_tensors
        
        torch.cuda.empty_cache() 
    
    sequences = sequences[:args.num_sequences]

    os.makedirs(args.output_path, exist_ok=True)

    with open(f"{args.output_path}/{args.output_file}", "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">generated_sequence_{i+1}\n")
            f.write(f"{seq}\n")


def get_args():
    parser = argparse.ArgumentParser(
        description="Generate protein sequences using ESM3 model."
    )
    parser.add_argument("--model", type=str, default="esm3_sm_open_v1", help="Name of the ESM3 model to use.")
    parser.add_argument("--device", type=str, default="cuda", help="Device to run the model on.")

    parser.add_argument("--pdb_id", type=str, default="1Y6E", help="PDB ID of the protein structure.")
    parser.add_argument("--chain_id", type=str, default="A", help="Chain ID of the protein structure.")
    parser.add_argument("--num_sequences", type=int, default=100, help="Number of sequences to generate.")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size for generation.")
    parser.add_argument("--output_path", type=str, default="./output", help="Output directory for generated sequences.")
    
    parser.add_argument("--output_file", type=str, default="generated_sequences.fasta", help="Output FASTA file for generated sequences.")

    parser.add_argument("--temperature", type=float, default=1.0, help="Sampling temperature for generation.")
    parser.add_argument("--max_constraint_rmsd", type=float, default=1.5, help="Maximum constraint RMSD for generation.")
    parser.add_argument("--min_backbone_rmsd", type=float, default=1.5, help="Minimum backbone RMSD for generation.")

    args = parser.parse_args()
    
    if args.output_file == "generated_sequences.fasta":
        args.output_file = f"generated_sequences_{time.time()}.fasta"
        
    return args

if __name__ == "__main__":
    args = get_args()

    main(args)