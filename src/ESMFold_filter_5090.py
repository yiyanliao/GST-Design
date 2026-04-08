#!/usr/bin/env python3
"""
Folding proteins with ESMFold to check foldability of generated sequences.
"""

import argparse
import torch
import os
import numpy as np
import gc

from transformers import EsmForProteinFolding, AutoTokenizer
# 抛弃不稳定的上层 API，直接从底层核心依赖中导入 PDB 构建工具
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

# 自动设置 Hugging Face 镜像源，解决国内无法直连的问题
os.environ["HF_ENDPOINT"] = "https://hf-mirror.com"

def convert_outputs_to_pdb(outputs):
    """
    官方推荐的 fallback 方案，手动将模型的底层张量组装为 PDB 文本格式。
    彻底绕开 transformers 版本差异导致的 output_to_pdb 导入失败问题。
    """
    final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
    outputs_numpy = {k: v.to("cpu").numpy() for k, v in outputs.items()}
    final_atom_positions = final_atom_positions.cpu().numpy()
    final_atom_mask = outputs_numpy["atom37_atom_exists"]
    
    pdbs = []
    for i in range(outputs_numpy["aatype"].shape[0]):
        aa = outputs_numpy["aatype"][i]
        pred_pos = final_atom_positions[i]
        mask = final_atom_mask[i]
        resid = outputs_numpy["residue_index"][i] + 1
        
        # 组装 PDB 需要的各层基础信息
        pred = OFProtein(
            aatype=aa,
            atom_positions=pred_pos,
            atom_mask=mask,
            residue_index=resid,
            b_factors=outputs_numpy["plddt"][i],
            chain_index=outputs_numpy["chain_index"][i] if "chain_index" in outputs_numpy else None,
        )
        pdbs.append(to_pdb(pred))
    return pdbs

def main(args):

    torch.set_num_threads(1)

    print(f"CUDA 是否可用: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"当前显卡: {torch.cuda.get_device_name(0)}")

    # 务必保留 low_cpu_mem_usage=True，这是拯救你 15GB 内存的唯一关键！
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", 
        device_map="cuda:0",
        low_cpu_mem_usage=True
    )
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = model.eval()
    
    if torch.cuda.is_available():
        # 严格对齐 Meta 官方的混合精度策略：语言模型主干强制 FP16
        # 解决 HuggingFace 默认全 FP32 导致的 pLDDT 分布异常偏低问题
        model.esm = model.esm.half()
        torch.backends.cuda.matmul.allow_tf32 = True

    # 限制 5090 显存占用，防止长序列导致 VRAM OOM
    model.trunk.set_chunk_size(64)

    os.makedirs(args.structure_output_path, exist_ok=True)

    if os.path.exists(args.filtered_sequences_file):
        os.remove(args.filtered_sequences_file)

    def read_fasta(path):
        with open(path, 'r') as f:
            seq_id, seq = None, []
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if seq_id: yield seq_id, "".join(seq)
                    seq_id, seq = line[1:], []
                else:
                    seq.append(line)
            if seq_id: yield seq_id, "".join(seq)

    for seq_id, sequence in read_fasta(args.sequences_file):
        if len(sequence) > 1000:
            print(f"Skipping {seq_id}: Sequence too long ({len(sequence)}aa)")
            continue

        with torch.no_grad():
            inputs = tokenizer([sequence], return_tensors='pt', add_special_tokens=False)
            inputs = {k: v.cuda() for k, v in inputs.items()}

            outputs = model(**inputs)
            plddt_val = outputs.plddt.mean().item() * 100.0
            ptm_val = outputs.ptm.item()
            pae_val = outputs.predicted_aligned_error.mean().item()

            pass_flag = False
            
            if (plddt_val >= args.plddt_threshold and
                ptm_val >= args.ptm_threshold and
                pae_val <= args.pae_threshold):
                pass_flag = True

                with open(args.filtered_sequences_file, 'a') as f:
                    f.write(f">{seq_id}\n")
                    f.write(f"{sequence}\n")

                # 显式暴露 PDB 写入错误（不静默吞异常）
                try:
                    pdb_strs = convert_outputs_to_pdb(outputs)
                    pdb_filename = os.path.join(args.structure_output_path, f"{seq_id}.pdb")
                    with open(pdb_filename, "w") as f:
                        f.write(pdb_strs[0])
                except Exception as e:
                    print(f"[{seq_id}] X 保存 PDB 失败，发生异常: {e}")

            status_icon = "✅" if pass_flag else "X"

            print(f"Sequence ID: {seq_id}, pLDDT: {plddt_val:.2f}, pTM: {ptm_val:.2f}, PAE: {pae_val:.2f}, Foldable: {status_icon}")
            
            del outputs, inputs
            gc.collect()
            torch.cuda.empty_cache()

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