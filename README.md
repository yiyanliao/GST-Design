# *De novo* design of GST using ESM3 - Biochemistry Lab

![Workflow](https://github.com/yiyanliao/GST-Design/blob/main/images/workflow.png)

## Installations

- Env for [ESM3](https://github.com/evolutionaryscale/esm):

```bash
conda create -n esm python=3.14
conda activate esm
pip install esm
```

- Env for [MMseqs2](https://github.com/soedinglab/MMseqs2):

```bash
conda create -n mmseqs2 -c conda-forge -c bioconda python=3.9 mmseqs2
conda activate mmseqs2
```

- Env for [ESMFold](https://github.com/facebookresearch/esm), please use `esmfold_environment.yaml` file:

```bash
conda env create -f esmfold_environment.yaml
conda activate esmfold
esm-fold -h
```

- Env for [AutoDock-Vina](https://autodock-vina.readthedocs.io/en/latest/installation.html):

```bash
conda create -n vina python=3.10
conda activate vina
conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme
pip install -U numpy scipy rdkit vina meeko gemmi prody
conda install bioconda::autodock-vina
```

## Downloads

- Due to network issues in China, we have to download checkpoints for ESM3 beforehand and "hack" ESM3 source code (no need if you have better internet on your server, refer to [ESM3 repo](https://github.com/evolutionaryscale/esm)):

ESM3 model download (you may try this several (or many) times):

```bash
hf auth login
hf download EvolutionaryScale/esm3-sm-open-v1 --local-dir esm3-sm-open-v1
```

"Hack" ESM3 source code:

change `/path/to/anaconda3/envs/esm/lib/python3.14/site-packages/esm/utils/constants/esm3.py` line `105`, from:

```python
path = Path(snapshot_download(repo_id="EvolutionaryScale/esm3-sm-open-v1"))
```

To:

```python
path = Path("/path/to/esm3-sm-open-v1") 
```

- Then you need to download `UniRef90` database (~ 40 GB):

```bash
conda activate mmseqs2
mmseqs databases UniRef90 uniref90 tmp
```

## Usages

- Generate 1000 potential candidates compiled with structural metrics:

```bash
conda activate esm
python src/ESM3_generate.py --num_sequences 1000 --batch_size 8
```

- Filter the sequences based on foldability with ESMFold:

```
conda activate esmfold
python src/ESMFold_filter.py --sequences_file ./output/generated_sequences_<time>.fasta
```

- Calculate the identity of candidates towards `UniRef90` and collect the ones `<50%` or have no hit:

```bash
conda activate mmseqs2
mmseqs easy-search --alignment-mode 3 -s 7 ./output/esmfold_filtered_sequences.fasta uniref90 ./output/esmfold_filtered_sequences.m8 tmp 
python src/UniRef90_filter.py --fasta_file ./output/esmfold_filtered_sequences.fasta --m8_file ./output/esmfold_filtered_sequences.m8 --identity_threshold 0.5
```

- Molecular docking:

```bash
conda activate vina
mkdir molecular_docking
cd molecular_docking
cp ../src/molecular_docking.ipynb .
```

then run `molecular_docking/molecular_docking.ipynb` step by step.

## Acknowledgements

This work is inspired by Prof. Qingsong Wang and Wenyuan Zhu, on behalf of Biochemistry Lab teaching team at Peking University.

We thank the Computing Platform of the Center for Life Science (Peking University) for providing resources for the GPU-based model inference. Part of the computation was performed on the computing platform of the Infinite Intelligence Pharma Ltd.

If you have any questions, please contact me through my e-mail: [yiyanliao@stu.pku.edu.cn](yiyanliao@stu.pku.edu.cn).

For more of my research, refer to [https://yiyanliao.github.io](https://yiyanliao.github.io).