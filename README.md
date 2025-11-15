# *De novo* design of GST using ESM3 - Biochemistry Lab

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

- Env for [ESMFold](https://github.com/facebookresearch/esm), please use `yaml` file from [One-command-install-ESMfold](https://github.com/RazzyChen/One-command-install-ESMfold):

```bash
conda env create -f environment.yaml
```

- Env for [AutoDock-Vina](https://autodock-vina.readthedocs.io/en/latest/installation.html):

```bash
conda create -n vina python=3.10
conda activate vina
conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme
pip install -U numpy scipy rdkit vina meeko gemmi prody
```

- You may install [GROMACS](https://www.gromacs.org) by yourself.

## Downloads

Due to network issues in China, we have to download checkpoints for ESM3 beforehand and "hack" ESM3 source code (no need if you have better internet on your server, refer to [ESM3 repo](https://github.com/evolutionaryscale/esm)):

- ESM3 model download (you may try this several (or many) times):

```bash
hf auth login
hf download EvolutionaryScale/esm3-sm-open-v1 --local-dir esm3-sm-open-v1
```

- "Hack" ESM3 source code:

  change `/path/to/anaconda3/envs/esm/lib/python3.14/site-packages/esm/utils/constants/esm3.py` line `105`, from:

  ```python
  path = Path(snapshot_download(repo_id="EvolutionaryScale/esm3-sm-open-v1"))
  ```

  To:

  ```python
  path = Path("/path/to/esm3-sm-open-v1") 
  ```

- The you need to download `UniRef90` database (~ 40 GB):

```bash
mmseqs databases UniRef90 uniref90 tmp
```

## Usages

- Generate 1000 potential candidates compiled with structural metrics:

```bash
python src/ESM3_generate.py --num_sequences 1000 --batch_size 8
```

- Calculate the identity of candidates towards `UniRef90` and collect the ones `<40%` or have no hit:

```bash
mmseqs easy-search --alignment-mode 3 -s 7 generated_sequences_<time>.fasta uniref90 generated_sequences_<time>.m8 tmp 
```

- Molecular docking:

```bash
mkdir molecular_docking
cd molecular_docking
wget https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/448041/record/SDF?record_type=2d&response_type=display
mv SDF\?record_type\=2d GSH.sdf
conda activate vina
cd molecular_docking
mk_prepare_ligand.py -i GSH.sdf -o GSH.pdbqt

```

## Acknowledgements

This work is inspired by Prof. Qingsong Wang and Wenyuan Zhu, on behalf of Biochemistry Lab teaching team at Peking University.

We thank the Computing Platform of the Center for Life Science (Peking University) for providing resources for the GPU-based model inference. Part of the computation was performed on the computing platform of the Infinite Intelligence Pharma Ltd.

If you have any questions, please contact me through my e-mail: [yiyanliao@stu.pku.edu.cn](yiyanliao@stu.pku.edu.cn).

For more of my research, refer to [https://yiyanliao.github.io](https://yiyanliao.github.io).