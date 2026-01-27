# HingeProt (Fortran) — Colab/Jupyter UI

This repository provides a **Google Colab / Jupyter Notebook** user interface (UI) for running the classic **HingeProt** pipeline (rigid parts + hinge prediction) based on **Elastic Network Models (GNM/ANM)**.

It supports:
- Running via **PDB code** (downloaded from RCSB) or **uploading your own PDB/ENT**
- **Multi-chain detection** and chain selection
- **CA-only** PDB inputs (coarse-grained workflow)
- Interactive **3D preview** (3Dmol.js / py3Dmol) and downloadable outputs

---

## 🚀 Run on Google Colab

Open the notebook here:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/enesemretas/hingeprot_fortran/blob/main/HingeProt.ipynb)

Direct link: https://colab.research.google.com/github/enesemretas/hingeprot_fortran/blob/main/HingeProt.ipynb

---

### Inputs
- **PDB Code** (e.g., `4cln`)  
  Downloads from RCSB and prepares the input.
- **Upload PDB/ENT**  
  Useful for custom structures or modified models.
- **Chain selection**
  Detects chain IDs and allows selecting one or more chains.
- **GNM/ANM cutoffs**
  Select from a list or enter custom values.

### Outputs
- **Rigid Parts table** (for slowest modes)
- **Hinge residues**
- **Short flexible fragments**
- **Mode files** (mode1/mode2) with download buttons
- **Mode trajectory viewer**

Outputs are stored under a run folder (default: `/content/hingeprot_runs` in Colab).

---

## 🧩 Requirements

### Recommended
- **Google Colab** (simplest; most dependencies are handled automatically)

### Local Jupyter (Linux recommended)
- Python 3.x
- Jupyter Notebook/Lab
- `pip install ipywidgets requests py3Dmol`
- `git`, `perl`
- The compiled HingeProt binaries available in the repository

The UI also includes an automated step to install `libg2c.so.0` on Debian/Ubuntu-like environments when needed.

---

🧠 Notes & limitations

Single-chain focus (original HingeProt behavior): the classic server/pipeline is designed primarily around single-chain analysis, although this UI allows selecting multiple chains.

CA-only inputs: supported and recommended for coarse-grained ENM workflows.

Missing residues: neglected in calculations; results reflect available coordinates.

Residue limit: HingeProt traditionally limits residue count (commonly cited as 999).


📌 Citation

If you use this work, please cite the original HingeProt paper:

Emekli U, Schneidman-Duhovny D, Wolfson HJ, Nussinov R, Haliloglu T. (2008) HingeProt: Automated Prediction of Hinges in Protein Structures. Proteins, 70(4):1219–1227.

@article{Emekli2008HingeProt,
  title={HingeProt: Automated Prediction of Hinges in Protein Structures},
  author={Emekli, U. and Schneidman-Duhovny, D. and Wolfson, H. J. and Nussinov, R. and Haliloglu, T.},
  journal={Proteins},
  volume={70},
  number={4},
  pages={1219--1227},
  year={2008}
}

