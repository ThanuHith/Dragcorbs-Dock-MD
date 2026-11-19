---
title: Dragcorbs
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: docker
app_port: 8080
pinned: false
---

# Dragcorbs — Molecular Docking & Molecular Dynamics Web Platform

**Dragcorbs** is a full-stack computational structural biology web application for automated protein-ligand docking, physicochemical profiling, and active-site molecular dynamics.

## Core Capabilities
- **RCSB PDB Auto-Fetcher & Co-Crystal Pocket Detection**: Automatically retrieves crystal structures by 4-character PDB ID and computes active-site bounding boxes from bound co-crystal ligands.
- **Physics-Based Molecular Docking**: Powered by **AutoDock Vina v1.2.7**, **RDKit** (ETKDGv3 3D conformer embedding & MMFF94 optimization), **Meeko**, and **Gemmi**.
- **3D Active-Site Contact & H-Bond Analysis**: Identifies interatomic hydrogen bonds ($\le 3.5\text{ \AA}$) and interacting pocket residues.
- **All-Atom Pocket Molecular Dynamics**: Evaluates Lennard-Jones 12-6 + Coulomb interactions and Langevin dynamics trajectories with interactive **3Dmol.js** playback and **PDF/DOCX/CSV** report export.
