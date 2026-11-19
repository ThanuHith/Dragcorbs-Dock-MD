import os
import platform
import shutil
import re
import subprocess
import tempfile
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

BACKEND_DIR = Path(__file__).resolve().parent.parent
PDB_CACHE_DIR = BACKEND_DIR / "cache" / "pdb"
PDB_CACHE_DIR.mkdir(parents=True, exist_ok=True)
def _resolve_vina_binary() -> Path:
    if platform.system().lower() == "windows":
        return BACKEND_DIR / "vina.exe"
    for candidate in (BACKEND_DIR / "vina", Path("/tmp/vina")):
        if candidate.exists() and os.access(candidate, os.X_OK):
            return candidate
    sys_vina = shutil.which("vina")
    if sys_vina:
        return Path(sys_vina)
    tmp_vina = Path("/tmp/vina")
    url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64"
    urllib.request.urlretrieve(url, tmp_vina)
    os.chmod(tmp_vina, 0o755)
    return tmp_vina

VINA_BINARY = _resolve_vina_binary()

EXCLUDED_HETATMS = {
    "HOH", "WAT", "TIP", "SO4", "PO4", "GOL", "EDO", "PEG", "DMS", "ACT",
    "CL", "NA", "MG", "ZN", "CA", "K", "FE", "MN", "IOD", "BR", "TRS", "MPD"
}

AROMATIC_ATOMS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "CD2", "CE1"},
}


def fetch_rcsb_pdb(pdb_id: str) -> str:
    clean_id = pdb_id.strip().upper()
    if not re.match(r"^[0-9][A-Z0-9]{3}$", clean_id):
        raise ValueError(f"Invalid 4-character PDB ID: '{pdb_id}'")

    cache_file = PDB_CACHE_DIR / f"{clean_id}.pdb"
    if cache_file.exists() and cache_file.stat().st_size > 100:
        return cache_file.read_text(encoding="utf-8", errors="ignore")

    url = f"https://files.rcsb.org/download/{clean_id}.pdb"
    try:
        with urllib.request.urlopen(url, timeout=15) as response:
            pdb_text = response.read().decode("utf-8", errors="ignore")
        cache_file.write_text(pdb_text, encoding="utf-8")
        return pdb_text
    except Exception as exc:
        raise RuntimeError(f"Failed to download PDB {clean_id} from RCSB: {exc}") from exc


def detect_binding_pocket(
    pdb_content: str, target_chain: Optional[str] = None, padding: float = 6.0
) -> Dict[str, float]:
    het_clusters: Dict[str, List[Tuple[float, float, float]]] = {}
    receptor_coords: List[Tuple[float, float, float]] = []

    for line in pdb_content.splitlines():
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            chain_id = line[21:22].strip()
            if target_chain and chain_id and chain_id != target_chain:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            if line.startswith("ATOM  "):
                receptor_coords.append((x, y, z))
            else:
                res_name = line[17:20].strip().upper()
                res_seq = line[22:26].strip()
                if res_name not in EXCLUDED_HETATMS:
                    key = f"{res_name}_{chain_id}_{res_seq}"
                    het_clusters.setdefault(key, []).append((x, y, z))

    valid_ligands = {k: v for k, v in het_clusters.items() if len(v) >= 6}
    if valid_ligands:
        best_key = max(valid_ligands, key=lambda k: len(valid_ligands[k]))
        coords = np.array(valid_ligands[best_key])
        min_pt = coords.min(axis=0)
        max_pt = coords.max(axis=0)
        center = (min_pt + max_pt) / 2.0
        span = np.clip((max_pt - min_pt) + padding, 18.0, 30.0)
        return {
            "center_x": round(float(center[0]), 3),
            "center_y": round(float(center[1]), 3),
            "center_z": round(float(center[2]), 3),
            "size_x": round(float(span[0]), 2),
            "size_y": round(float(span[1]), 2),
            "size_z": round(float(span[2]), 2),
            "pocket_source": f"co_crystal_ligand:{best_key}",
        }

    if not receptor_coords:
        raise ValueError("Receptor PDB contains no valid ATOM coordinates.")

    coords = np.array(receptor_coords)
    center = coords.mean(axis=0)
    return {
        "center_x": round(float(center[0]), 3),
        "center_y": round(float(center[1]), 3),
        "center_z": round(float(center[2]), 3),
        "size_x": 22.0,
        "size_y": 22.0,
        "size_z": 22.0,
        "pocket_source": "receptor_centroid",
    }


def prepare_ligand_pdbqt(input_data: str, is_smiles: bool = True) -> str:
    if is_smiles:
        mol = Chem.MolFromSmiles(input_data.strip())
    else:
        mol = Chem.MolFromMolBlock(input_data, removeHs=False)
        if mol is None:
            mol = Chem.MolFromPDBBlock(input_data, removeHs=False)
    if mol is None:
        raise ValueError("Invalid ligand input: could not parse SMILES or structure block.")

    mol = Chem.AddHs(mol, addCoords=True)
    if is_smiles or mol.GetNumConformers() == 0:
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) != 0:
            AllChem.EmbedMolecule(mol, useRandomCoords=True)
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
    if not is_ok:
        raise RuntimeError(f"Meeko PDBQT conversion failed: {err_msg}")
    return pdbqt_string


def clean_receptor_to_pdbqt(pdb_content: str, target_chain: Optional[str] = None) -> str:
    pdbqt_lines = []
    for line in pdb_content.splitlines():
        if not line.startswith("ATOM  "):
            continue
        chain_id = line[21:22].strip()
        if target_chain and chain_id and chain_id != target_chain:
            continue

        alt_loc = line[16:17]
        if alt_loc not in (" ", "A", "1"):
            continue

        atom_name = line[12:16].strip()
        res_name = line[17:20].strip().upper()
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            element = re.sub(r"[0-9]", "", atom_name)[:1].upper()

        if element == "C":
            ad_type = "A" if atom_name in AROMATIC_ATOMS.get(res_name, set()) else "C"
        elif element == "N":
            ad_type = "NA" if (res_name == "HIS" and atom_name in ("ND1", "NE2")) else "N"
        elif element == "O":
            ad_type = "OA"
        elif element == "S":
            ad_type = "SA"
        elif element == "H":
            if atom_name.startswith(("HN", "HO", "HS", "HG1", "HH", "HE", "HD")) or atom_name == "H":
                ad_type = "HD"
            else:
                continue
        else:
            ad_type = element[:2]

        base_line = line[:54].ljust(54)
        pdbqt_lines.append(f"{base_line}  1.00  0.00    +0.000 {ad_type:<2}")

    if not pdbqt_lines:
        raise ValueError("No valid receptor ATOM lines remained after cleaning.")
    return "\n".join(pdbqt_lines) + "\n"


def pdbqt_pose_to_pdb_block(output_pdbqt: str, model_index: int = 1) -> str:
    in_target_model = False
    current_model = 0
    pdb_lines = []
    serial = 1

    for line in output_pdbqt.splitlines():
        if line.startswith("MODEL"):
            current_model += 1
            if current_model == model_index:
                in_target_model = True
            continue
        if line.startswith("ENDMDL"):
            if in_target_model:
                break
            continue

        if current_model == 0 and (line.startswith("ATOM  ") or line.startswith("HETATM")):
            in_target_model = True

        if in_target_model and (line.startswith("ATOM  ") or line.startswith("HETATM")):
            ad_type = line[77:79].strip() if len(line) >= 78 else line.split()[-1]
            elem_map = {"OA": "O", "NA": "N", "SA": "S", "HD": "H", "A": "C"}
            element = elem_map.get(ad_type, ad_type[:2].capitalize())
            atom_name = line[12:16]
            xyz_part = line[30:54]
            pdb_line = (
                f"HETATM{serial:5d} {atom_name} LIG L   1    "
                f"{xyz_part}  1.00  0.00          {element:>2}"
            )
            pdb_lines.append(pdb_line)
            serial += 1

    return "\n".join(pdb_lines)


def run_autodock_vina(
    receptor_pdb: str,
    ligand_input: str,
    is_smiles: bool = True,
    center_x: Optional[float] = None,
    center_y: Optional[float] = None,
    center_z: Optional[float] = None,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    exhaustiveness: int = 8,
    target_chain: Optional[str] = None,
) -> Dict:
    if not VINA_BINARY.exists():
        raise FileNotFoundError(f"AutoDock Vina binary not found at {VINA_BINARY}")

    pocket_info = detect_binding_pocket(receptor_pdb, target_chain=target_chain)
    cx = center_x if center_x is not None else pocket_info["center_x"]
    cy = center_y if center_y is not None else pocket_info["center_y"]
    cz = center_z if center_z is not None else pocket_info["center_z"]
    sx = size_x if center_x is not None else pocket_info["size_x"]
    sy = size_y if center_y is not None else pocket_info["size_y"]
    sz = size_z if center_z is not None else pocket_info["size_z"]

    receptor_pdbqt = clean_receptor_to_pdbqt(receptor_pdb, target_chain=target_chain)
    ligand_pdbqt = prepare_ligand_pdbqt(ligand_input, is_smiles=is_smiles)

    with tempfile.TemporaryDirectory() as tmpdir:
        rec_path = os.path.join(tmpdir, "receptor.pdbqt")
        lig_path = os.path.join(tmpdir, "ligand.pdbqt")
        out_path = os.path.join(tmpdir, "out.pdbqt")

        with open(rec_path, "w", encoding="utf-8") as f:
            f.write(receptor_pdbqt)
        with open(lig_path, "w", encoding="utf-8") as f:
            f.write(ligand_pdbqt)

        cmd = [
            str(VINA_BINARY),
            "--receptor", rec_path,
            "--ligand", lig_path,
            "--center_x", str(cx),
            "--center_y", str(cy),
            "--center_z", str(cz),
            "--size_x", str(sx),
            "--size_y", str(sy),
            "--size_z", str(sz),
            "--exhaustiveness", str(exhaustiveness),
            "--out", out_path,
        ]

        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise RuntimeError(f"Vina execution failed:\n{proc.stderr or proc.stdout}")

        output_pdbqt = ""
        if os.path.exists(out_path):
            with open(out_path, "r", encoding="utf-8") as f:
                output_pdbqt = f.read()

    poses = []
    for line in output_pdbqt.splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            parts = line.replace("REMARK VINA RESULT:", "").split()
            if len(parts) >= 3:
                poses.append({
                    "mode": len(poses) + 1,
                    "affinity_kcal_mol": float(parts[0]),
                    "rmsd_lb": float(parts[1]),
                    "rmsd_ub": float(parts[2]),
                })

    best_affinity = poses[0]["affinity_kcal_mol"] if poses else 0.0
    top_pose_pdb = pdbqt_pose_to_pdb_block(output_pdbqt, model_index=1)

    return {
        "best_affinity_kcal_mol": best_affinity,
        "poses": poses,
        "grid_box": {
            "center_x": cx, "center_y": cy, "center_z": cz,
            "size_x": sx, "size_y": sy, "size_z": sz,
            "pocket_source": pocket_info["pocket_source"],
        },
        "output_pdbqt": output_pdbqt,
        "top_pose_pdb": top_pose_pdb,
        "vina_stdout": proc.stdout,
    }


