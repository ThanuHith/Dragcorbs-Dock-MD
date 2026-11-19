"""Input Hub component for the Dragcorbs dashboard."""
from __future__ import annotations

import io
import pandas as pd
import streamlit as st

DEFAULT_STRUCTURE_ID = "5V5Z"
DEFAULT_LIGAND_SMILES = "OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F"

BENCHMARK_PRESETS = {
    "Custom / Manual Input": None,
    "5V5Z — CYP51 Lanosterol 14α-Demethylase + Fluconazole": (
        "5V5Z",
        "DC-FLUCONAZOLE",
        "OC(Cn1cncn1)(Cn2cncn2)c3ccc(F)cc3F",
    ),
    "1IEP — Human Abl Kinase Domain + Imatinib": (
        "1IEP",
        "DC-IMATINIB",
        "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc4nccc(n4)c5cccnc5",
    ),
    "1IEP — Human Abl Kinase Domain + Aspirin": (
        "1IEP",
        "DC-ASPIRIN",
        "CC(=O)OC1=CC=CC=C1C(=O)O",
    ),
    "3ERT — Estrogen Receptor Alpha + Tamoxifen": (
        "3ERT",
        "DC-TAMOXIFEN",
        "CC/C(=C(\\c1ccccc1)/c2ccc(OCCN(C)C)cc2)/c3ccccc3",
    ),
}


def _parse_uploaded_ligand(uploaded_file) -> tuple[int, list[dict]]:
    """Extract ligand count and structured ligand payloads (compound_id, smiles, ligand_pdb_block)."""
    ext = uploaded_file.name.rsplit(".", 1)[-1].lower()
    stem = uploaded_file.name.rsplit(".", 1)[0]
    raw_bytes = uploaded_file.getvalue()
    text = raw_bytes.decode("utf-8", errors="ignore")

    if ext == "csv":
        try:
            df = pd.read_csv(io.BytesIO(raw_bytes))
            smiles_col = next(
                (c for c in df.columns if c.strip().lower() in ("smiles", "canonical_smiles", "smi")),
                None,
            )
            id_col = next(
                (c for c in df.columns if c.strip().lower() in ("id", "compound_id", "name", "title")),
                None,
            )
            ligands = []
            for idx, row in df.head(25).iterrows():
                smi = str(row[smiles_col]).strip() if smiles_col else DEFAULT_LIGAND_SMILES
                cid = str(row[id_col]).strip() if id_col else f"{stem}-{idx + 1}"
                ligands.append({"compound_id": cid, "smiles": smi, "ligand_pdb_block": None})
            return len(df), ligands
        except Exception:
            return 0, []

    if ext == "sdf":
        blocks = [b.strip() for b in text.split("$$$$") if b.strip()]
        ligands = []
        for idx, block in enumerate(blocks[:20], start=1):
            first_line = block.splitlines()[0].strip() if block.splitlines() else ""
            cid = first_line if (first_line and len(first_line) <= 24) else f"{stem}-{idx}"
            ligands.append(
                {
                    "compound_id": cid,
                    "smiles": st.session_state.get("active_ligand_smiles", DEFAULT_LIGAND_SMILES),
                    "ligand_pdb_block": block + "\n",
                }
            )
        return len(blocks), ligands

    return 1, [
        {
            "compound_id": stem,
            "smiles": st.session_state.get("active_ligand_smiles", DEFAULT_LIGAND_SMILES),
            "ligand_pdb_block": text,
        }
    ]


def render_input_hub() -> dict:
    """Render the form container, benchmark selector, file uploader, and summary metric cards."""
    st.subheader("Target & Input Hub")

    if "active_structure_id" not in st.session_state:
        st.session_state["active_structure_id"] = DEFAULT_STRUCTURE_ID
    if "structure_status" not in st.session_state:
        st.session_state["structure_status"] = f"Ready ({DEFAULT_STRUCTURE_ID})"
    if "uploaded_records" not in st.session_state:
        st.session_state["uploaded_records"] = []
    if "active_ligand_smiles" not in st.session_state:
        st.session_state["active_ligand_smiles"] = DEFAULT_LIGAND_SMILES
    if "active_compound_id" not in st.session_state:
        st.session_state["active_compound_id"] = "DC-101"

    with st.container():
        preset_choice = st.selectbox(
            "Quick Validated Crystal Benchmark Presets (Optional)",
            options=list(BENCHMARK_PRESETS.keys()),
            index=0,
        )
        if BENCHMARK_PRESETS[preset_choice] is not None:
            p_pdb, p_cid, p_smi = BENCHMARK_PRESETS[preset_choice]
            if st.button("Apply Benchmark Preset", use_container_width=False):
                st.session_state["active_structure_id"] = p_pdb
                st.session_state["structure_status"] = f"RCSB Benchmark ({p_pdb})"
                st.session_state["active_compound_id"] = p_cid
                st.session_state["active_ligand_smiles"] = p_smi
                st.session_state.pop("custom_structure_text", None)
                st.session_state.pop("detected_grid_center", None)
                st.success(f"Loaded benchmark target {p_pdb} with ligand {p_cid}!")

        col_left, col_right = st.columns(2, gap="large")

        with col_left:
            st.markdown("#### 1. RCSB PDB Structure Identifier")
            structure_id = (
                st.text_input(
                    "Structure ID (Any 4-character RCSB PDB code, e.g., 5V5Z, 1IEP, 3ERT)",
                    value=st.session_state["active_structure_id"],
                    max_chars=4,
                    help="Automatically downloads the full crystal structure from RCSB PDB.",
                )
                .strip()
                .upper()
            )

            if st.button("Load Structure", use_container_width=True):
                if len(structure_id) == 4 and structure_id.isalnum():
                    st.session_state["active_structure_id"] = structure_id
                    st.session_state["structure_status"] = f"RCSB Ready ({structure_id})"
                    st.session_state.pop("custom_structure_text", None)
                    st.session_state.pop("detected_grid_center", None)
                    st.success(f"Target structure {structure_id} queued for automatic RCSB fetch & pocket detection.")
                else:
                    st.error("Please enter a valid 4-character alphanumeric PDB ID.")

            custom_structure = st.file_uploader(
                "Or upload a custom receptor structure file (.pdb)",
                type=["pdb"],
                key="custom_structure_uploader",
            )
            if custom_structure is not None:
                st.session_state["active_structure_id"] = custom_structure.name[:4].upper()
                st.session_state["structure_status"] = f"Custom ({custom_structure.name})"
                st.session_state["custom_structure_text"] = custom_structure.getvalue().decode(
                    "utf-8", errors="ignore"
                )

        with col_right:
            st.markdown("#### 2. Input Library & Ligand SMILES")
            cid_val = st.text_input(
                "Active Compound ID",
                value=st.session_state.get("active_compound_id", "DC-101"),
            ).strip()
            if cid_val:
                st.session_state["active_compound_id"] = cid_val

            smiles_val = st.text_input(
                "Active Ligand SMILES (for AutoDock Vina & Complex MD)",
                value=st.session_state["active_ligand_smiles"],
            ).strip()
            if smiles_val:
                st.session_state["active_ligand_smiles"] = smiles_val

            uploaded_files = st.file_uploader(
                "Upload batch ligand files (.sdf, .pdb, .csv)",
                type=["sdf", "pdb", "csv"],
                accept_multiple_files=True,
                key="library_file_uploader",
            )

            records_summary = []
            parsed_ligands: list[dict] = []
            if uploaded_files:
                for item in uploaded_files:
                    entry_count, lig_list = _parse_uploaded_ligand(item)
                    parsed_ligands.extend(lig_list)
                    records_summary.append(
                        {
                            "Filename": item.name,
                            "Format": item.name.rsplit(".", 1)[-1].upper(),
                            "Size (KB)": round(item.size / 1024, 2),
                            "Entries": entry_count,
                        }
                    )
                st.session_state["uploaded_records"] = records_summary
                st.session_state["uploaded_ligands"] = parsed_ligands
                st.dataframe(
                    pd.DataFrame(records_summary),
                    use_container_width=True,
                    hide_index=True,
                )
            else:
                st.session_state["uploaded_records"] = []
                st.session_state["uploaded_ligands"] = [
                    {
                        "compound_id": st.session_state.get("active_compound_id", "DC-101"),
                        "smiles": st.session_state["active_ligand_smiles"],
                        "ligand_pdb_block": None,
                    }
                ]

    total_entries = sum(r["Entries"] for r in st.session_state["uploaded_records"]) or 1
    total_files = len(st.session_state["uploaded_records"])

    st.markdown("#### 3. Summary Status")
    card_col1, card_col2 = st.columns(2, gap="large")

    with card_col1:
        st.markdown(
            f"""
            <div class="dc-card">
                <h4>Target Structure Status</h4>
                <p style="font-size: 1.35rem; font-weight: 700; color: #0052CC; margin: 4px 0;">
                    {st.session_state["structure_status"]}
                </p>
                <span style="font-size: 0.82rem; color: #5A6478;">
                    Active ID: {st.session_state["active_structure_id"]}
                </span>
            </div>
            """,
            unsafe_allow_html=True,
        )

    with card_col2:
        st.markdown(
            f"""
            <div class="dc-card">
                <h4>Active Ligand Library</h4>
                <p style="font-size: 1.35rem; font-weight: 700; color: #0052CC; margin: 4px 0;">
                    {total_entries} Compound(s) Ready
                </p>
                <span style="font-size: 0.82rem; color: #5A6478;">
                    Primary Ligand: {st.session_state.get("active_compound_id", "DC-101")}
                </span>
            </div>
            """,
            unsafe_allow_html=True,
        )

    return {
        "structure_id": st.session_state["active_structure_id"],
        "structure_status": st.session_state["structure_status"],
        "total_files": total_files,
        "total_entries": total_entries,
        "files": st.session_state["uploaded_records"],
    }
