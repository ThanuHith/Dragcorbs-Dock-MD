"""Analytics & Export component for the Dragcorbs dashboard."""
from __future__ import annotations

import io
import json
import zipfile
import pandas as pd
import streamlit as st

from reports import generate_pdf_report, generate_word_report


DEFAULT_SCORES = [
    {
        "Pose Rank": 1,
        "Compound ID": "DC-101",
        "Binding Score (kcal/mol)": -9.4,
        "RMSD Lower (Å)": 0.00,
        "RMSD Upper (Å)": 0.00,
        "H-Bonds": 4,
        "Interacting Residues": "GLU286(A), THR315(A), MET318(A), ASP381(A)",
        "Status": "Top Hit",
    },
    {
        "Pose Rank": 2,
        "Compound ID": "DC-104",
        "Binding Score (kcal/mol)": -8.8,
        "RMSD Lower (Å)": 1.12,
        "RMSD Upper (Å)": 1.85,
        "H-Bonds": 3,
        "Interacting Residues": "GLU286(A), VAL289(A), MET290(A), ILE293(A)",
        "Status": "Candidate",
    },
    {
        "Pose Rank": 3,
        "Compound ID": "DC-109",
        "Binding Score (kcal/mol)": -8.3,
        "RMSD Lower (Å)": 1.45,
        "RMSD Upper (Å)": 2.30,
        "H-Bonds": 2,
        "Interacting Residues": "VAL289(A), ILE293(A), LEU298(A), ALA380(A)",
        "Status": "Candidate",
    },
]

DEFAULT_PROPERTIES = [
    {"title": "Molecular Weight", "value": "384.42 Da", "rule": "Threshold <= 500 Da", "status": "Optimal"},
    {"title": "LogP (Lipophilicity)", "value": "2.85", "rule": "Threshold <= 5.0", "status": "Optimal"},
    {"title": "TPSA", "value": "74.60 Å²", "rule": "Range 20 - 130 Å²", "status": "Optimal"},
    {"title": "H-Bond Donors / Acceptors", "value": "2 / 6", "rule": "HBD <= 5, HBA <= 10", "status": "Pass"},
    {"title": "Rotatable Bonds", "value": "5", "rule": "Threshold <= 10", "status": "Pass"},
    {"title": "Rule-of-Five Compliance", "value": "0 Violations", "rule": "Max 1 Violation", "status": "Compliant"},
]


def _build_export_bundle(scores_df: pd.DataFrame, properties: list[dict], target_id: str) -> bytes:
    """Create a ZIP bundle containing CSV, JSON, and Markdown summary reports."""
    buffer = io.BytesIO()
    summary_md = [
        f"# Dragcorbs Analysis Report — Target {target_id}",
        "",
        "## 1. Calculated Binding Affinity Scores",
        scores_df.to_markdown(index=False) if hasattr(scores_df, "to_markdown") else scores_df.to_csv(index=False),
        "",
        "## 2. Physicochemical Property Profile",
    ]
    for prop in properties:
        summary_md.append(f"- **{prop['title']}**: {prop['value']} ({prop['rule']} — {prop['status']})")

    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("binding_scores.csv", scores_df.to_csv(index=False))
        zf.writestr(
            "physicochemical_profile.json",
            json.dumps({"target": target_id, "properties": properties}, indent=2),
        )
        zf.writestr("report_summary.md", "\n".join(summary_md))
        zf.writestr(
            "dragcorbs_report.pdf",
            generate_pdf_report(
                f"Dragcorbs ({target_id})",
                {"Active Target": target_id},
                properties,
                scores_df.to_dict(orient="records"),
            ),
        )
        zf.writestr(
            "dragcorbs_report.docx",
            generate_word_report(
                f"Dragcorbs ({target_id})",
                {"Active Target": target_id},
                properties,
                scores_df.to_dict(orient="records"),
            ),
        )

    buffer.seek(0)
    return buffer.getvalue()


def render_analytics() -> dict:
    """Render the analytics tables, physicochemical profile cards, and export actions."""
    st.subheader("Analytics & Property Profiling")

    target_id = st.session_state.get("active_structure_id", "5V5Z")
    pocket_src = st.session_state.get("pocket_source", "auto-detected")
    scores_df = pd.DataFrame(st.session_state.get("scores_data", DEFAULT_SCORES))
    properties = st.session_state.get("properties_data", DEFAULT_PROPERTIES)

    st.markdown(f"#### 1. AutoDock Vina v1.2.7 Binding Affinity Scores (Target: `{target_id}` | Pocket: `{pocket_src}`)")
    st.dataframe(
        scores_df,
        use_container_width=True,
        hide_index=True,
    )

    st.markdown("#### 2. RDKit Physicochemical Property Profile")
    grid_cols = st.columns(3, gap="medium")
    for idx, prop in enumerate(properties):
        with grid_cols[idx % 3]:
            st.markdown(
                f"""
                <div class="dc-card">
                    <h4>{prop["title"]}</h4>
                    <p style="font-size: 1.35rem; font-weight: 700; color: #0052CC; margin: 4px 0;">
                        {prop["value"]}
                    </p>
                    <div style="display: flex; justify-content: space-between; align-items: center;">
                        <span style="font-size: 0.78rem; color: #5A6478;">{prop["rule"]}</span>
                        <span class="dc-badge">{prop["status"]}</span>
                    </div>
                </div>
                """,
                unsafe_allow_html=True,
            )

    st.markdown("#### 3. Report & Data Export")
    export_cols = st.columns(4, gap="medium")

    params = {
        "Active Structure ID": target_id,
        "Pocket Source": pocket_src,
        "Structure Status": st.session_state.get("structure_status", f"Ready ({target_id})"),
        "Execution Mode": st.session_state.get("last_execution_config", {}).get(
            "execution_mode", "Single Evaluation"
        ),
    }
    table_records = scores_df.to_dict(orient="records")

    pdf_bytes = generate_pdf_report(
        project_title=f"Dragcorbs ({target_id})",
        parameters=params,
        metrics=properties,
        table_rows=table_records,
    )
    export_cols[0].download_button(
        label="Download PDF (.PDF)",
        data=pdf_bytes,
        file_name=f"dragcorbs_{target_id}_report.pdf",
        mime="application/pdf",
        use_container_width=True,
    )

    docx_bytes = generate_word_report(
        project_title=f"Dragcorbs ({target_id})",
        parameters=params,
        metrics=properties,
        table_rows=table_records,
    )
    export_cols[1].download_button(
        label="Download Word (.DOCX)",
        data=docx_bytes,
        file_name=f"dragcorbs_{target_id}_report.docx",
        mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
        use_container_width=True,
    )

    csv_bytes = scores_df.to_csv(index=False).encode("utf-8")
    export_cols[2].download_button(
        label="Download Scores (.CSV)",
        data=csv_bytes,
        file_name=f"dragcorbs_{target_id}_scores.csv",
        mime="text/csv",
        use_container_width=True,
    )

    bundle_bytes = _build_export_bundle(scores_df, properties, target_id)
    export_cols[3].download_button(
        label="Export Archive (.ZIP)",
        data=bundle_bytes,
        file_name=f"dragcorbs_{target_id}_report_bundle.zip",
        mime="application/zip",
        use_container_width=True,
    )

    return {
        "target_id": target_id,
        "top_score": float(scores_df["Binding Score (kcal/mol)"].min()),
        "num_poses": len(scores_df),
    }
