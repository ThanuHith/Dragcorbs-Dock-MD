"""Main Streamlit entry point for the Dragcorbs web application."""
from __future__ import annotations

from pathlib import Path
import streamlit as st

def _inject_head_seo_and_verification() -> None:
    """Inject Google Search Console verification, title, and SEO metadata into Streamlit's <head>."""
    try:
        index_path = Path(st.__file__).resolve().parent / "static" / "index.html"
        if not index_path.exists():
            return
        html = index_path.read_text(encoding="utf-8")
        token = "gqmCSWYuHBZSyz0bIU2ilGKEU2p921xwrL9vyZLErIc"
        if token in html:
            return
        head_tags = (
            f'<meta name="google-site-verification" content="{token}" />\n'
            '<meta name="description" content="Dragcorbs is a web-based structural biology platform for online molecular docking (AutoDock Vina v1.2.7), RDKit physicochemical profiling, and active-site molecular dynamics simulation." />\n'
            '<meta name="keywords" content="Dragcorbs, molecular docking web app, online AutoDock Vina, molecular dynamics simulation, protein ligand docking, RDKit, drug discovery" />\n'
            '<meta property="og:title" content="Dragcorbs — Molecular Docking & Molecular Dynamics Web Platform" />\n'
            '<meta property="og:description" content="Run genuine AutoDock Vina v1.2.7 molecular docking, RCSB PDB active-site detection, and all-atom molecular dynamics directly in your browser." />\n'
            '<meta property="og:type" content="website" />\n'
            '<script type="application/ld+json">'
            '{"@context":"https://schema.org","@type":"SoftwareApplication","name":"Dragcorbs",'
            '"applicationCategory":"ScientificApplication","operatingSystem":"Web Browser",'
            '"url":"https://dragcorbs-dock-md.streamlit.app/",'
            '"description":"Online molecular docking, physicochemical descriptor profiling, and active-site molecular dynamics web application powered by AutoDock Vina and RDKit.",'
            '"offers":{"@type":"Offer","price":"0","priceCurrency":"USD"}}'
            '</script>\n'
        )
        html = html.replace(
            "<title>Streamlit</title>",
            "<title>Dragcorbs — Online Molecular Docking & Molecular Dynamics Web App</title>",
        )
        if "</head>" in html:
            html = html.replace("</head>", f"{head_tags}</head>", 1)
            index_path.write_text(html, encoding="utf-8")
    except Exception:
        pass


_inject_head_seo_and_verification()

from theme import apply_theme
from components.branding import render_branding, _logo_data_uri
from components.input_hub import render_input_hub
from components.controls import render_controls
from components.viewer import render_viewer
from components.analytics import render_analytics, DEFAULT_SCORES, DEFAULT_PROPERTIES
from reports import generate_combined_archive

LOGO_FILE = Path(__file__).resolve().parent / "assets" / "logo.png"


def _init_session_state() -> None:
    """Initialize default session state keys so state flows cleanly across tabs."""
    defaults = {
        "active_structure_id": "5V5Z",
        "structure_status": "Ready (5V5Z)",
        "uploaded_records": [],
        "execution_completed": False,
        "last_execution_config": {
            "execution_mode": "Single Evaluation",
            "search_iterations": 16,
            "output_poses": 9,
            "compute_backend": "CUDA (GPU)",
        },
        "scores_data": DEFAULT_SCORES,
        "properties_data": DEFAULT_PROPERTIES,
    }
    for key, default_val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = default_val


def _render_top_header() -> None:
    """Render the top navigation header with the company logo and live status badge."""
    logo_uri = _logo_data_uri()
    img_html = (
        f'<img src="{logo_uri}" alt="Dragcorbs" '
        'style="width:64px;height:64px;border-radius:12px;object-fit:cover;'
        'box-shadow:0 0 14px rgba(0,82,204,0.25);"/>'
        if logo_uri
        else ""
    )
    exec_badge = (
        "Execution Complete" if st.session_state.get("execution_completed") else "Workspace Ready"
    )
    active_id = st.session_state.get("active_structure_id", "5V5Z")

    st.markdown(
        f"""
        <div style="display:flex;align-items:center;justify-content:space-between;flex-wrap:wrap;
                    padding:14px 18px;margin-bottom:18px;border-radius:14px;background:#FFFFFF;
                    border:1px solid #E3E8F0;box-shadow:0 0 20px rgba(0,82,204,0.12);">
            <div style="display:flex;align-items:center;gap:16px;">
                {img_html}
                <div>
                    <h1 style="margin:0;font-size:1.85rem;color:#0A0F1C;line-height:1.15;">
                        Dragcorbs
                    </h1>
                    <div style="font-size:0.78rem;color:#5A6478;letter-spacing:1.2px;font-weight:600;">
                        BIOTECH &bull; AI &bull; STRUCTURAL ANALYTICS PLATFORM
                    </div>
                </div>
            </div>
            <div style="display:flex;align-items:center;gap:10px;">
                <span class="dc-badge">Active Target: {active_id}</span>
                <span class="dc-badge">&#9679; {exec_badge}</span>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _render_sidebar_summary() -> str:
    """Render sidebar navigation options, live workspace summary, and quick report export."""
    render_branding(
        status="Active" if st.session_state.get("execution_completed") else "Online"
    )

    st.sidebar.markdown("#### Workspace Layout")
    layout_mode = st.sidebar.radio(
        "View Mode",
        options=["Modular Tabs", "Single-Page Dashboard"],
        index=0,
        label_visibility="collapsed",
    )

    st.sidebar.markdown("---")
    st.sidebar.markdown("#### Session Summary")
    active_id = st.session_state.get("active_structure_id", "5V5Z")
    files_count = len(st.session_state.get("uploaded_records", []))
    entries_count = sum(r.get("Entries", 0) for r in st.session_state.get("uploaded_records", []))
    last_cfg = st.session_state.get("last_execution_config", {})

    st.sidebar.markdown(
        f"""
        <div class="dc-card">
            <h4>Current Pipeline State</h4>
            <p style="margin:2px 0;font-size:0.88rem;"><b>Target ID:</b> {active_id}</p>
            <p style="margin:2px 0;font-size:0.88rem;"><b>Library:</b> {entries_count} entries ({files_count} files)</p>
            <p style="margin:2px 0;font-size:0.88rem;"><b>Mode:</b> {last_cfg.get("execution_mode", "Single Evaluation")}</p>
            <p style="margin:2px 0;font-size:0.88rem;"><b>Backend:</b> {last_cfg.get("compute_backend", "CUDA (GPU)")}</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    archive_bytes = generate_combined_archive(
        project_title=f"Dragcorbs ({active_id})",
        parameters={
            "Active Target": active_id,
            "Library Entries": entries_count,
            "Execution Mode": last_cfg.get("execution_mode", "Single Evaluation"),
            "Backend": last_cfg.get("compute_backend", "CUDA (GPU)"),
        },
        metrics=st.session_state.get("properties_data", DEFAULT_PROPERTIES),
        table_rows=st.session_state.get("scores_data", DEFAULT_SCORES),
    )
    st.sidebar.download_button(
        label="One-Click PDF & Word (.ZIP)",
        data=archive_bytes,
        file_name=f"dragcorbs_{active_id}_reports.zip",
        mime="application/zip",
        use_container_width=True,
    )

    return layout_mode


def main() -> None:
    """Configure page, apply styling, and render the Dragcorbs dashboard."""
    st.set_page_config(
        page_title="Dragcorbs — Online Molecular Docking & Molecular Dynamics Web App",
        page_icon=str(LOGO_FILE) if LOGO_FILE.exists() else "🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    _init_session_state()
    apply_theme()
    st.markdown(
        """
        <meta name="description" content="Dragcorbs is a web-based computational structural biology platform for online molecular docking (AutoDock Vina v1.2.7), RDKit physicochemical profiling, and active-site molecular dynamics simulation." />
        <meta name="keywords" content="Dragcorbs, molecular docking web app, online AutoDock Vina, molecular dynamics simulation, protein ligand docking, RDKit, drug discovery" />
        <meta property="og:title" content="Dragcorbs — Molecular Docking & Molecular Dynamics Web Platform" />
        <meta property="og:description" content="Run genuine AutoDock Vina v1.2.7 molecular docking, RCSB PDB active-site detection, and all-atom molecular dynamics directly in your browser." />
        <meta property="og:type" content="website" />
        <script type="application/ld+json">
        {
          "@context": "https://schema.org",
          "@type": "SoftwareApplication",
          "name": "Dragcorbs",
          "applicationCategory": "ScientificApplication",
          "operatingSystem": "Web Browser",
          "description": "Online molecular docking, physicochemical descriptor profiling, and active-site molecular dynamics web application powered by AutoDock Vina and RDKit.",
          "offers": { "@type": "Offer", "price": "0", "priceCurrency": "USD" }
        }
        </script>
        """,
        unsafe_allow_html=True,
    )

    layout_mode = _render_sidebar_summary()
    _render_top_header()

    if layout_mode == "Modular Tabs":
        tab_hub, tab_ctrl, tab_view, tab_analytics = st.tabs(
            [
                "1. Target & Ligand Hub",
                "2. Execution Control Panel",
                "3. 3D Molecular Viewer",
                "4. Analytics & Report Export",
            ]
        )
        with tab_hub:
            render_input_hub()
        with tab_ctrl:
            render_controls()
        with tab_view:
            render_viewer()
        with tab_analytics:
            render_analytics()
    else:
        render_input_hub()
        st.markdown("---")
        render_controls()
        st.markdown("---")
        render_viewer()
        st.markdown("---")
        render_analytics()


if __name__ == "__main__":
    main()



