"""Sidebar branding block: logo, name, tagline and status badge."""
import base64
from pathlib import Path

import streamlit as st

LOGO_PATH = Path(__file__).resolve().parent.parent / "assets" / "logo.png"


def _logo_data_uri() -> str:
    if not LOGO_PATH.exists():
        return ""
    encoded = base64.b64encode(LOGO_PATH.read_bytes()).decode()
    return f"data:image/jpeg;base64,{encoded}"


def render_branding(status: str = "Online") -> None:
    logo = _logo_data_uri()
    img_tag = f'<img src="{logo}" alt="Dragcorbs logo"/>' if logo else ""
    st.sidebar.markdown(
        f"""
        <div class="dc-brand">
            {img_tag}
            <div>
                <p class="dc-name">Dragcorbs</p>
                <div class="dc-tag">BIOTECH &bull; AI &bull; DRUG DISCOVERY</div>
                <span class="dc-badge">&#9679; {status}</span>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.sidebar.markdown("---")
