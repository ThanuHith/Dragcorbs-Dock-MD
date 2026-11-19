"""Brand theme for the Dragcorbs dashboard."""
import streamlit as st

WHITE = "#FFFFFF"
BLUE = "#0052CC"
BLACK = "#0A0F1C"

CSS = f"""
<style>
html, body, [data-testid="stAppViewContainer"] {{
    background-color: {WHITE};
    color: {BLACK};
}}
[data-testid="stSidebar"] {{
    background-color: {WHITE};
    border-right: 1px solid #E3E8F0;
}}
h1, h2, h3, h4 {{ color: {BLACK}; font-weight: 700; }}
.stButton > button, .stDownloadButton > button {{
    background-color: {BLUE};
    color: {WHITE};
    border: none;
    border-radius: 8px;
    font-weight: 600;
}}
.stButton > button:hover, .stDownloadButton > button:hover {{
    background-color: #003F9E;
    color: {WHITE};
    box-shadow: 0 0 12px rgba(0, 82, 204, 0.35);
}}
.dc-brand {{
    display: flex; align-items: center; gap: 12px;
    padding: 10px; border-radius: 12px;
    box-shadow: 0 0 18px rgba(0, 82, 204, 0.18);
}}
.dc-brand img {{ width: 56px; height: 56px; border-radius: 10px; }}
.dc-name {{ font-size: 1.5rem; font-weight: 800; color: {BLACK}; margin: 0; }}
.dc-tag {{ font-size: 0.7rem; color: #5A6478; letter-spacing: 1px; }}
.dc-badge {{
    display: inline-block; margin-top: 4px;
    padding: 2px 10px; border-radius: 999px;
    background: rgba(0, 82, 204, 0.08); color: {BLUE};
    font-size: 0.7rem; font-weight: 700;
}}
.dc-card {{
    border: 1px solid #E3E8F0; border-left: 4px solid {BLUE};
    border-radius: 10px; padding: 14px; margin-bottom: 10px;
    background: {WHITE};
}}
.dc-card h4 {{ margin: 0 0 6px 0; }}
.dc-viewer {{
    border: 2px dashed {BLUE}; border-radius: 12px;
    min-height: 420px; display: flex; align-items: center;
    justify-content: center; color: #5A6478;
}}
</style>
"""


def apply_theme() -> None:
    st.markdown(CSS, unsafe_allow_html=True)
