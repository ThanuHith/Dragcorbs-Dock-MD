"""3D Visualization Viewport & Trajectory Player for the Dragcorbs dashboard."""
from __future__ import annotations

import json
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from components.controls import _generate_local_md_fallback


def _ensure_default_trajectory(active_id: str) -> None:
    """Seed a default trajectory in session state if none has been run yet."""
    if not st.session_state.get("md_trajectory_frames"):
        fallback = _generate_local_md_fallback(
            pdb_id=active_id,
            force_field="amber14-all.xml",
            water_model="amber14/tip3pfb.xml",
            steps=2500,
            num_frames=15,
            temperature_k=300.0,
            friction_coeff=1.0,
            minimize_energy=True,
        )
        st.session_state["md_result"] = fallback
        st.session_state["md_trajectory_frames"] = fallback["trajectory_frames_pdb"]
        st.session_state["md_telemetry"] = fallback["telemetry"]


def _build_3dmol_trajectory_html(
    pdb_id: str,
    style_mode: str,
    color_scheme: str,
    show_surface: bool,
    focus_pocket: bool,
    active_frame_pdb: str,
    all_frames_pdb: list[str],
    is_playing: bool,
) -> str:
    """Build an embedded 3Dmol.js viewer that renders the protein-ligand complex and animates MD frames."""
    surface_js = (
        'viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.45, color: "#0052CC"}, {hetflag: false});'
        if show_surface
        else ""
    )
    style_map = {
        "Cartoon": f'{{cartoon: {{color: "{color_scheme}"}}}}',
        "Sticks": f'{{stick: {{colorscheme: "{color_scheme}"}}}}',
        "Spheres": f'{{sphere: {{colorscheme: "{color_scheme}", scale: 0.35}}}}',
    }
    style_js = style_map.get(style_mode, style_map["Cartoon"])
    combined_multimodel = "\n".join(all_frames_pdb)
    frame_payload_js = json.dumps(combined_multimodel if is_playing else active_frame_pdb)
    animate_js = (
        'viewer.animate({loop: "forward", reps: 0, interval: 140});'
        if is_playing
        else ""
    )
    zoom_target_js = 'viewer.zoomTo({resn: "LIG"}); viewer.zoom(0.75);' if focus_pocket else "viewer.zoomTo();"
    has_full_receptor = "true" if active_frame_pdb.count("ATOM  ") > 50 else "false"

    return f"""
    <div id="mol-viewport" style="width: 100%; height: 460px; position: relative; border-radius: 10px; border: 2px solid #0052CC; background: #FFFFFF;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      (function() {{
        const elem = document.getElementById("mol-viewport");
        if (!window.$3Dmol) {{
          elem.innerHTML = '<div style="padding:24px;color:#0A0F1C;font-family:sans-serif;">3D Molecular Viewport Active (Offline Mode)</div>';
          return;
        }}
        const viewer = $3Dmol.createViewer(elem, {{ backgroundColor: "white" }});
        const mdFrameData = {frame_payload_js};
        const hasFullReceptor = {has_full_receptor};

        function styleComplex() {{
          viewer.setStyle({{hetflag: false}}, {style_js});
          // Highlight binding-pocket amino acid sidechains within 5.0 Å of the docked ligand
          viewer.addStyle(
            {{byres: true, within: {{distance: 5.0, sel: {{resn: "LIG"}}}}, hetflag: false}},
            {{stick: {{radius: 0.14, colorscheme: "orangeCarbon"}}}}
          );
          // Highlight the docked ligand in high-contrast cyan ball-and-stick
          viewer.setStyle(
            {{resn: "LIG"}},
            {{stick: {{radius: 0.24, colorscheme: "cyanCarbon"}}, sphere: {{scale: 0.30, colorscheme: "cyanCarbon"}}}}
          );
          {surface_js}
          {animate_js}
          {zoom_target_js}
          viewer.render();
        }}

        if (hasFullReceptor && mdFrameData && mdFrameData.trim().length > 0) {{
          viewer.addModelsAsFrames(mdFrameData, "pdb");
          styleComplex();
        }} else {{
          $3Dmol.download("pdb:{pdb_id}", viewer, {{}}, function() {{
            if (mdFrameData && mdFrameData.trim().length > 0) {{
              viewer.addModelsAsFrames(mdFrameData, "pdb");
            }}
            styleComplex();
          }});
        }}
      }})();
    </script>
    """


def render_viewer() -> dict:
    """Render the 3D visualization viewport, dynamic MD frame scrubber, and live telemetry."""
    st.subheader("3D Molecular Viewport & Complex Trajectory Scrubber")

    active_id = st.session_state.get("active_structure_id", "5V5Z")
    _ensure_default_trajectory(active_id)

    md_result: dict = st.session_state.get("md_result", {})
    frames_pdb: list[str] = st.session_state.get("md_trajectory_frames", [])
    telemetry: list[dict] = st.session_state.get("md_telemetry", [])
    total_frames = max(1, len(frames_pdb))

    if "viewer_frame" not in st.session_state or st.session_state["viewer_frame"] > total_frames:
        st.session_state["viewer_frame"] = 1
    if "viewer_playing" not in st.session_state:
        st.session_state["viewer_playing"] = False

    with st.container():
        ctrl_col1, ctrl_col2, ctrl_col3 = st.columns([2, 2, 3], gap="medium")

        with ctrl_col1:
            style_mode = (
                st.segmented_control(
                    "Representation Style",
                    options=["Cartoon", "Sticks", "Spheres"],
                    default="Cartoon",
                )
                or "Cartoon"
            )

        with ctrl_col2:
            color_theme = st.selectbox(
                "Color Theme",
                options=["spectrum", "chain", "ssJmol", "default"],
                index=0,
            )

        with ctrl_col3:
            toggles = st.columns(3)
            show_surface = toggles[0].toggle("Surface", value=False)
            focus_pocket = toggles[1].toggle("Zoom Pocket", value=True)
            live_3d = toggles[2].toggle("Interactive 3D", value=True)

        st.markdown("#### Complex MD Trajectory Playback Controls")
        pb_cols = st.columns([1, 1, 1, 1, 1, 4], gap="small")

        if pb_cols[0].button("⏮ First", use_container_width=True):
            st.session_state["viewer_frame"] = 1
            st.session_state["viewer_playing"] = False

        if pb_cols[1].button("⏪ Step -1", use_container_width=True):
            st.session_state["viewer_frame"] = max(1, st.session_state["viewer_frame"] - 1)
            st.session_state["viewer_playing"] = False

        if pb_cols[2].button("▶ Play", use_container_width=True):
            st.session_state["viewer_playing"] = True

        if pb_cols[3].button("⏸ Pause", use_container_width=True):
            st.session_state["viewer_playing"] = False

        if pb_cols[4].button("⏩ Step +1", use_container_width=True):
            st.session_state["viewer_frame"] = min(total_frames, st.session_state["viewer_frame"] + 1)
            st.session_state["viewer_playing"] = False

        frame_val = pb_cols[5].slider(
            "Scrub MD Trajectory Frame",
            min_value=1,
            max_value=total_frames,
            value=st.session_state["viewer_frame"],
        )
        st.session_state["viewer_frame"] = frame_val

        active_frame_pdb = frames_pdb[frame_val - 1] if frames_pdb else ""
        frame_metrics = (
            telemetry[frame_val - 1]
            if telemetry and (frame_val - 1) < len(telemetry)
            else {
                "frame": frame_val,
                "time_ps": round(frame_val * 0.5, 2),
                "potential_energy_kj_mol": -43756.0,
                "temperature_k": 300.0,
                "rmsd_angstrom": 0.05,
            }
        )
        playback_state = "Playing Loop" if st.session_state["viewer_playing"] else "Scrubbing Frame"
        engine_name = md_result.get("engine", "Dragcorbs All-Atom Pocket Langevin Engine")
        init_pe = md_result.get("initial_energy_kj_mol", -43753.0)
        min_pe = md_result.get("minimized_energy_kj_mol", -43756.0)

        st.markdown(
            f"""
            <div class="dc-card" style="border: 2px solid #0052CC; background-color: #FFFFFF; padding: 18px; margin-top: 12px;">
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap;">
                    <div>
                        <h4 style="color: #0A0F1C; margin: 0;">3D Molecular Viewport • {engine_name}</h4>
                        <p style="color: #5A6478; font-size: 0.88rem; margin: 4px 0 0 0;">
                            Docking pose (cyan sticks) & 5.0 Å interacting pocket residues (orange sticks) |
                            Minimization: <b>{init_pe} → {min_pe} kJ/mol</b>
                        </p>
                    </div>
                    <span class="dc-badge">
                        Target: {active_id} • Frame: {frame_val}/{total_frames} ({frame_metrics['time_ps']} ps) • {playback_state}
                    </span>
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        if live_3d:
            components.html(
                _build_3dmol_trajectory_html(
                    pdb_id=active_id,
                    style_mode=style_mode,
                    color_scheme=color_theme,
                    show_surface=show_surface,
                    focus_pocket=focus_pocket,
                    active_frame_pdb=active_frame_pdb,
                    all_frames_pdb=frames_pdb,
                    is_playing=st.session_state["viewer_playing"],
                ),
                height=475,
            )
        else:
            st.markdown(
                f"""
                <div class="dc-viewer">
                    <div style="text-align: center; padding: 24px;">
                        <h3 style="color: #0A0F1C; margin-bottom: 8px;">3D Molecular Viewport Active</h3>
                        <p style="color: #0052CC; font-weight: 600; margin-bottom: 6px;">
                            Target: {active_id} | Frame: {frame_val}/{total_frames} | RMSD: {frame_metrics['rmsd_angstrom']} Å
                        </p>
                    </div>
                </div>
                """,
                unsafe_allow_html=True,
            )

        st.markdown("#### Live Frame Telemetry & Time-Series")
        m_cols = st.columns(4, gap="medium")
        m_cols[0].metric("Simulation Time", f"{frame_metrics['time_ps']} ps", f"Frame {frame_val}/{total_frames}")
        m_cols[1].metric("Potential Energy", f"{frame_metrics['potential_energy_kj_mol']} kJ/mol")
        m_cols[2].metric("System Temperature", f"{frame_metrics['temperature_k']} K")
        m_cols[3].metric("Ligand Pose RMSD", f"{frame_metrics['rmsd_angstrom']} Å")

        if telemetry:
            df_tel = pd.DataFrame(telemetry).set_index("time_ps")
            chart_cols = st.columns(3, gap="medium")
            with chart_cols[0]:
                st.caption("Potential Energy (kJ/mol) vs. Time (ps)")
                st.line_chart(df_tel[["potential_energy_kj_mol"]], height=180)
            with chart_cols[1]:
                st.caption("Temperature (K) vs. Time (ps)")
                st.line_chart(df_tel[["temperature_k"]], height=180)
            with chart_cols[2]:
                st.caption("RMSD (Å) vs. Time (ps)")
                st.line_chart(df_tel[["rmsd_angstrom"]], height=180)

        with st.expander("Inspect Active Frame PDB Coordinates"):
            st.code(active_frame_pdb, language="text")

    return {
        "active_structure": active_id,
        "style_mode": style_mode,
        "color_theme": color_theme,
        "show_surface": show_surface,
        "current_frame": st.session_state["viewer_frame"],
        "total_frames": total_frames,
        "is_playing": st.session_state["viewer_playing"],
    }
