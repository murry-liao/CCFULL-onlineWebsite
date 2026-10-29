#!/usr/bin/env python3
"""Physical-model workspace for the retained real-potential CCFULL programs."""

from __future__ import annotations

import csv
import io
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

from physical_model.interface.experimental_data import discover_datasets, load_dataset
from physical_model.interface.physical_runner import CODE_SPECS, run_selected_code
from physical_model.interface.input_parser import input_text, parse_input

ROOT = PROJECT_ROOT
DEFAULT_INPUT_PATH = ROOT / "data" / "input_templates" / "o16_sm154.inp"
BARRIER_DELTA_E = 2.5

st.set_page_config(
    page_title="CCFULL Workbench",
    page_icon="⚛️",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
    :root { --cc-ink:#14213d; --cc-muted:#5e6b7d; --cc-blue:#2457d6; --cc-line:rgba(36,87,214,.14); }
    [data-testid="stAppViewContainer"] { background:radial-gradient(circle at 8% 0%,rgba(36,87,214,.10),transparent 28rem),#f7f9fc; }
    [data-testid="stHeader"] { background:transparent; }
    .block-container { max-width:1420px; padding-top:.65rem; padding-bottom:4rem; }
    .ccfull-hero { text-align:center; padding:.15rem 1rem .7rem; }
    .ccfull-hero h1 { color:var(--cc-ink); font-size:clamp(2rem,4vw,3.15rem); font-weight:760; letter-spacing:-.045em; margin:.2rem 0 .45rem; padding:.25rem 0 0; }
    .ccfull-hero p { color:var(--cc-muted); font-size:1.04rem; margin:0 auto; max-width:760px; }
    .panel-label { border-bottom:1px solid var(--cc-line); color:var(--cc-blue); font-size:.76rem; font-weight:760; letter-spacing:.13em; margin:0 0 1rem; padding:.1rem 0 .7rem; text-transform:uppercase; }
    [data-testid="stExpander"] { background:rgba(255,255,255,.88); border-color:var(--cc-line); border-radius:.95rem; }
    [data-testid="stMetric"] { background:rgba(255,255,255,.88); border:1px solid var(--cc-line); border-radius:.9rem; min-height:7rem; padding:.9rem 1rem; }
    [data-testid="stMetricLabel"] { white-space:normal !important; overflow:visible !important; }
    [data-testid="stMetricValue"] { font-size:1.55rem !important; white-space:normal !important; overflow:visible !important; text-overflow:clip !important; }
    [data-testid="stNumberInputStepUp"], [data-testid="stNumberInputStepDown"] { display:none; }
    div.stButton > button[kind="primary"] { background:linear-gradient(105deg,var(--cc-blue),#3277e5); border:none; }
    h1,h2,h3 { color:var(--cc-ink); }
    </style>
    """,
    unsafe_allow_html=True,
)


def render_page_header() -> None:
    st.markdown(
        """
        <div class="ccfull-hero">
          <h1>CCFULL Workbench</h1>
          <p>Coupled-channels fusion calculations.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )


def rows_to_csv(rows: list[dict]) -> bytes:
    if not rows:
        return b""
    handle = io.StringIO()
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue().encode()


def experimental_data_controls(
    apro: int, zpro: int, atar: int, ztar: int, key: str,
    *, default_enabled: bool = False,
    default_paths: tuple[Path, ...] = (),
) -> tuple[list[dict], float]:
    datasets = discover_datasets(apro, zpro, atar, ztar)
    show = st.checkbox(
        "Include experimental fusion cross-section data",
        value=default_enabled,
        key=f"{key}_show_experiment",
    )
    if not show:
        return [], 1.0
    if not datasets:
        st.info("No matching experimental dataset was found.")
        return [], 1.0
    barrier_de = st.number_input(
        "Experimental barrier distribution dE (MeV)",
        min_value=0.05,
        value=BARRIER_DELTA_E,
        step=0.05,
        format="%.2f",
        key=f"{key}_barrier_de",
        help=(
            "Uses a centred three-point difference for D(E). When E ± dE is not measured, "
            "E·sigma is linearly interpolated. The quoted cross-section uncertainty is "
            "propagated assuming independent measurements."
        ),
    )
    labels = {item["label"]: item for item in datasets}
    preferred = {Path(path).resolve() for path in default_paths}
    defaults = [
        label for label, item in labels.items()
        if Path(item["path"]).resolve() in preferred
    ]
    if not defaults:
        defaults = [next(iter(labels))]
    selected = st.multiselect(
        "Experimental datasets", list(labels), default=defaults,
        key=f"{key}_experiment_sets",
    )
    rows = []
    for label in selected:
        rows.extend(load_dataset(labels[label]["path"]))
    return rows, float(barrier_de)


def _cross_section_value(row: dict) -> float:
    key = next(key for key in ("fusion_cross_section_mb", "absorption_cross_section_mb") if key in row)
    return float(row[key])


def cross_section_figure(direct=None, experiment=None, y_axis_title="Fusion cross section (mb)", direct_name="Physical model", baseline=None, baseline_name="Hagino original CCFULL baseline"):
    figure = go.Figure()
    if direct:
        figure.add_trace(go.Scatter(x=[row["energy_mev"] for row in direct], y=[_cross_section_value(row) for row in direct], mode="lines+markers", name=direct_name, line={"color":"#276fbf","width":2.5}))
    if baseline:
        figure.add_trace(go.Scatter(x=[row["energy_mev"] for row in baseline], y=[_cross_section_value(row) for row in baseline], mode="lines+markers", name=baseline_name, line={"color":"#009e73","width":2.5,"dash":"dot"}))
    if experiment:
        datasets = list(dict.fromkeys(row["dataset"] for row in experiment))
        for index, dataset in enumerate(datasets):
            rows = [row for row in experiment if row["dataset"] == dataset]
            figure.add_trace(go.Scatter(x=[row["energy_cm_mev"] for row in rows], y=[row["fusion_cross_section_mb"] for row in rows], error_y={"type":"data","array":[row["uncertainty_upper_mb"] for row in rows],"arrayminus":[row["uncertainty_lower_mb"] for row in rows]}, mode="markers", name=f"Experiment · {dataset}", marker={"color":["#222","#009e73","#cc79a7","#e69f00"][index%4],"size":7}))
    figure.update_layout(height=500, margin={"l":20,"r":20,"t":35,"b":20}, xaxis_title="E_cm (MeV)", yaxis_title=y_axis_title, yaxis_type="log", hovermode="x unified", legend={"orientation":"h"})
    return figure


def barrier_distribution_rows(direct: list[dict]) -> list[dict]:
    """Estimate D(E) = d²[E*sigma(E)]/dE² at interior energy points."""
    ordered = sorted(direct or [], key=lambda row: float(row["energy_mev"]))
    if len(ordered) < 3:
        return []
    energies = [float(row["energy_mev"]) for row in ordered]
    weighted_cross_sections = [
        energy * _cross_section_value(row)
        for energy, row in zip(energies, ordered)
    ]
    distribution = []
    for index in range(1, len(ordered) - 1):
        h_previous = energies[index] - energies[index - 1]
        h_next = energies[index + 1] - energies[index]
        if h_previous <= 0 or h_next <= 0:
            continue
        second_derivative = 2.0 * (
            weighted_cross_sections[index - 1]
            / (h_previous * (h_previous + h_next))
            - weighted_cross_sections[index] / (h_previous * h_next)
            + weighted_cross_sections[index + 1]
            / (h_next * (h_previous + h_next))
        )
        distribution.append({
            "energy_mev": energies[index],
            "barrier_distribution_mb_per_mev": second_derivative,
        })
    return distribution


def second_derivative(E, sigma, sigma_error, dE=BARRIER_DELTA_E):
    """Return centred D(E) values and propagated uncertainties."""
    E = np.asarray(E, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    sigma_error = np.asarray(sigma_error, dtype=float)
    if not (E.shape == sigma.shape == sigma_error.shape):
        raise ValueError("E, sigma and sigma_error must have the same length")
    if E.ndim != 1:
        raise ValueError("The input data must be one-dimensional")
    if np.any(np.diff(E) <= 0):
        raise ValueError("E must be strictly increasing")
    if dE <= 0:
        raise ValueError("dE must be positive")

    weighted_cross_section = E * sigma
    weighted_error = E * sigma_error
    distribution = np.full(E.shape, np.nan, dtype=float)
    uncertainty = np.full(E.shape, np.nan, dtype=float)

    def add_interpolation_weights(coefficients, target_E, scale):
        right = np.searchsorted(E, target_E)
        if (
            right < len(E)
            and np.isclose(E[right], target_E, rtol=0.0, atol=1.0e-10)
        ):
            coefficients[right] += scale
            return
        left = right - 1
        if left < 0 or right >= len(E):
            raise ValueError("Interpolation position is outside the data range")
        fraction = (target_E - E[left]) / (E[right] - E[left])
        coefficients[left] += scale * (1.0 - fraction)
        coefficients[right] += scale * fraction

    for index, center_E in enumerate(E):
        lower_E = center_E - dE
        upper_E = center_E + dE
        if lower_E < E[0] or upper_E > E[-1]:
            continue
        coefficients = np.zeros(len(E), dtype=float)
        add_interpolation_weights(coefficients, lower_E, 1.0 / dE**2)
        coefficients[index] += -2.0 / dE**2
        add_interpolation_weights(coefficients, upper_E, 1.0 / dE**2)
        distribution[index] = coefficients @ weighted_cross_section
        uncertainty[index] = np.sqrt(np.sum((coefficients * weighted_error) ** 2))
    return distribution, uncertainty


def experimental_barrier_distribution_rows(
    experiment: list[dict], barrier_de: float,
) -> list[dict]:
    """Derive experimental D(E) and uncertainty at measured energy centres."""
    if barrier_de <= 0:
        return []
    distribution = []
    datasets = list(dict.fromkeys(row["dataset"] for row in (experiment or [])))
    for dataset in datasets:
        source_rows = sorted(
            (row for row in experiment if row["dataset"] == dataset),
            key=lambda row: float(row["energy_cm_mev"]),
        )
        if len(source_rows) < 3:
            continue
        grouped = {}
        for row in source_rows:
            energy = float(row["energy_cm_mev"])
            sigma = float(row["fusion_cross_section_mb"])
            sigma_error = 0.5 * (
                abs(float(row["uncertainty_lower_mb"]))
                + abs(float(row["uncertainty_upper_mb"]))
            )
            grouped.setdefault(energy, []).append((sigma, sigma_error))
        energies = sorted(grouped)
        values = np.asarray([
            sum(value for value, _ in grouped[energy]) / len(grouped[energy])
            for energy in energies
        ])
        errors = np.asarray([
            np.sqrt(sum(error ** 2 for _, error in grouped[energy])) / len(grouped[energy])
            for energy in energies
        ])
        barrier_values, barrier_errors = second_derivative(
            energies, values, errors, dE=barrier_de,
        )
        for energy, value, error in zip(energies, barrier_values, barrier_errors):
            if not (np.isfinite(value) and np.isfinite(error)):
                continue
            distribution.append({
                "energy_mev": float(energy),
                "barrier_distribution_mb_per_mev": float(value),
                "barrier_distribution_uncertainty_mb_per_mev": float(error),
                "dataset": dataset,
            })
    return distribution


def barrier_distribution_figure(rows: list[dict], experimental: list[dict] | None = None):
    figure = go.Figure()
    if rows:
        figure.add_trace(go.Scatter(
            x=[row["energy_mev"] for row in rows],
            y=[row["barrier_distribution_mb_per_mev"] for row in rows],
            mode="lines",
            name="Physical model D(E)",
            line={"color": "#111111", "width": 2.2},
            #marker={"color": "#111111", "size": 7, "symbol": "circle"},
        ))
    if experimental:
        colors = ["#e69f00", "#009e73", "#56b4e9", "#d55e00"]
        datasets = list(dict.fromkeys(row["dataset"] for row in experimental))
        for index, dataset in enumerate(datasets):
            dataset_rows = [row for row in experimental if row["dataset"] == dataset]
            color = colors[index % len(colors)]
            figure.add_trace(go.Scatter(
                x=[row["energy_mev"] for row in dataset_rows],
                y=[row["barrier_distribution_mb_per_mev"] for row in dataset_rows],
                mode="markers",
                name=f"Experiment · {dataset}",
                marker={"color": color, "size": 7, "symbol": "diamond"},
                error_y={
                    "type": "data",
                    "array": [
                        row["barrier_distribution_uncertainty_mb_per_mev"]
                        for row in dataset_rows
                    ],
                    "visible": True,
                    "thickness": 0.8,
                    "width": 2,
                },
            ))
    figure.add_hline(y=0, line_dash="dot", line_color="#6b7280", line_width=1)
    figure.update_layout(
        height=420,
        margin={"l": 20, "r": 20, "t": 35, "b": 20},
        xaxis_title="E_cm (MeV)",
        yaxis_title="D(E) (mb / MeV)",
        hovermode="x unified",
        legend={"orientation": "h"},
    )
    return figure


def partial_figure(direct, energy):
    rows = [row for row in (direct or []) if abs(row["energy_mev"]-energy)<1e-9]
    figure = make_subplots(specs=[[{"secondary_y":True}]])
    if rows:
        figure.add_trace(go.Bar(x=[row["j"] for row in rows], y=[row["delta_sigma_dynamic_iwbc_mb"] for row in rows], name="Physical model Δσ_J", marker_color="#276fbf", opacity=.5), secondary_y=False)
        figure.add_trace(go.Scatter(x=[row["j"] for row in rows], y=[row["penetrability_dynamic_iwbc"] for row in rows], name="Physical model P_J", mode="lines"), secondary_y=True)
    figure.update_xaxes(title_text="J")
    figure.update_yaxes(title_text="Δσ_J (mb)", secondary_y=False)
    figure.update_yaxes(title_text="P_J", type="log", secondary_y=True)
    figure.update_layout(height=410, margin={"l":20,"r":20,"t":35,"b":20}, legend={"orientation":"h"})
    return figure


def four_inputs(title: str, values: list[float], labels: list[str], prefix: str,
                integer_positions=()) -> list[float | int]:
    st.markdown(f"**{title}**")
    columns = st.columns(len(values))
    result = []
    for index, (column, label, value) in enumerate(zip(columns, labels, values)):
        if index in integer_positions:
            result.append(column.number_input(label, value=int(value), step=1, key=f"{prefix}_{index}"))
        else:
            result.append(column.number_input(label, value=float(value), format="%.8g", key=f"{prefix}_{index}"))
    return result

def render_physical_mode():
    parameter_panel, result_panel = st.columns([0.50, 0.50], gap="large")
    with result_panel:
        st.markdown('<div class="panel-label">Calculation results</div>', unsafe_allow_html=True)
    with parameter_panel:
        st.markdown('<div class="panel-label">Parameter settings</div>', unsafe_allow_html=True)
        physical_code = "real_numerov"
        code_label = CODE_SPECS[physical_code]["label"]
        base = parse_input(DEFAULT_INPUT_PATH)
        prefix = "phys"

        with st.container(border=True):
            # Reaction system and nuclear radii
            c = st.columns(4)
            apro = c[0].number_input("A_pro", min_value=1, value=base["apro"], step=1, key=prefix+"_apro")
            zpro = c[1].number_input("Z_pro", min_value=1, value=base["zpro"], step=1, key=prefix+"_zpro")
            atar = c[2].number_input("A_tar", min_value=1, value=base["atar"], step=1, key=prefix+"_atar")
            ztar = c[3].number_input("Z_tar", min_value=1, value=base["ztar"], step=1, key=prefix+"_ztar")
            c = st.columns(4)
            r0p = c[0].number_input("R0P (fm)", value=base["r0p"], format="%.8g", key=prefix+"_r0p")
            projectile_mode = c[1].selectbox("Projectile coupling", [-1, 0, 1],
                                             index=[-1, 0, 1].index(base["projectile_mode"]),
                                             format_func=lambda x: {-1: "-1 · inert", 0: "0 · vibration", 1: "1 · rotation"}[x],
                                             key=prefix+"_pmode")
            r0t = c[2].number_input("R0T (fm)", value=base["r0t"], format="%.8g", key=prefix+"_r0t")
            target_modes = [-1, 0, 1]
            target_mode = c[3].selectbox("Target coupling", target_modes,
                                         index=target_modes.index(base["target_mode"]),
                                         format_func=lambda x: {-1: "-1 · inert", 0: "0 · vibration", 1: "1 · rotation"}[x],
                                         key=prefix+"_tmode")

            # Coupling parameters
            if target_mode == 1:
                target1_labels = [
                    "E2T (MeV)", "Beta2T", "Beta4T", "Beta6T", "NrotT"]
                target1_int = (4,)
                target1_defaults = list(base["target1"])
                if len(target1_defaults) == 4:
                    target1_defaults = (
                        target1_defaults[:3] + [0.0, target1_defaults[3]])
            elif target_mode == 0:
                target1_labels = ["OmegaT (MeV)", "BetaT", "LambdaT", "NphononT"]
                target1_int = (2, 3)
                target1_defaults = list(base["target1"][:4])
            else:
                target1_labels = ["target record 1", "target record 2", "target record 3", "target record 4"]
                target1_int = ()
                target1_defaults = list(base["target1"][:4])
            target1 = four_inputs(
                "Target Primary Excitation / Rotational Band",
                target1_defaults, target1_labels, prefix+"_t1", target1_int)
            target2 = four_inputs("Target Secondary Vibration (may be built on a rotational target)", base["target2"],
                                  ["OmegaT2 (MeV)", "BetaT2", "LambdaT2", "NphononT2"], prefix+"_t2", (2, 3))
            if int(target2[2]) == 3 and int(target2[3]) > 0:
                if target_mode == 1:
                    st.warning(
                        "LambdaT2 = 3 is present in the input, but the selected standalone Numerov/Hagino "
                        "physical programs do not activate a secondary vibration on a rotational target. "
                        "That coupling is not active in the selected standalone physical programs."
                    )
                else:
                    st.info("Octupole vibration enabled: LambdaT2 = 3.")
            if projectile_mode == 1:
                projectile_title = "Projectile Rotational Band"
                projectile_labels = [
                    "E2P (MeV · 2+ energy)",
                    "Beta2P (quadrupole)",
                    "Beta4P (hexadecapole)",
                    "NrotP (band levels)",
                ]
                projectile_int = (3,)
            elif projectile_mode == 0:
                projectile_title = "Projectile Vibrational Excitation"
                projectile_labels = [
                    "OmegaP (MeV · phonon energy)",
                    "BetaP (deformation)",
                    "LambdaP (multipolarity)",
                    "NphononP (phonons)",
                ]
                projectile_int = (2, 3)
            else:
                projectile_title = "Projectile Coupling Record · Inactive"
                projectile_labels = [
                    "OmegaP (MeV · inactive)",
                    "BetaP (inactive)",
                    "LambdaP (inactive)",
                    "NphononP (inactive)",
                ]
                projectile_int = (2, 3)
            projectile = four_inputs(
                projectile_title, base["projectile"], projectile_labels,
                prefix+"_p", projectile_int,
            )

            # Transfer coupling, nuclear potential and numerical grid
            c = st.columns(5)
            qtrans = c[0].number_input("Qtrans (MeV)", value=base["qtrans"], format="%.8g", key=prefix+"_q")
            ftr = c[1].number_input("Ftr", value=base["ftr"], format="%.8g", key=prefix+"_f")
            v0 = c[2].number_input("V0 (MeV)", min_value=1e-8, value=base["v0"], format="%.8g", key=prefix+"_v0")
            r0 = c[3].number_input("r0 (fm)", min_value=1e-8, value=base["r0"], format="%.8g", key=prefix+"_potr0")
            a0 = c[4].number_input("a0 (fm)", min_value=1e-8, value=base["a0"], format="%.8g", key=prefix+"_a0")
            c = st.columns(5)
            emin = c[0].number_input("Emin (MeV)", value=base["emin"], format="%.8g", key=prefix+"_emin")
            emax = c[1].number_input("Emax (MeV)", value=base["emax"], format="%.8g", key=prefix+"_emax")
            de = c[2].number_input("dE (MeV)", min_value=1e-6, value=base["de"], format="%.8g", key=prefix+"_de")
            rmax = c[3].number_input("Rmax (fm)", min_value=1e-6, value=base["rmax"], format="%.8g", key=prefix+"_rmax")
            dr = c[4].number_input("dr (fm)", min_value=1e-4, value=base["dr"], format="%.8g", key=prefix+"_dr")
            lmax, tail_window, tail_tolerance = base["lmax"], base["tail_window"], base["tail_tolerance"]
            st.caption("The selected standalone program uses its internal partial-wave loop (J ≤ 200) and convergence criterion.")

        config = dict(base, apro=apro, zpro=zpro, atar=atar, ztar=ztar, r0p=r0p, r0t=r0t,
                      projectile_mode=projectile_mode, target_mode=target_mode,
                      target1=target1, target2=target2, projectile=projectile,
                      qtrans=qtrans, ftr=ftr, v0=v0, r0=r0, a0=a0, emin=emin, emax=emax,
                      de=de, rmax=rmax, dr=dr, lmax=lmax, tail_window=tail_window,
                      tail_tolerance=tail_tolerance)
        experiment, experiment_barrier_de = experimental_data_controls(apro, zpro, atar, ztar, "physical")
        try:
            preview = input_text(config)
            download_columns = st.columns(1)
            download_columns[0].download_button("Download current ccfull.inp", preview.encode(), "ccfull.inp", "text/plain")
            input_ok = True
        except Exception as exc:
            st.error(str(exc)); input_ok = False
        if st.button("Run Physical Model", type="primary", disabled=not input_ok, key="run_physical"):
            try:
                with st.spinner(f"Running {code_label}…"):
                    direct, partials, info = run_selected_code(config, physical_code)
                    baseline = baseline_info = None
                    baseline, _, baseline_info = run_selected_code(config, "hagino_original")
                st.session_state["physical_result"] = {
                    "config": config, "direct": direct, "partials": partials, "info": info,
                    "baseline": baseline, "baseline_info": baseline_info,
                }
            except Exception as exc:
                st.error(f"Physical-model calculation failed: {exc}")
    with result_panel:
        result = st.session_state.get("physical_result")
        if result and result["info"].get("code") == physical_code:
            info = result["info"]
            baseline = result.get("baseline")
            baseline_info = result.get("baseline_info")
            relative = []
            if baseline:
                baseline_by_energy = {round(row["energy_mev"], 10): _cross_section_value(row) for row in baseline}
                for row in result["direct"]:
                    reference = baseline_by_energy.get(round(row["energy_mev"], 10))
                    if reference is not None:
                        relative.append(abs(_cross_section_value(row) - reference) / max(abs(reference), 1e-12))
            with st.container(border=True):
                st.plotly_chart(cross_section_figure(
                    direct=result["direct"], experiment=experiment, y_axis_title="Fusion cross section (mb)",
                    direct_name=info["code_label"], baseline=baseline,
                ), width="stretch", config={"displaylogo": False})
            barrier_rows = barrier_distribution_rows(result["direct"])
            experimental_barrier_rows = experimental_barrier_distribution_rows(
                experiment, experiment_barrier_de,
            )
            if barrier_rows or experimental_barrier_rows:
                with st.container(border=True):
                    st.plotly_chart(
                        barrier_distribution_figure(barrier_rows, experimental_barrier_rows),
                        width="stretch", config={"displaylogo": False},
                    )
                    captions = [
                        "Physical-model D(E) uses D(E) = d²[E·σ(E)]/dE²; endpoints are omitted and no smoothing is applied."
                    ]
                    if experimental_barrier_rows:
                        captions.append(
                            f"Experimental curves use linear interpolation at dE = {experiment_barrier_de:g} MeV; "
                            "their error bars propagate the quoted cross-section uncertainties assuming independent measurements."
                        )
                    st.caption(" ".join(captions))
            else:
                st.info("At least three energy points are required to derive a barrier distribution.")
            c = st.columns(2)
            c[0].metric("Solver time", f"{info['wall_s']:.3f} s")
            c[1].metric("Baseline time", "—" if not baseline_info else f"{baseline_info['wall_s']:.3f} s")
            c = st.columns(2)
            c[0].metric("Max relative error", "—" if not relative else f"{100*max(relative):.3f}%")
            c[1].metric("Energy points", len(result["direct"]))
            st.caption(f"Physical code used: {info['code_label']} · Observable: {info['observable']}")
            table = pd.DataFrame(result["direct"])
            if baseline:
                baseline_table = pd.DataFrame({
                    "energy_mev": [row["energy_mev"] for row in baseline],
                    "hagino_baseline_mb": [_cross_section_value(row) for row in baseline],
                })
                table = table.merge(baseline_table, on="energy_mev", how="left")
                selected_values = [_cross_section_value(row) for row in result["direct"]]
                table["selected_solver_mb"] = selected_values
                table["relative_difference_vs_hagino"] = (
                    (table["selected_solver_mb"] - table["hagino_baseline_mb"]).abs()
                    / table["hagino_baseline_mb"].abs().clip(lower=1e-12)
                )
            st.dataframe(table, width="stretch", hide_index=True)
            if result["partials"]:
                energy = st.selectbox("Partial-wave diagnostic energy", [row["energy_mev"] for row in result["direct"]], key="physical_partial_energy")
                st.plotly_chart(partial_figure(result["partials"], energy), width="stretch", config={"displaylogo": False})
            else:
                st.info("The selected standalone program reports the total fusion cross section, mean angular momentum, and s-wave probability; it does not export individual partial-wave rows.")
            c = st.columns(4 if baseline else 3)
            c[0].download_button("Cross section CSV", rows_to_csv(result["direct"]), "physical_cross_sections.csv", "text/csv", width="stretch")
            next_column = 1
            if baseline:
                c[next_column].download_button("Hagino baseline CSV", rows_to_csv(baseline), "hagino_baseline.csv", "text/csv", width="stretch")
                next_column += 1
            c[next_column].download_button("Partial-wave CSV", rows_to_csv(result["partials"]), "physical_partial_waves.csv", "text/csv", width="stretch", disabled=not result["partials"])
            c[next_column + 1].download_button("Input Used", info["input"].encode(), "ccfull_used.inp", "text/plain", width="stretch")
        else:
            st.info("Set the physical parameters on the left, then run the solver to display the calculation here.")


render_page_header()
render_physical_mode()
