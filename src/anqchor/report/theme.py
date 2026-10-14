"""Colors and the Plotly base layout shared by every report.

Values follow a validated reference palette (categorical order is the CVD-safety
mechanism; never cycle it). Text never wears a series color; status colors are
reserved for verdicts (PASS/WARN/FAIL) and always come with a label.
"""

from __future__ import annotations

from typing import Any

CATEGORICAL_LIGHT = [
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#6250d6",
    "#e34948",
]
CATEGORICAL_DARK = [
    "#3987e5",
    "#d95926",
    "#199e70",
    "#c98500",
    "#d55181",
    "#008300",
    "#9085e9",
    "#e66767",
]
"""Fixed categorical order. Scatter-like charts (all pairs visible) use at most the first 3."""

SEQUENTIAL_BLUE = [
    "#cde2fb",
    "#b7d3f6",
    "#9ec5f4",
    "#86b6ef",
    "#6da7ec",
    "#5598e7",
    "#3987e5",
    "#2a78d6",
    "#256abf",
    "#1c5cab",
    "#184f95",
    "#104281",
    "#0d366b",
]
"""Blue ramp 100 → 700 for magnitude (heatmaps)."""

DIVERGING = {"neg": "#e34948", "mid_light": "#f0efec", "mid_dark": "#383835", "pos": "#2a78d6"}

STATUS = {"good": "#0ca30c", "warning": "#fab219", "serious": "#ec835a", "critical": "#d03b3b"}
"""Reserved for verdicts; always paired with a label or icon."""

WARNING_LINE = {"light": "#b07a00", "dark": STATUS["warning"]}
"""Thin warning-colored lines (a threshold rule): the amber fill color has too little
contrast for a 1-2 px line on the light surface (1.8:1); this darker amber has 3.6:1."""

LIGHT = {
    "surface": "#fcfcfb",
    "page": "#f9f9f7",
    "ink": "#0b0b0b",
    "ink2": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
}
DARK = {
    "surface": "#1a1a19",
    "page": "#0d0d0d",
    "ink": "#f0efec",
    "ink2": "#c3c2b7",
    "muted": "#898781",
    "grid": "#2c2c2a",
    "axis": "#383835",
}

FONT = 'system-ui, -apple-system, "Segoe UI", sans-serif'

# Match verdict groups → status color. NO_TEST is neutral (muted ink), never a status color.
GROUP_COLORS = {
    "PASS": STATUS["good"],
    "WARN": STATUS["warning"],
    "FAIL": STATUS["critical"],
    "NO_TEST": LIGHT["muted"],
}


def series_color(i: int) -> str:
    """Light-mode categorical color for series ``i`` (0-based). Dark mode is swapped in
    the browser via :data:`CATEGORICAL_LIGHT` → :data:`CATEGORICAL_DARK`."""
    if i >= len(CATEGORICAL_LIGHT):
        raise ValueError(
            "more than 8 series: fold extra series into 'Other' or use small multiples"
        )
    return CATEGORICAL_LIGHT[i]


def base_layout(**overrides: Any) -> dict[str, Any]:
    """Plotly layout dict in light mode; the page script re-themes it for dark mode."""
    axis = {
        "showgrid": True,
        "gridcolor": LIGHT["grid"],
        "gridwidth": 1,
        "zeroline": False,
        "linecolor": LIGHT["axis"],
        "showline": True,
        "ticks": "",
        "tickfont": {"color": LIGHT["muted"], "size": 12},
        "title": {"font": {"color": LIGHT["ink2"], "size": 13}},
        "automargin": True,
    }
    layout: dict[str, Any] = {
        "font": {"family": FONT, "color": LIGHT["ink"], "size": 13},
        "paper_bgcolor": LIGHT["surface"],
        "plot_bgcolor": LIGHT["surface"],
        "margin": {"l": 56, "r": 16, "t": 16, "b": 48},
        "xaxis": dict(axis),
        "yaxis": dict(axis),
        "hoverlabel": {
            "bgcolor": LIGHT["surface"],
            "bordercolor": LIGHT["axis"],
            "font": {"color": LIGHT["ink"], "family": FONT},
        },
        "legend": {
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "left",
            "x": 0,
            "font": {"color": LIGHT["ink2"]},
        },
        "bargap": 0.35,
        "hovermode": "closest",
    }
    for k, v in overrides.items():
        if isinstance(v, dict) and isinstance(layout.get(k), dict):
            layout[k] = {**layout[k], **v}
        else:
            layout[k] = v
    return layout


def css_variables() -> str:
    """CSS custom properties for both modes (OS preference and explicit toggle)."""

    def block(d: dict[str, str], cat: list[str]) -> str:
        lines = [f"  --{k}: {v};" for k, v in d.items()]
        lines += [f"  --series-{i + 1}: {c};" for i, c in enumerate(cat)]
        lines += [f"  --status-{k}: {v};" for k, v in STATUS.items()]
        mode = "light" if d is LIGHT else "dark"
        lines.append(f"  --status-warning-line: {WARNING_LINE[mode]};")
        return "\n".join(lines)

    return (
        ":root {\n  color-scheme: light;\n" + block(LIGHT, CATEGORICAL_LIGHT) + "\n}\n"
        '@media (prefers-color-scheme: dark) {\n  :root:not([data-theme="light"]) {\n  color-scheme: dark;\n'
        + block(DARK, CATEGORICAL_DARK)
        + "\n  }\n}\n"
        ':root[data-theme="dark"] {\n  color-scheme: dark;\n'
        + block(DARK, CATEGORICAL_DARK)
        + "\n}\n"
    )
