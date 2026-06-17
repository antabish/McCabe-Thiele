"""
McCabe-Thiele Analysis: Benzene / Toluene Distillation
=======================================================
Example 7.1 from:
  Seader, J.D., Henley, E.J., & Roper, D.K.
  Separation Process Principles, 3rd Edition

Generates four publication-quality figures:
  Figure A — Minimum stages at total reflux
  Figure B — Minimum reflux condition (pinch point)
  Figure C — Building the McCabe-Thiele design framework
  Figure D — Actual operating design with stage stepping

Exports: PNG (600 dpi), SVG, and CSV of stage coordinates for PowerPoint.

Usage:
  python Example7.1.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

# Ensure UTF-8 output on Windows consoles that default to cp1252
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

# ─────────────────────────────────────────────────────────────────────────────
# OUTPUT DIRECTORY — same folder as this script
# ─────────────────────────────────────────────────────────────────────────────

OUT_DIR = Path(__file__).parent

# ─────────────────────────────────────────────────────────────────────────────
# GLOBAL PLOT STYLE — applied once; all figures inherit these settings
# ─────────────────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family":       "serif",
    "font.size":         18,
    "axes.linewidth":    1.2,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "lines.antialiased": True,
    "figure.dpi":        150,
})

# ─────────────────────────────────────────────────────────────────────────────
# COLOUR PALETTE
# ─────────────────────────────────────────────────────────────────────────────

C = dict(
    eq     = "#1a56db",   # deep blue  — equilibrium curve
    diag   = "#111111",   # near-black — 45° line
    qline  = "#d97706",   # amber      — q-line
    rect   = "#16a34a",   # green      — rectifying line
    strip  = "#9333ea",   # purple     — stripping line
    stages = "#dc2626",   # red        — stage steps
    pinch  = "#e11d48",   # rose       — pinch point marker
    spec   = "#64748b",   # slate      — xB / xF / xD verticals
)

# ─────────────────────────────────────────────────────────────────────────────
# 1.  VLE DATA — Benzene / Toluene at 1 atm  (Seader et al., Table 7.1)
# ─────────────────────────────────────────────────────────────────────────────

# Experimental data with pure-component endpoints (0,0) and (1,1) added.
_X_EXP = np.array([
    0.000, 0.014, 0.100, 0.200, 0.300, 0.400,
    0.500, 0.600, 0.700, 0.800, 0.962, 1.000,
])
_Y_EXP = np.array([
    0.000, 0.032, 0.207, 0.374, 0.509, 0.620,
    0.713, 0.790, 0.855, 0.911, 0.985, 1.000,
])

# Cubic spline through the experimental VLE data (monotone; one root per y)
_SPLINE = CubicSpline(_X_EXP, _Y_EXP)

# Dense grid for smooth curve rendering
_X_CURVE = np.linspace(0.0, 1.0, 1000)
_Y_CURVE = _SPLINE(_X_CURVE)


def y_eq(x: float | np.ndarray) -> float | np.ndarray:
    """Return equilibrium vapour composition y* for liquid composition x."""
    return _SPLINE(x)


def x_eq(y: float) -> float:
    """
    Invert the equilibrium curve: find x such that y_eq(x) == y.
    Uses brentq; valid because y_eq is strictly monotone on [0, 1].
    """
    return float(brentq(lambda xv: float(y_eq(xv)) - y, 0.0, 1.0, xtol=1e-10))


# ─────────────────────────────────────────────────────────────────────────────
# 2.  DESIGN SPECIFICATIONS  (Example 7.1)
# ─────────────────────────────────────────────────────────────────────────────

XF = 0.60   # feed composition           (mol fraction benzene)
XD = 0.95   # distillate specification   (mol fraction benzene)
XB = 0.05   # bottoms specification      (mol fraction benzene)
Q  = 0.61   # feed quality  (0 < Q < 1 → partially vaporised feed)

# ─────────────────────────────────────────────────────────────────────────────
# 3.  PROCESS CALCULATIONS
# ─────────────────────────────────────────────────────────────────────────────

def q_line(x: float | np.ndarray) -> float | np.ndarray:
    """
    q-line: y = [q/(q-1)]*x − xF/(q-1).
    Passes through (xF, xF); slope is negative for 0 < q < 1 (partial vapour).
    """
    return (Q / (Q - 1.0)) * x - XF / (Q - 1.0)


def _find_pinch_point() -> tuple[float, float]:
    """
    Pinch point = intersection of q-line with equilibrium curve.
    At minimum reflux the rectifying line passes through this point.
    brentq is used because the residual changes sign exactly once in (XB, XD).
    """
    x_p = brentq(
        lambda x: float(y_eq(x)) - float(q_line(x)),
        XB + 1e-6, XD - 1e-6,
    )
    return float(x_p), float(y_eq(x_p))


X_PINCH, Y_PINCH = _find_pinch_point()

# Slope of the minimum rectifying line through (xD, xD) and the pinch point
_M_RECT_MIN = (Y_PINCH - XD) / (X_PINCH - XD)

# Rmin derived from: slope = R/(R+1)  →  R = slope / (1 − slope)
RMIN = _M_RECT_MIN / (1.0 - _M_RECT_MIN)

# Operating reflux ratio (1.3 × Rmin as specified)
R = 1.3 * RMIN


def rectifying_line(x: float | np.ndarray) -> float | np.ndarray:
    """
    Rectifying (enriching section) operating line.
    Passes through (xD, xD) with slope R/(R+1).
    """
    return (R / (R + 1.0)) * x + XD / (R + 1.0)


def _find_feed_intersection() -> tuple[float, float]:
    """
    Intersection of the rectifying line with the q-line.
    Defines the column mid-point between rectifying and stripping sections.
    """
    x_f = brentq(
        lambda x: float(rectifying_line(x)) - float(q_line(x)),
        XB + 1e-6, XD - 1e-6,
    )
    return float(x_f), float(rectifying_line(x_f))


X_INT, Y_INT = _find_feed_intersection()

# Stripping line slope and intercept: passes through (xB, xB) and (X_INT, Y_INT)
_M_STRIP = (Y_INT - XB) / (X_INT - XB)
_B_STRIP = XB - _M_STRIP * XB          # ensures y = xB when x = xB


def stripping_line(x: float | np.ndarray) -> float | np.ndarray:
    """
    Stripping section operating line.
    Passes through (xB, xB) and the feed intersection point (X_INT, Y_INT).
    """
    return _M_STRIP * x + _B_STRIP


# ─────────────────────────────────────────────────────────────────────────────
# 4.  STAGE STEPPING
# ─────────────────────────────────────────────────────────────────────────────

# Segment list type: each entry is (x1, y1, x2, y2).
# Even-indexed entries (0, 2, …) are horizontal (equilibrium steps).
# Odd-indexed entries  (1, 3, …) are vertical   (operating-line steps).
Segments = list[tuple[float, float, float, float]]


def step_minimum_stages() -> tuple[Segments, int]:
    """
    Step off minimum stages at total reflux.

    At total reflux the 45° line (y = x) acts as the operating line.
    Stages are stepped between the equilibrium curve and the diagonal.

    Returns
    -------
    segs  : list of (x1, y1, x2, y2) segments
    n_min : minimum number of theoretical stages
    """
    segs: Segments = []
    x = XD
    n = 0

    while x > XB + 1e-6:
        y = x                                      # on the 45° line
        x_next = x_eq(y)                           # horizontal → equilibrium

        segs.append((x,      y, x_next, y     ))   # horizontal step
        segs.append((x_next, y, x_next, x_next))   # vertical step (to 45°)

        n += 1
        x = x_next

        if n > 200:
            raise RuntimeError(
                "Minimum-stage stepping did not converge. Check VLE data."
            )

    return segs, n


def step_actual_stages() -> tuple[Segments, int, int]:
    """
    Step off theoretical stages at the actual operating reflux ratio R.

    Switching rule (optimal feed stage):
      Use the rectifying line above the feed intersection (x > X_INT).
      Switch to the stripping line as soon as x_next ≤ X_INT.
      The stripping line is lower in this region, giving larger steps
      and therefore fewer total stages.

    Returns
    -------
    segs       : list of (x1, y1, x2, y2) segments
    n_stages   : total theoretical stages
    feed_stage : stage number at which the feed is introduced
    """
    segs: Segments = []
    x = XD
    y = XD          # rectifying_line(xD) == xD by construction
    n = 0
    feed_stage: int | None = None
    on_strip = False

    while x > XB + 1e-6:
        x_next = x_eq(y)                           # horizontal → equilibrium
        segs.append((x, y, x_next, y))
        n += 1

        # Switch to stripping once the equilibrium step crosses the feed point
        if not on_strip and x_next <= X_INT:
            on_strip = True
            feed_stage = n

        y_next = float(
            stripping_line(x_next) if on_strip else rectifying_line(x_next)
        )
        segs.append((x_next, y, x_next, y_next))  # vertical → operating line

        x = x_next
        y = y_next

        if n > 200:
            raise RuntimeError(
                "Actual-design stage stepping did not converge. "
                "Check design specifications."
            )

    return segs, n, (feed_stage if feed_stage is not None else n)


# ─────────────────────────────────────────────────────────────────────────────
# 5.  PLOTTING UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

def _make_axes(title: str) -> tuple[plt.Figure, plt.Axes]:
    """
    Create a square figure pre-populated with the VLE equilibrium curve,
    the 45° diagonal, and the experimental data points.
    All three figures share this common base.
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot(_X_CURVE, _Y_CURVE, color=C["eq"], lw=2.5,
            label="Equilibrium curve", zorder=3)
    ax.plot([0, 1], [0, 1], color=C["diag"], lw=1.5, ls="--",
            label="45° line", zorder=2)
    ax.scatter(_X_EXP, _Y_EXP, s=28, color=C["eq"], zorder=5,
               label="Experimental data")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Liquid mole fraction of benzene,  $x$", fontsize=16)
    ax.set_ylabel("Vapour mole fraction of benzene,  $y$", fontsize=16)
    """ax.set_title(title, fontsize=15, fontweight="bold", pad=12)"""
    ax.tick_params(labelsize=14)

    return fig, ax


def _draw_spec_verticals(ax: plt.Axes) -> None:
    """Dotted vertical reference lines and labels at xB, xF, xD."""
    for xv, lbl in [(XB, "$x_B$"), (XF, "$x_F$"), (XD, "$x_D$")]:
        ax.axvline(xv, ls=":", color=C["spec"], lw=1.0, zorder=1)
        ax.text(xv + 0.008, 0.03, lbl, fontsize=11,
                color=C["spec"], va="bottom")


def _draw_steps(ax: plt.Axes, segs: Segments,
                number_stages: bool = True) -> None:
    """
    Draw McCabe-Thiele step segments and optionally label each stage.

    Stage numbers are placed just above the point where the horizontal
    equilibrium step meets the equilibrium curve (end of each horiz. seg).
    """
    for i, (x1, y1, x2, y2) in enumerate(segs):
        lbl = "Theoretical stages" if i == 0 else "_nolegend_"
        ax.plot([x1, x2], [y1, y2], color=C["stages"], lw=1.8,
                label=lbl, zorder=4)

        # Horizontal segments are even-indexed; label at equilibrium-curve end
        if number_stages and i % 2 == 0:
            stage_num = i // 2 + 1
            ax.text(x2 + 0.010, y2 + 0.012, str(stage_num),
                    fontsize=16, color=C["stages"],
                    va="bottom", ha="left")


def _save_figure(fig: plt.Figure, stem: str) -> None:
    """Save figure as PNG at 300 dpi and as SVG to OUT_DIR."""
    for ext, kw in [(".png", {"dpi": 300})]:  # , (".svg", {})]:
        path = OUT_DIR / (stem + ext)
        fig.savefig(path, bbox_inches="tight", **kw)
        print(f"  Saved -> {path.name}")


# ─────────────────────────────────────────────────────────────────────────────
# 6.  FIGURE A — MINIMUM STAGES (TOTAL REFLUX)
# ─────────────────────────────────────────────────────────────────────────────

def plot_figure_a(segs: Segments, n_min: int) -> None:
    """
    Figure A: minimum number of theoretical stages at total reflux.
    Stages are stepped between the equilibrium curve and the 45° diagonal.
    """
    fig, ax = _make_axes(
        f"Figure A — Minimum Stages at Total Reflux\n"
        f"$N_{{\\mathrm{{min}}}}$ = {n_min}"
    )
    _draw_spec_verticals(ax)
    _draw_steps(ax, segs, number_stages=True)

    ax.legend(loc="upper left", fontsize=11, framealpha=0.9)
    fig.tight_layout()
    _save_figure(fig, "FigA_MinimumStages")
#    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# 7.  FIGURE B — MINIMUM REFLUX CONDITION
# ─────────────────────────────────────────────────────────────────────────────

def plot_figure_b() -> None:
    """
    Figure B: minimum reflux condition.
    Shows the q-line, the rectifying line at Rmin, and the pinch point.
    The pinch point is where the q-line intersects the equilibrium curve;
    passing the rectifying line through this point gives Rmin.
    """
    def rect_min_line(x: np.ndarray) -> np.ndarray:
        return _M_RECT_MIN * (x - XD) + XD

    # Clip q-line and Rmin rectifying line to the visible axis range [0,1]²
    x_all = np.linspace(0.0, 1.0, 500)

    y_q = q_line(x_all)
    mask_q = (y_q >= 0.0) & (y_q <= 1.0)

    y_rm = rect_min_line(x_all)
    mask_rm = (y_rm >= 0.0) & (y_rm <= 1.0) & (x_all <= XD)

    fig, ax = _make_axes(
        f"Figure B — Minimum Reflux Condition\n"
        f"$R_{{\\mathrm{{min}}}}$ = {RMIN:.3f}"
    )
    _draw_spec_verticals(ax)

    ax.plot(x_all[mask_q], y_q[mask_q], color=C["qline"], lw=2.0,
            label=f"$q$-line  ($q$ = {Q:.2f})", zorder=3)

    ax.plot(x_all[mask_rm], y_rm[mask_rm], color=C["rect"], lw=2.0, ls="--",
            label=f"Rectifying line at $R_{{\\mathrm{{min}}}}$", zorder=3)

    ax.plot(X_PINCH, Y_PINCH, "o", color=C["pinch"], ms=9, zorder=6,
            label=f"Pinch point  ({X_PINCH:.3f}, {Y_PINCH:.3f})")

    ax.legend(loc="upper left", fontsize=11, framealpha=0.9)
    fig.tight_layout()
    _save_figure(fig, "FigB_MinimumReflux")
#    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# 8.  FIGURE C — BUILDING THE McCABE-THIELE DESIGN FRAMEWORK
# ─────────────────────────────────────────────────────────────────────────────

def _label_framework_curves(ax: plt.Axes) -> None:
    """
    Place direct educational labels on each design element in Figure C.
    Arrow annotations are used so students can immediately connect each
    label to the corresponding curve without a separate legend.
    """
    _lkw = dict(
        fontsize=14, fontweight="bold", ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.88),
    )
    _akw = dict(arrowstyle="->", lw=1.3, connectionstyle="arc3,rad=0.15")

    # Equilibrium curve
    ax.annotate(
        "Equilibrium\nCurve",
        xy=(0.28, float(y_eq(0.28))),
        xytext=(0.09, 0.61),
        color=C["eq"],
        arrowprops=dict(**_akw, color=C["eq"]),
        **_lkw,
    )

    # 45° diagonal — small rotated label sitting on the line
    ax.text(
        0.83, 0.80, "45° Line",
        fontsize=14, color=C["diag"], ha="center",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.78),
    )

    # q-line
    ax.annotate(
        "$q$-Line",
        xy=(0.50, float(q_line(0.50))),
        xytext=(0.33, 0.89),
        color=C["qline"],
        arrowprops=dict(**_akw, color=C["qline"]),
        **_lkw,
    )

    # Rectifying operating line
    ax.annotate(
        "Rectifying\nLine",
        xy=(0.72, float(rectifying_line(0.72))),
        xytext=(0.60, 0.94),
        color=C["rect"],
        arrowprops=dict(**_akw, color=C["rect"]),
        **_lkw,
    )

    # Stripping operating line
    ax.annotate(
        "Stripping\nLine",
        xy=(0.22, float(stripping_line(0.22))),
        xytext=(0.06, 0.40),
        color=C["strip"],
        arrowprops=dict(**_akw, color=C["strip"]),
        **_lkw,
    )

    # Feed intersection point
    ax.annotate(
        "Feed\nIntersection",
        xy=(X_INT, Y_INT),
        xytext=(X_INT + 0.16, Y_INT - 0.13),
        color=C["qline"],
        arrowprops=dict(**_akw, color=C["qline"]),
        **_lkw,
    )


def plot_figure_c() -> None:
    """
    Figure C: Building the McCabe-Thiele Design Framework.

    Displays all curves needed before stage stepping can begin:
    equilibrium curve, 45° line, q-line, rectifying and stripping operating
    lines, and the feed intersection point. Stage stepping is absent by
    design — this slide answers "what information is required first?"

    Intended for a lecture slide: "Translating Distillation Theory
    into Design Decisions."
    """
    x_all  = np.linspace(0.0, 1.0, 500)
    y_q    = q_line(x_all)
    mask_q = (y_q >= 0.0) & (y_q <= 1.0)

    x_rect  = np.linspace(X_INT, XD, 200)
    x_strip = np.linspace(XB, X_INT, 200)

    fig, ax = _make_axes(
        "Figure C — Building the McCabe-Thiele Design Framework"
    )
    _draw_spec_verticals(ax)

    # q-line (clipped to axis range)
    ax.plot(x_all[mask_q], y_q[mask_q], color=C["qline"], lw=2.0, zorder=3)

    # Rectifying line: feed intersection → distillate
    ax.plot(x_rect, rectifying_line(x_rect), color=C["rect"], lw=2.0, zorder=3)

    # Stripping line: bottoms → feed intersection
    ax.plot(x_strip, stripping_line(x_strip), color=C["strip"], lw=2.0, zorder=3)

    # Feed intersection marker
    ax.plot(X_INT, Y_INT, "D", color=C["qline"], ms=8, zorder=6)

    # Direct curve labels — no legend, cleaner for a lecture slide
    _label_framework_curves(ax)

    fig.tight_layout()
    _save_figure(fig, "FigC_DesignFramework")
#    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# 9.  FIGURE D — ACTUAL OPERATING DESIGN
# ─────────────────────────────────────────────────────────────────────────────

def plot_figure_d(segs: Segments, n_stages: int, feed_stage: int) -> None:
    """
    Figure D: actual operating design.
    Shows the q-line, rectifying and stripping operating lines, theoretical
    stage stepping, feed stage annotation, and a results summary textbox.
    """
    x_all = np.linspace(0.0, 1.0, 500)
    y_q = q_line(x_all)
    mask_q = (y_q >= 0.0) & (y_q <= 1.0)

    x_rect  = np.linspace(X_INT, XD, 200)
    x_strip = np.linspace(XB, X_INT, 200)

    fig, ax = _make_axes(
        "Figure D — Actual Operating Design\n"
        f"$N$ = {n_stages} stages,   Feed stage = {feed_stage}"
    )
    _draw_spec_verticals(ax)

    # q-line (clipped to axis range)
    ax.plot(x_all[mask_q], y_q[mask_q], color=C["qline"], lw=1.8,
            label=f"$q$-line  ($q$ = {Q:.2f})", zorder=3)

    # Rectifying line: feed intersection → distillate
    ax.plot(x_rect, rectifying_line(x_rect), color=C["rect"], lw=2.0,
            label=f"Rectifying line  ($R$ = {R:.3f})", zorder=3)

    # Stripping line: bottoms → feed intersection
    ax.plot(x_strip, stripping_line(x_strip), color=C["strip"], lw=2.0,
            label="Stripping line", zorder=3)

    # Feed intersection marker
    ax.plot(X_INT, Y_INT, "D", color=C["qline"], ms=7, zorder=6,
            label=f"Feed intersection  ({X_INT:.3f}, {Y_INT:.3f})")

    # Stage stepping with numbered labels
    _draw_steps(ax, segs, number_stages=True)

    # Feed stage annotation arrow
    feed_h_idx = 2 * (feed_stage - 1)          # horizontal seg of feed stage
    if feed_h_idx < len(segs):
        x1, y1, x2, _ = segs[feed_h_idx]
        x_arrow = (x1 + x2) / 2.0
        y_arrow = y1
        ax.annotate(
            f"Feed stage {feed_stage}",
            xy=(x_arrow, y_arrow),
            xytext=(x_arrow + 0.12, y_arrow - 0.09),
            fontsize=11,
            color=C["stages"],
            arrowprops=dict(arrowstyle="->", color=C["stages"], lw=1.4),
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.85),
        )

    # Results textbox — monospace so columns align
    box_lines = [
        f"xF  =  {XF:.2f}",
        f"xD  =  {XD:.2f}",
        f"xB  =  {XB:.2f}",
        f"q   =  {Q:.2f}",
        "─" * 20,
        f"Rmin  =  {RMIN:.3f}",
        f"R     =  {R:.3f}",
        "─" * 20,
        f"Stages      =  {n_stages}",
        f"Feed stage  =  {feed_stage}",
    ]
    ax.text(
        0.020, 0.975, "\n".join(box_lines),
        transform=ax.transAxes,
        fontsize=14,
        va="top", ha="left",
        family="monospace",
        linespacing=1.55,
        bbox=dict(boxstyle="round,pad=0.55", fc="white",
                  ec="#94a3b8", alpha=0.93),
        zorder=10,
    )

    ax.legend(loc="lower right", fontsize=12, framealpha=0.9)
    fig.tight_layout()
    _save_figure(fig, "FigD_ActualDesign")
    #    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# 10. CSV EXPORT FOR POWERPOINT ANIMATION
# ─────────────────────────────────────────────────────────────────────────────

def export_csv(segs: Segments, path: Path) -> None:
    """
    Write stage step coordinates to a CSV file for PowerPoint animation.

    Each row is one line segment (horizontal or vertical).

    Columns
    -------
    Stage        : theoretical stage number (1-based)
    SegmentType  : 'horizontal' (equilibrium step) or 'vertical' (op-line step)
    x_start, y_start, x_end, y_end : segment endpoint coordinates
    """
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["Stage", "SegmentType", "x_start", "y_start", "x_end", "y_end"]
        )
        stage = 0
        for i, (x1, y1, x2, y2) in enumerate(segs):
            seg_type = "horizontal" if i % 2 == 0 else "vertical"
            if seg_type == "horizontal":
                stage += 1
            writer.writerow([
                stage, seg_type,
                f"{x1:.6f}", f"{y1:.6f}",
                f"{x2:.6f}", f"{y2:.6f}",
            ])
    print(f"  Saved  →  {path.name}")


# ─────────────────────────────────────────────────────────────────────────────
# 11. CONSOLE SUMMARY
# ─────────────────────────────────────────────────────────────────────────────

def print_summary(n_min: int, n_stages: int, feed_stage: int) -> None:
    """Print a formatted summary of all calculated design values."""
    sep = "─" * 48
    print(f"\n{sep}")
    print("  Example 7.1 — Benzene / Toluene Distillation")
    print(sep)
    print(f"  Feed composition           xF  =  {XF:.3f}")
    print(f"  Distillate specification   xD  =  {XD:.3f}")
    print(f"  Bottoms specification      xB  =  {XB:.3f}")
    print(f"  Feed quality                q  =  {Q:.3f}")
    print(sep)
    print(f"  Pinch point                 x  =  {X_PINCH:.4f}")
    print(f"                              y  =  {Y_PINCH:.4f}")
    print(f"  Minimum reflux           Rmin  =  {RMIN:.4f}")
    print(f"  Operating reflux            R  =  {R:.4f}   (= 1.3 × Rmin)")
    print(sep)
    print(f"  Feed intersection           x  =  {X_INT:.4f}")
    print(f"                              y  =  {Y_INT:.4f}")
    print(sep)
    print(f"  Minimum stages           Nmin  =  {n_min}")
    print(f"  Theoretical stages          N  =  {n_stages}")
    print(f"  Optimal feed stage             =  {feed_stage}")
    print(f"{sep}\n")


# ─────────────────────────────────────────────────────────────────────────────
# 12. MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    """Run all calculations, generate all figures, and export all files."""
    print("\nExample 7.1 — McCabe-Thiele Analysis")
    print("Benzene / Toluene at 1 atm\n")

    # Stage calculations
    segs_min, n_min             = step_minimum_stages()
    segs_act, n_stages, n_feed  = step_actual_stages()

    # Figures
    print("Generating figures …")
    plot_figure_a(segs_min, n_min)
    plot_figure_b()
    plot_figure_c()
    plot_figure_d(segs_act, n_stages, n_feed)

    # CSV
    export_csv(segs_act, OUT_DIR / "stage_steps.csv")

    # Summary
    print_summary(n_min, n_stages, n_feed)


if __name__ == "__main__":
    main()
