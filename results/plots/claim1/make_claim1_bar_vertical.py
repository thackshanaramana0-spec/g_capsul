"""
Journal-ready bar graph: archive size as % of raw input (ratio_pct), one bar
per tool (CAPSULE, SPRING, Genozip) per dataset, from
benchmark/results/claim1_T1.1_T1.2.csv.

Datasets on the X axis, bars projecting upward (shorter bar = smaller
archive). Vertical grouped bars, sorted by dataset size. Interpretive notes
(sample size, "lossless", direction of a win) belong in the figure caption,
not on the plot -- kept off this chart per standard journal figure style.

PLOS ONE figure rules applied here exactly (https://journals.plos.org/plosone/s/figures):
  - file format: TIFF, LZW-compressed, flattened (no alpha), RGB
  - resolution: 300 dpi
  - width: 789-2250 px at 300 dpi (this figure is capped at 2250 px / 7.5 in)
  - height: <=2625 px at 300 dpi
  - font: Arial only, 8-12 pt (used 9 pt throughout)
  - line weight: 0.2 mm (~0.57 pt)
  - file size: <10 MB

Usage:
    python make_claim1_bar_vertical.py
Outputs (next to this script):
    claim1_margin_vertical.pdf   -- vector, for Nature/Science-style submission
    claim1_margin_vertical.png   -- 600 dpi raster, for quick preview
    claim1_margin_vertical.tif   -- PLOS ONE-compliant submission file
"""
import csv
import pathlib
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image

HERE = pathlib.Path(__file__).parent
CSV_PATH = HERE.parent.parent / "claim1" / "claim1_T1.1_T1.2.csv"

# ---- palette (dataviz skill reference palette, categorical slots 1-3) ----
BLUE = "#2a78d6"    # Genozip
ORANGE = "#eb6834"  # SPRING
AQUA = "#1baf7a"    # CAPSULE
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"

TOOLS = ["Genozip", "SPRING", "CAPSULE"]
COLORS = {"CAPSULE": AQUA, "SPRING": ORANGE, "Genozip": BLUE}

# ---- load ----
rows = {}
with open(CSV_PATH, newline="") as f:
    for r in csv.DictReader(f):
        rows.setdefault(r["dataset"], {})[r["tool"]] = float(r["ratio_pct"])

datasets = sorted(rows, key=lambda d: (not d.startswith("HG"), -rows[d]["CAPSULE"]))
values = {tool: [rows[d][tool] for d in datasets] for tool in TOOLS}

# ---- style: Arial, bold, 8-12pt (PLOS: "Arial, Times, or Symbol ... 8-12 point") ----
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.weight": "bold",
    "font.size": 9,
    "axes.linewidth": 0.57,   # PLOS: 0.2 mm line weight
    "axes.labelweight": "bold",
    "xtick.major.width": 0.57,
    "ytick.major.width": 0.57,
    "pdf.fonttype": 42,   # embed fonts as editable text, not curves
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

n = len(datasets)
bar_w = 0.22
x = range(n)
offsets = [-bar_w, 0, bar_w]  # CAPSULE, SPRING, Genozip

# Width capped so 300 dpi output never exceeds PLOS's 2250 px max width.
fig_h_in = 90 / 25.4
fig_w_in = min(7.5, max(7.0, 0.37 * n + 0.6))
fig, ax = plt.subplots(figsize=(fig_w_in, fig_h_in))

for tool, off in zip(TOOLS, offsets):
    ax.bar(
        [i + off for i in x], values[tool], width=bar_w,
        color=COLORS[tool], label=tool, zorder=3,
    )

ax.axhline(0, color=BASELINE, linewidth=0.8, zorder=2)
ax.set_xticks(list(x))
ax.set_xticklabels(datasets, fontsize=9, fontweight="bold", rotation=45, ha="right", rotation_mode="anchor")
ax.set_ylim(bottom=0)
ax.set_xlabel("Dataset", fontsize=9, fontweight="bold", color=INK, labelpad=8)

ax.set_ylabel(
    "Archive size as % of raw input",
    fontsize=9, fontweight="bold", color=INK,
)
ax.tick_params(axis="y", colors=INK, labelsize=9)
ax.tick_params(axis="x", colors=INK, length=0)
for label in ax.get_yticklabels():
    label.set_fontweight("bold")

for spine in ("top", "right", "left"):
    ax.spines[spine].set_visible(False)
ax.spines["bottom"].set_color(BASELINE)

ax.yaxis.grid(True, color=GRID, linewidth=0.5, zorder=0)
ax.set_axisbelow(True)

legend = ax.legend(
    loc="upper right", frameon=False, fontsize=9,
    handlelength=1.0, handleheight=1.0, borderaxespad=0.4,
    prop={"weight": "bold", "size": 9},
)
for text in legend.get_texts():
    text.set_color(INK)

ax.set_title(
    "Archive size vs raw input: CAPSULE, SPRING, Genozip",
    fontsize=10.5, fontweight="bold", color=INK, loc="center", pad=10,
)

fig.tight_layout()
fig.subplots_adjust(bottom=0.36)
fig.savefig(HERE / "claim1_margin_vertical.pdf")
fig.savefig(HERE / "claim1_margin_vertical.png", dpi=600)
print("wrote", HERE / "claim1_margin_vertical.pdf")
print("wrote", HERE / "claim1_margin_vertical.png")

# ---- PLOS ONE-compliant TIFF: 300 dpi, RGB, flattened, LZW-compressed ----
tiff_dpi = 300
w_px = round(fig_w_in * tiff_dpi)
h_px = round(fig_h_in * tiff_dpi)
assert 789 <= w_px <= 2250, f"width {w_px}px outside PLOS 789-2250px range"
assert h_px <= 2625, f"height {h_px}px exceeds PLOS 2625px max"

png_buf = HERE / "_claim1_margin_vertical_tmp.png"
fig.savefig(png_buf, dpi=tiff_dpi, facecolor="white")
img = Image.open(png_buf).convert("RGB")  # flatten alpha -> RGB, no channel layers
tiff_path = HERE / "claim1_margin_vertical.tif"
img.save(tiff_path, format="TIFF", compression="tiff_lzw", dpi=(tiff_dpi, tiff_dpi))
png_buf.unlink()

size_mb = tiff_path.stat().st_size / (1024 * 1024)
assert size_mb < 10, f"TIFF is {size_mb:.1f} MB, exceeds PLOS 10 MB limit"
print(f"wrote {tiff_path}  ({w_px}x{h_px}px @ {tiff_dpi}dpi, {size_mb:.2f} MB, RGB, LZW)")
