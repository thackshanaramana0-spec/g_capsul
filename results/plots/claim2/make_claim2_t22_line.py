"""
Journal-ready line chart: precision / recall / F1 vs. sequencing depth for
CAPSULE's het-SNV calls on HG002, from benchmark/results/claim2_T2.2_coverage_sweep.csv.

Single individual, four depth points (10x/15x/20x/30x), CAPSULE only -- no
DiscoSNP++/Kmer2SNP comparator at each depth in this CSV, so this is an
illustrative sweep of CAPSULE's own behavior, not a head-to-head comparison.
Sweep form matches DiscoSNP++'s own Figure 5 (precision/recall vs. a swept
variable) and Kmer2SNP's Figure 4 (metrics vs. coverage).

Same PLOS ONE styling as the other Claim 1/2 charts: Arial 9pt bold, 0.2mm
line weight, RGB/LZW TIFF capped at 2250px width / 300dpi.

Usage:
    python make_claim2_t22_line.py
Outputs (next to this script):
    claim2_t22_coverage_sweep.pdf
    claim2_t22_coverage_sweep.png   -- 600 dpi preview
    claim2_t22_coverage_sweep.tif   -- PLOS ONE-compliant submission file
"""
import csv
import pathlib
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image

HERE = pathlib.Path(__file__).parent
CSV_PATH = HERE.parent.parent / "claim2" / "claim2_T2.2_coverage_sweep.csv"

# ---- palette (dataviz skill reference palette, categorical slots 1-3) ----
BLUE = "#2a78d6"    # precision
ORANGE = "#eb6834"  # recall
AQUA = "#1baf7a"    # F1
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"

METRICS = ["precision", "recall", "f1"]
COLORS = {"precision": BLUE, "recall": ORANGE, "f1": AQUA}
LABELS = {"precision": "Precision", "recall": "Recall", "f1": "F1"}

# ---- load ----
points = []
with open(CSV_PATH, newline="") as f:
    for r in csv.DictReader(f):
        if r["individual"] == "HG002":
            points.append((int(r["depth_x"]), r["precision"], r["recall"], r["f1"]))
points.sort(key=lambda p: p[0])
depths = [p[0] for p in points]
series = {
    "precision": [float(p[1]) for p in points],
    "recall": [float(p[2]) for p in points],
    "f1": [float(p[3]) for p in points],
}

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
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

fig_h_in = 90 / 25.4
fig_w_in = 6.0
fig, ax = plt.subplots(figsize=(fig_w_in, fig_h_in))

for metric in METRICS:
    ax.plot(
        depths, series[metric], color=COLORS[metric], label=LABELS[metric],
        linewidth=2.0, marker="o", markersize=5, markeredgecolor="white",
        markeredgewidth=0.6, zorder=3,
    )

ax.set_xticks(depths)
ax.set_xticklabels([f"{d}×" for d in depths], fontsize=9, fontweight="bold")
ax.set_ylim(0, 1.0)
ax.set_xlabel("Sequencing depth", fontsize=9, fontweight="bold", color=INK, labelpad=8)
ax.set_ylabel("Score (HG002, het-SNV)", fontsize=9, fontweight="bold", color=INK)
ax.tick_params(axis="y", colors=INK, labelsize=9)
ax.tick_params(axis="x", colors=INK, length=0)
for label in ax.get_yticklabels():
    label.set_fontweight("bold")

for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.spines["left"].set_color(BASELINE)
ax.spines["bottom"].set_color(BASELINE)

ax.yaxis.grid(True, color=GRID, linewidth=0.5, zorder=0)
ax.set_axisbelow(True)

legend = ax.legend(
    loc="lower center", bbox_to_anchor=(0.5, 1.12), ncol=3,
    frameon=False, fontsize=9,
    handlelength=1.4, columnspacing=1.2,
    prop={"weight": "bold", "size": 9},
)
for text in legend.get_texts():
    text.set_color(INK)

ax.set_title(
    "CAPSULE het-SNV accuracy vs. sequencing depth",
    fontsize=10.5, fontweight="bold", color=INK, loc="center", pad=42,
)

fig.tight_layout()
fig.subplots_adjust(bottom=0.18, top=0.80)
fig.savefig(HERE / "claim2_t22_coverage_sweep.pdf")
fig.savefig(HERE / "claim2_t22_coverage_sweep.png", dpi=600)
print("wrote", HERE / "claim2_t22_coverage_sweep.pdf")
print("wrote", HERE / "claim2_t22_coverage_sweep.png")

# ---- PLOS ONE-compliant TIFF: 300 dpi, RGB, flattened, LZW-compressed ----
tiff_dpi = 300
w_px = round(fig_w_in * tiff_dpi)
h_px = round(fig_h_in * tiff_dpi)
assert 789 <= w_px <= 2250, f"width {w_px}px outside PLOS 789-2250px range"
assert h_px <= 2625, f"height {h_px}px exceeds PLOS 2625px max"

png_buf = HERE / "_claim2_t22_tmp.png"
fig.savefig(png_buf, dpi=tiff_dpi, facecolor="white")
img = Image.open(png_buf).convert("RGB")
tiff_path = HERE / "claim2_t22_coverage_sweep.tif"
img.save(tiff_path, format="TIFF", compression="tiff_lzw", dpi=(tiff_dpi, tiff_dpi))
png_buf.unlink()

size_mb = tiff_path.stat().st_size / (1024 * 1024)
assert size_mb < 10, f"TIFF is {size_mb:.1f} MB, exceeds PLOS 10 MB limit"
print(f"wrote {tiff_path}  ({w_px}x{h_px}px @ {tiff_dpi}dpi, {size_mb:.2f} MB, RGB, LZW)")
