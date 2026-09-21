"""
Journal-ready bar graph: het-SNV F1 score, one bar per tool (Kmer2SNP,
DiscoSNP++, CAPSULE) per individual, from benchmark/results/claim2_T2.1_snv.csv.

Same vertical-grouped-bar form and PLOS ONE styling as
benchmark/results/plots/make_claim1_bar_vertical.py: Arial 9pt bold, 0.2mm
line weight, RGB/LZW TIFF capped at 2250px width / 300dpi.

Usage:
    python make_claim2_t21_bar.py
Outputs (next to this script):
    claim2_t21_snv_f1.pdf
    claim2_t21_snv_f1.png   -- 600 dpi preview
    claim2_t21_snv_f1.tif   -- PLOS ONE-compliant submission file
"""
import csv
import pathlib
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image

HERE = pathlib.Path(__file__).parent
CSV_PATH = HERE.parent.parent / "claim2" / "claim2_T2.1_snv.csv"

# ---- palette (dataviz skill reference palette, categorical slots 1-3) ----
BLUE = "#2a78d6"    # Kmer2SNP (worst)
ORANGE = "#eb6834"  # DiscoSNP++ (middle)
AQUA = "#1baf7a"    # CAPSULE (best)
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"

TOOLS = ["Kmer2SNP", "DiscoSNP++", "CAPSULE"]
COLORS = {"CAPSULE": AQUA, "DiscoSNP++": ORANGE, "Kmer2SNP": BLUE}

# ---- load ----
rows = {}
with open(CSV_PATH, newline="") as f:
    for r in csv.DictReader(f):
        rows.setdefault(r["individual"], {})[r["tool"]] = float(r["f1"])

individuals = sorted(rows)  # HG002, HG003, HG004, HG005
values = {tool: [rows[d][tool] for d in individuals] for tool in TOOLS}

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

n = len(individuals)
bar_w = 0.22
x = range(n)
offsets = [-bar_w, 0, bar_w]  # Kmer2SNP, DiscoSNP++, CAPSULE

fig_h_in = 90 / 25.4
fig_w_in = min(7.5, max(4.0, 0.9 * n + 1.5))
fig, ax = plt.subplots(figsize=(fig_w_in, fig_h_in))

for tool, off in zip(TOOLS, offsets):
    ax.bar(
        [i + off for i in x], values[tool], width=bar_w,
        color=COLORS[tool], label=tool, zorder=3,
    )

ax.axhline(0, color=BASELINE, linewidth=0.8, zorder=2)
ax.set_xticks(list(x))
ax.set_xticklabels(individuals, fontsize=9, fontweight="bold")
ax.set_ylim(0, 1.0)
ax.set_xlabel("Individual", fontsize=9, fontweight="bold", color=INK, labelpad=8)
ax.set_ylabel("het-SNV F1 score", fontsize=9, fontweight="bold", color=INK)
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
    loc="lower center", bbox_to_anchor=(0.5, 1.12), ncol=3,
    frameon=False, fontsize=9,
    handlelength=1.0, handleheight=1.0, columnspacing=1.2,
    prop={"weight": "bold", "size": 9},
)
for text in legend.get_texts():
    text.set_color(INK)

ax.set_title(
    "het-SNV calling accuracy: CAPSULE, DiscoSNP++, Kmer2SNP",
    fontsize=10.5, fontweight="bold", color=INK, loc="center", pad=42,
)

fig.tight_layout()
fig.subplots_adjust(bottom=0.18, top=0.80)
fig.savefig(HERE / "claim2_t21_snv_f1.pdf")
fig.savefig(HERE / "claim2_t21_snv_f1.png", dpi=600)
print("wrote", HERE / "claim2_t21_snv_f1.pdf")
print("wrote", HERE / "claim2_t21_snv_f1.png")

# ---- PLOS ONE-compliant TIFF: 300 dpi, RGB, flattened, LZW-compressed ----
tiff_dpi = 300
w_px = round(fig_w_in * tiff_dpi)
h_px = round(fig_h_in * tiff_dpi)
assert 789 <= w_px <= 2250, f"width {w_px}px outside PLOS 789-2250px range"
assert h_px <= 2625, f"height {h_px}px exceeds PLOS 2625px max"

png_buf = HERE / "_claim2_t21_tmp.png"
fig.savefig(png_buf, dpi=tiff_dpi, facecolor="white")
img = Image.open(png_buf).convert("RGB")
tiff_path = HERE / "claim2_t21_snv_f1.tif"
img.save(tiff_path, format="TIFF", compression="tiff_lzw", dpi=(tiff_dpi, tiff_dpi))
png_buf.unlink()

size_mb = tiff_path.stat().st_size / (1024 * 1024)
assert size_mb < 10, f"TIFF is {size_mb:.1f} MB, exceeds PLOS 10 MB limit"
print(f"wrote {tiff_path}  ({w_px}x{h_px}px @ {tiff_dpi}dpi, {size_mb:.2f} MB, RGB, LZW)")
