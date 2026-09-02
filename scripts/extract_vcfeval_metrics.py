#!/usr/bin/env python3
"""
extract_vcfeval_metrics.py
Parse rtg vcfeval summary.txt files for multiple individuals/regions and emit
a single CSV with TP, FP, FN, Precision, Recall, F1, plus runtime if a timing
file is present alongside the summary.

Usage:
    python3 extract_vcfeval_metrics.py <root_dir> [--out metrics.csv]

Expected directory layout (produced by run_5x5_snv_bench.sh):
    <root_dir>/
        HG001_r2/e_arcs_snv/summary.txt
        HG001_r2/e_disco_snv/summary.txt
        HG001_r2/e_k2s_snv/summary.txt
        HG001_r2/arcs_time.txt          (optional, from /usr/bin/time)
        HG001_r3/...
        HG002_na/...
        ...

summary.txt format (rtg vcfeval):
    Threshold  True-pos-baseline  True-pos-call  False-pos  False-neg  Precision  Sensitivity  F-measure
    None       <tp>               <tp_c>         <fp>       <fn>       <P>        <R>          <F1>
"""
import os, sys, csv, re, argparse

def parse_summary(path):
    """Return dict with TP, FP, FN, P, R, F1 from the 'None' threshold row."""
    with open(path) as f:
        lines = [l.rstrip() for l in f if l.strip()]
    for line in lines:
        if line.strip().startswith("None"):
            parts = line.split()
            # cols: None tp_base tp_call fp fn precision recall f1
            if len(parts) >= 8:
                return {
                    "TP": int(parts[1]),
                    "FP": int(parts[3]),
                    "FN": int(parts[4]),
                    "Precision": float(parts[5]),
                    "Recall":    float(parts[6]),
                    "F1":        float(parts[7]),
                }
    return None

def parse_time(path):
    """Extract wall-clock seconds from /usr/bin/time -v or bash time output."""
    if not os.path.exists(path):
        return ""
    with open(path) as f:
        txt = f.read()
    # /usr/bin/time -v: "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.23"
    m = re.search(r'Elapsed.*?:\s+(?:(\d+):)?(\d+):(\d+(?:\.\d+)?)', txt)
    if m:
        h = int(m.group(1) or 0)
        mins = int(m.group(2))
        secs = float(m.group(3))
        return f"{h*3600 + mins*60 + secs:.2f}"
    # bash time: "real\t0m1.230s"
    m = re.search(r'real\s+(\d+)m([\d.]+)s', txt)
    if m:
        return f"{int(m.group(1))*60 + float(m.group(2)):.2f}"
    return ""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="Root directory containing per-region subdirs")
    ap.add_argument("--out", default="metrics.csv")
    ap.add_argument("--variant", default="snv", choices=["snv","indel"],
                    help="Which variant type summary to collect")
    args = ap.parse_args()

    rows = []
    root = args.root
    vt = args.variant

    for name in sorted(os.listdir(root)):
        d = os.path.join(root, name)
        if not os.path.isdir(d):
            continue
        # parse individual/region from directory name e.g. HG001_r3, HG002_na
        parts = name.split("_", 1)
        ind = parts[0]
        region = parts[1] if len(parts) > 1 else "?"

        for tool, edir in [("ARCS", f"e_arcs_{vt}"),
                            ("DiscoSNP++", f"e_disco_{vt}"),
                            ("Kmer2SNP",   f"e_k2s_{vt}")]:
            sfile = os.path.join(d, edir, "summary.txt")
            if not os.path.exists(sfile):
                continue
            m = parse_summary(sfile)
            if m is None:
                print(f"WARN: could not parse {sfile}", file=sys.stderr)
                continue
            tfile = os.path.join(d, f"{tool.lower().replace('+','').replace(' ','_')}_time.txt")
            rows.append({
                "Individual": ind,
                "Region":     region,
                "Tool":       tool,
                "TP":         m["TP"],
                "FP":         m["FP"],
                "FN":         m["FN"],
                "Precision":  f"{m['Precision']:.3f}",
                "Recall":     f"{m['Recall']:.3f}",
                "F1":         f"{m['F1']:.3f}",
                "Runtime_s":  parse_time(tfile),
            })

    if not rows:
        print("No summary.txt files found under", root, file=sys.stderr)
        sys.exit(1)

    fields = ["Individual","Region","Tool","TP","FP","FN","Precision","Recall","F1","Runtime_s"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.out}")

    # Also print a quick console summary grouped by individual
    print("\n--- Quick summary (F1) ---")
    print(f"{'Individual':<12} {'Region':<8} {'ARCS':>6} {'DiscoSNP++':>11} {'Kmer2SNP':>9}")
    by_ind_reg = {}
    for r in rows:
        k = (r["Individual"], r["Region"])
        by_ind_reg.setdefault(k, {})[r["Tool"]] = r["F1"]
    for k in sorted(by_ind_reg):
        d2 = by_ind_reg[k]
        print(f"{k[0]:<12} {k[1]:<8} {d2.get('ARCS','--'):>6} "
              f"{d2.get('DiscoSNP++','--'):>11} {d2.get('Kmer2SNP','--'):>9}")

if __name__ == "__main__":
    main()
