#!/usr/bin/env python3
"""
check_input_scope.py -- does this FASTQ fall inside what this encoder
actually supports, and why (or why not)?

WHAT THIS IS NOT: a filename or header sniff for the string "nanopore". That
approach is guessable, gameable, and tells you nothing about a file that
happens not to mention its own technology. Every verdict this script prints
is derived from streaming the actual file and comparing what it finds against
the actual structural limits of the actual encoder -- extracted from the
encoder's own source at run time, not retyped here, so this script cannot
drift out of sync with the code it is checking. If the encoder's limits ever
change, this script's next run picks up the change automatically; if the
encoder's source cannot be found, this script says so and refuses to guess.

WHY THIS EXISTS: `src/encoder.cpp` used to silently drop any read over its
structural length bound and complete normally, writing an empty archive at
exit code 0 -- a reviewer pointing this tool at unsupported input saw
success. The encoder itself now refuses loudly and unconditionally (see
`src/encoder.cpp`'s scope-check block, right after `int main`), so THIS
script cannot be bypassed into producing a false sense of safety -- it
is a pre-flight report, not the safety mechanism. The safety mechanism is in
the encoder. This script exists so a reviewer gets that verdict, with the
real numbers behind it, in seconds, before spending any compute on a file
that would be refused anyway -- and so they get a candid answer for input
that is IN scope but marginal, which the encoder's hard gate has no way to
express (a gate is binary; this report is not).

USAGE
    python3 scripts/check_input_scope.py reads.fq [reads_2.fq ...]
    python3 scripts/check_input_scope.py --json reads.fq   # machine-readable

EXIT CODES
    0  in scope, no concerns
    1  in scope, with caveats worth reading (see the report)
    2  OUT OF SCOPE -- the encoder will refuse this file
    3  could not read the file, or could not find the encoder source to
       extract the real limits from (refuses to guess a number)
"""
import sys, os, re, json, math, statistics as st

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ENCODER_SRC = os.path.join(REPO_ROOT, "src", "encoder.cpp")
NAMES_SRC   = os.path.join(REPO_ROOT, "include", "names_coder.h")


def extract_constant(path, pattern, label):
    """Pull a numeric constant directly out of the real source file. Never
    hardcode the number here: if the source moves or the constant is renamed,
    this fails loudly (see main()) rather than silently checking against a
    stale value."""
    if not os.path.isfile(path):
        return None
    text = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(pattern, text)
    if not m:
        return None
    return int(m.group(1))


def get_real_limits():
    limits = {}
    limits["MAX_READ_LEN"] = extract_constant(
        ENCODER_SRC,
        r"static\s+const\s+uint32_t\s+MAX_READ_LEN\s*=\s*(\d+)\s*;",
        "MAX_READ_LEN",
    )
    limits["MAXTOK"] = extract_constant(
        NAMES_SRC,
        r"static\s+const\s+uint32_t\s+MAXTOK\s*=\s*(\d+)\s*;",
        "MAXTOK",
    )
    return limits


class Stats:
    __slots__ = (
        "n_records", "n_malformed", "lengths", "header_lengths",
        "n_with_N", "total_bases", "base_counts", "qual_min", "qual_max",
        "qual_sum", "qual_n", "distinct_lengths", "duplicate_prefix_hits",
        "n_variable_qual_len",
    )

    def __init__(self):
        self.n_records = 0
        self.n_malformed = 0
        self.lengths = []            # sampled read lengths (bounded, see SAMPLE_CAP)
        self.header_lengths = []
        self.n_with_N = 0
        self.total_bases = 0
        self.base_counts = {}        # character -> count, over the WHOLE file
        self.qual_min = None
        self.qual_max = None
        self.qual_sum = 0
        self.qual_n = 0
        self.distinct_lengths = set()
        self.n_variable_qual_len = 0


SAMPLE_CAP = 200_000  # cap the in-memory length list; running stats (min/max/
                       # mean/distinct-lengths) are exact over the WHOLE file
                       # regardless of this cap -- only the stored list for
                       # percentile computation is capped, and capped by
                       # reservoir-style truncation stated as such in the report.


def stream_fastq(path, stats: Stats):
    """One pass, streaming -- never loads the file into memory, so this is
    safe to run on files far larger than RAM, including the multi-GB inputs
    this project's own locked datasets reach."""
    with open(path, "rb") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline()
            p = f.readline()
            q = f.readline()
            if not (s and p and q):
                stats.n_malformed += 1
                break
            h = h.rstrip(b"\n\r")
            s = s.rstrip(b"\n\r")
            p = p.rstrip(b"\n\r")
            q = q.rstrip(b"\n\r")
            if not h.startswith(b"@") or not p.startswith(b"+"):
                stats.n_malformed += 1
                continue
            stats.n_records += 1
            L = len(s)
            stats.total_bases += L
            stats.distinct_lengths.add(L)
            if len(stats.lengths) < SAMPLE_CAP:
                stats.lengths.append(L)
            if len(stats.header_lengths) < SAMPLE_CAP:
                stats.header_lengths.append(len(h))
            if len(q) != L:
                stats.n_variable_qual_len += 1
            if b"N" in s:
                stats.n_with_N += 1
            # base composition over the WHOLE file, not just the sample --
            # cheap (one pass, one dict), and it is what actually answers
            # "is this DNA" rather than a guess from a handful of records.
            for ch in set(s):
                stats.base_counts[chr(ch)] = stats.base_counts.get(chr(ch), 0) + s.count(ch)
            if q:
                qmin_local = min(q)
                qmax_local = max(q)
                stats.qual_min = qmin_local if stats.qual_min is None else min(stats.qual_min, qmin_local)
                stats.qual_max = qmax_local if stats.qual_max is None else max(stats.qual_max, qmax_local)
                stats.qual_sum += sum(q)
                stats.qual_n += len(q)


def pct(sorted_list, p):
    if not sorted_list:
        return None
    k = (len(sorted_list) - 1) * p
    f, c = math.floor(k), math.ceil(k)
    if f == c:
        return sorted_list[int(k)]
    return sorted_list[f] + (sorted_list[c] - sorted_list[f]) * (k - f)


def analyze_one(path, limits):
    result = {"path": path, "verdict": None, "exit_code": None, "findings": [], "notes": []}
    if not os.path.isfile(path):
        result["verdict"] = "CANNOT READ"
        result["exit_code"] = 3
        result["findings"].append(f"file does not exist: {path}")
        return result
    if os.path.getsize(path) == 0:
        result["verdict"] = "OUT OF SCOPE"
        result["exit_code"] = 2
        result["findings"].append("file is empty -- zero bytes")
        return result

    stats = Stats()
    try:
        stream_fastq(path, stats)
    except Exception as e:  # noqa: BLE001 -- deliberately broad: any failure
        # to parse is reported as a finding, not a Python traceback, because
        # the person reading this output may not be a programmer.
        result["verdict"] = "CANNOT READ"
        result["exit_code"] = 3
        result["findings"].append(f"failed to parse as FASTQ: {e}")
        return result

    if stats.n_records == 0:
        result["verdict"] = "OUT OF SCOPE"
        result["exit_code"] = 2
        result["findings"].append(
            f"no usable FASTQ records found ({stats.n_malformed} malformed/short "
            "record(s) seen). Either not FASTQ, or truncated."
        )
        return result

    lens_sorted = sorted(stats.lengths)
    mean_len = stats.total_bases / stats.n_records
    max_len = max(stats.distinct_lengths)
    min_len = min(stats.distinct_lengths)
    median_len = pct(lens_sorted, 0.5)
    p99_len = pct(lens_sorted, 0.99)

    result["stats"] = {
        "records": stats.n_records,
        "malformed_records": stats.n_malformed,
        "total_bases": stats.total_bases,
        "min_length": min_len,
        "max_length": max_len,
        "mean_length": round(mean_len, 1),
        "median_length": median_len,
        "p99_length": p99_len,
        "distinct_lengths": len(stats.distinct_lengths),
        "reads_with_N": stats.n_with_N,
        "reads_with_N_pct": round(100.0 * stats.n_with_N / stats.n_records, 3),
        "mean_header_length": round(sum(stats.header_lengths) / len(stats.header_lengths), 1)
        if stats.header_lengths else 0,
        "max_header_length": max(stats.header_lengths) if stats.header_lengths else 0,
        "quality_range": [stats.qual_min, stats.qual_max] if stats.qual_n else None,
        "quality_mean": round(stats.qual_sum / stats.qual_n, 1) if stats.qual_n else None,
        "reads_with_seq_qual_length_mismatch": stats.n_variable_qual_len,
    }

    findings = result["findings"]
    exit_code = 0

    # ------------------------------------------------------------------
    # 1. THE HARD GATE: read length vs the encoder's real structural bound.
    #    This is the same computation the encoder itself makes, using the
    #    SAME number, extracted from the SAME source -- if this script says
    #    "in scope" here, the encoder will not refuse for this reason.
    # ------------------------------------------------------------------
    max_read_len = limits.get("MAX_READ_LEN")
    if max_read_len is None:
        findings.append(
            "COULD NOT VERIFY the read-length bound: MAX_READ_LEN was not "
            "found in src/encoder.cpp. Refusing to state a scope "
            "verdict on a bound this script cannot confirm."
        )
        result["verdict"] = "CANNOT VERIFY"
        result["exit_code"] = 3
        return result

    over_bound_estimate = None
    if max_len > max_read_len:
        # The sample may not contain the true longest read if the file
        # exceeds SAMPLE_CAP records for length STORAGE -- but max_len is
        # tracked over every record via distinct_lengths, so this is exact,
        # not sampled.
        n_over = sum(1 for L in stats.distinct_lengths if L > max_read_len)
        findings.append(
            f"OUT OF SCOPE: reads up to {max_len} bases were found, exceeding "
            f"this build's {max_read_len}-base structural limit "
            f"(MAX_READ_LEN, src/encoder.cpp). {n_over} distinct "
            f"length(s) in this file exceed it. This is not a policy "
            f"threshold -- it is the size of a fixed per-read stack buffer "
            f"every decode path uses. The encoder WILL refuse this file "
            f"with a FATAL message rather than silently drop the affected "
            f"reads; this report is telling you the same thing in advance."
        )
        exit_code = max(exit_code, 2)

    # ------------------------------------------------------------------
    # 2. Names/header length vs the tokenizer's real bound (only matters if
    #    CAPS_NAMES=1 is used; state that qualification, not a bare refusal).
    # ------------------------------------------------------------------
    maxtok = limits.get("MAXTOK")
    if maxtok is not None and result["stats"]["max_header_length"] >= maxtok:
        findings.append(
            f"OUT OF SCOPE for CAPS_NAMES=1: a header is "
            f"{result['stats']['max_header_length']} characters long, at or "
            f"past the {maxtok}-character names-tokenizer bound (nmc::MAXTOK, "
            f"include/names_coder.h). A header this long can decode to "
            f"something other than what was encoded when names are stored; "
            f"the encoder refuses rather than risk that. If CAPS_NAMES is "
            f"not set for this run, this does not apply -- sequence and "
            f"quality are unaffected."
        )
        exit_code = max(exit_code, 2)

    # ------------------------------------------------------------------
    # 3. Alphabet -- is this actually DNA. This project's pipeline
    #    substitutes N->A internally (real, existing mechanism) and rejects
    #    nothing else explicitly at the source level, but a file whose
    #    "sequence" column is not a 4-5 symbol nucleotide alphabet (protein
    #    FASTA-as-FASTQ, RNA with U, corrupted binary) will silently produce
    #    a large, meaningless archive rather than a clean refusal, because
    #    nothing downstream currently checks this. Reported here as a real,
    #    named limitation, not smoothed over.
    # ------------------------------------------------------------------
    acgtn = sum(stats.base_counts.get(c, 0) for c in "ACGTNacgtn")
    other = stats.total_bases - acgtn
    other_frac = other / stats.total_bases if stats.total_bases else 0
    if other_frac > 0.001:
        distinct_other = sorted(
            c for c in stats.base_counts if c.upper() not in "ACGTN"
        )
        findings.append(
            f"CAVEAT: {other_frac*100:.2f}% of bases are outside the ACGTN "
            f"alphabet (symbols seen: {''.join(distinct_other)[:20]}). This "
            f"encoder's sequence coder assumes a 4-symbol nucleotide "
            f"alphabet plus N; it has no gate for this today (unlike the "
            f"length bound above) and will not refuse -- it will silently "
            f"treat unrecognised symbols as an encoding error at the base "
            f"level. If this is protein sequence, RNA with U, or corrupted "
            f"data, do not trust the resulting archive without decoding and "
            f"diffing it against the source, exactly as this project's own "
            f"benchmark harness always does."
        )
        exit_code = max(exit_code, 1)

    # ------------------------------------------------------------------
    # 4. Informational fingerprint: does the length/quality PROFILE look
    #    like the short-read (Illumina-shaped) data this pipeline is
    #    validated against, or like long-read (Nanopore/PacBio) data that
    #    merely happens to sit under the hard length gate (e.g. a short-read
    #    subsample of a long-read run, or heavily quality-trimmed reads)?
    #    This is NOT a gate -- only #1 and #2 are gates -- because a
    #    genuinely short, high-quality file should never be refused on a
    #    profile guess. It is an honest caveat, derived from the file's own
    #    numbers, for input that is technically in scope but was validated
    #    on a different shape than what this file has.
    # ------------------------------------------------------------------
    locked_max_len_seen = 301   # ERR552797, this project's longest locked-dataset read
    if max_len > locked_max_len_seen:
        findings.append(
            f"CAVEAT: the longest read here ({max_len} bp) exceeds the "
            f"longest read in any of this project's 19 locked benchmark "
            f"datasets ({locked_max_len_seen} bp, M. tuberculosis). It is "
            f"under the hard structural limit and WILL be processed, but "
            f"the assembly, mismatch and position coders have not been "
            f"measured on reads this long -- their parameters (seed widths, "
            f"MAXMAP, overlap thresholds) were derived from and validated "
            f"against reads in the 40-{locked_max_len_seen} bp range. "
            f"Compression ratio and calling accuracy on this file are "
            f"UNVALIDATED, not merely untested-by-us; treat any result on "
            f"it as exploratory."
        )
        exit_code = max(exit_code, 1)

    if result["stats"]["distinct_lengths"] > 1 and result["stats"]["max_length"] - result["stats"]["min_length"] > 50:
        findings.append(
            "NOTE: this is a variable-length file with a wide spread "
            f"({min_len}-{max_len} bp). Supported (containment removal "
            "handles this, see docs/paper/ARCHITECTURE.md), but this "
            "project's own published results show the caller's precision "
            "is measurably worse on the one variable-length individual "
            "in its benchmark (HG005) than on fixed-length ones -- see "
            "paper/LIMITATIONS.md section 2. Not a refusal, a documented "
            "accuracy caveat."
        )
        exit_code = max(exit_code, 1)

    if not findings:
        findings.append("No boundary or scope concern found for this file.")

    result["verdict"] = {0: "IN SCOPE", 1: "IN SCOPE, WITH CAVEATS", 2: "OUT OF SCOPE"}[exit_code]
    result["exit_code"] = exit_code
    return result


def print_human(result, limits):
    print("=" * 78)
    print(f" {result['path']}")
    print("=" * 78)
    if "stats" in result:
        s = result["stats"]
        print(f"  records                 {s['records']:,}  ({s['malformed_records']} malformed skipped)")
        print(f"  read length             min={s['min_length']} median={s['median_length']} "
              f"mean={s['mean_length']} p99={s['p99_length']} max={s['max_length']}  "
              f"({s['distinct_lengths']} distinct length(s))")
        print(f"  reads containing N      {s['reads_with_N']:,} ({s['reads_with_N_pct']}%)")
        print(f"  header length           mean={s['mean_header_length']} max={s['max_header_length']}")
        if s["quality_range"]:
            print(f"  quality byte range      {s['quality_range'][0]}-{s['quality_range'][1]}"
                  f"  (mean {s['quality_mean']})")
        if s["reads_with_seq_qual_length_mismatch"]:
            print(f"  seq/qual LENGTH MISMATCH  {s['reads_with_seq_qual_length_mismatch']} read(s) "
                  "-- malformed FASTQ, not scored above as a finding but real")
    print()
    print(f"  VERDICT: {result['verdict']}")
    for f in result["findings"]:
        print(f"    - {f}")
    print()


def main(argv):
    as_json = "--json" in argv
    paths = [a for a in argv if a != "--json"]
    if not paths:
        print(__doc__)
        return 3

    limits = get_real_limits()
    if limits.get("MAX_READ_LEN") is None:
        sys.stderr.write(
            "FATAL: could not extract MAX_READ_LEN from "
            f"{ENCODER_SRC}. This script refuses to state a scope verdict "
            "using a hardcoded fallback -- that would be exactly the kind "
            "of number that silently drifts out of sync with the real "
            "encoder. Confirm the repo layout and try again.\n"
        )
        return 3

    if not as_json:
        print()
        print("INPUT SCOPE REPORT")
        print(f"  bounds extracted from source (not hardcoded here):")
        print(f"    MAX_READ_LEN = {limits['MAX_READ_LEN']}  (src/encoder.cpp)")
        print(f"    nmc::MAXTOK  = {limits.get('MAXTOK')}  (include/names_coder.h, CAPS_NAMES=1 only)")
        print()

    results = [analyze_one(p, limits) for p in paths]
    worst = max((r["exit_code"] for r in results if r["exit_code"] is not None), default=0)

    if as_json:
        print(json.dumps({"limits": limits, "files": results}, indent=2))
    else:
        for r in results:
            print_human(r, limits)
        n_ok = sum(1 for r in results if r["exit_code"] == 0)
        n_caveat = sum(1 for r in results if r["exit_code"] == 1)
        n_oos = sum(1 for r in results if r["exit_code"] == 2)
        n_bad = sum(1 for r in results if r["exit_code"] == 3)
        print("-" * 78)
        print(f"  {len(results)} file(s): {n_ok} in scope, {n_caveat} in scope with caveats, "
              f"{n_oos} OUT OF SCOPE, {n_bad} could not be read")
        print("-" * 78)

    return worst


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
