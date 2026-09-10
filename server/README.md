# server/ — rebuild this project on a fresh machine

Two documents, deliberately:

| file | when to read it |
|---|---|
| **`SETUP_SHORT.md`** | you want it running. Copy-paste, ~2 h wall, mostly downloads. |
| **`SETUP_FULL.md`** | something failed, or you need to know *why* a step is the way it is. Every gotcha this project actually hit, with the exact error text so you can search for it. |
| **`DATASETS.md`** | what data to fetch, where it goes, how big, and how to verify it. |
| **`TOOLS.md`** | every external tool, its exact version here, and how it was installed. |

**Read `SETUP_FULL.md` before improvising.** Several steps look wrong and are
not: the SRA toolkit must be symlinked rather than copied, SPRING needs a
source patch under GCC 13, Genozip needs a network-reachable `/dev/stdout`, and
Kmer2SNP must bypass its own dependency wrappers. Each of those cost a debugging
session; they are written down so they cost you nothing.

## What this project is

Two repositories, one paper:

- **ARCS** (`/root/arcs-clean`, `github.com/thackshanaramana0-spec/ARCS`) — the
  paper project. Contains the locked dataset definitions, the earlier benchmark
  harness, and `method_c/` (a clone of PgRC2, used as a comparison baseline,
  GPL-3, never vendored into our source).
- **G_CAPSUL / Capsule** (`/root/arcs-clean/c_star_pg_advance`,
  `github.com/thackshanaramana0-spec/g_capsul`) — the compressor, the caller,
  the final benchmark harness, and every result quoted in the paper. This is
  where the work is.

ARCS's `.gitignore` excludes `c_star_pg_advance/`, so they are independent
checkouts that happen to nest. Clone both.

## The five-minute version

    tmux new -s work                      # do this FIRST, see SETUP_FULL 1
    git clone https://github.com/thackshanaramana0-spec/ARCS.git /root/arcs-clean
    git clone https://github.com/thackshanaramana0-spec/g_capsul.git \
              /root/arcs-clean/c_star_pg_advance
    cd /root/arcs-clean/c_star_pg_advance
    bash server/install_tools.sh          # everything apt-installable + built
    # then: Genozip activation is INTERACTIVE, see SETUP_FULL 3.7
    bash server/fetch_data.sh             # ~55 GB, hours
    bash scripts/benchmark_0_preflight.sh # must print VERDICT: GO
    bash scripts/benchmark_1_run.sh       # the full sweep, ~5 h

## Machine requirements, measured not guessed

| | required | why |
|---|---|---|
| cores | 8+ | 12 used here; fewer works and is slower, output is identical |
| RAM | **32 GB min, 64 GB comfortable** | peak 10.9 GB (compress SRR10676752), 19.9 GB (caller on HG005) |
| disk | **250 GB** | 55 GB FASTQ + ~7 GB archives + ~36 GB peak transient |
| OS | Ubuntu 24.04 | GCC 13; the SPRING patch in SETUP_FULL is specific to it |

`benchmark_0_preflight.sh` fails the run if RAM < 32 GB or free disk < 60 GB.
That threshold is not arbitrary: a full disk once silently truncated a k-mer
spill and produced a *fake* F1 regression that cost a day. See
`SETUP_FULL.md` §6.

## Where the results live

Nothing in `server/` produces a number. Results, their provenance, and the
scripts that generated them are in `../benchmark/`, and every reported value is
traced end to end in `../benchmark/documentation/RESULT_CODE.md`.
