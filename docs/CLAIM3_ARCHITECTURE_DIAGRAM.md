# Claim 3 (ADDRESSABLE) — architecture diagram, source

**Status: DRAFT CONTENT ONLY. Not wired into the manuscript, not finalized.**
Written so the diagram can be rendered/refined without re-deriving the code
trace. Every box below is anchored to an exact function or line range in
`stages/capsule_decode.cpp`, verified this session, not recalled from memory.

---

## The real shape (not five parallel boxes)

Five published tables (T3.1-T3.5) come out of **four actual code paths**,
because T3.3 and T3.4 are the same function with different inputs, and T3.5
is not in the C++ at all — it is an orchestration pattern built by calling
T3.4's primitive twice.

```mermaid
flowchart TD
    ARC[("*.capsule archive")]

    ARC -->|"mode==export"| EXP["EXPORT<br/>contig_spans stream<br/>(LEB128 gap,len pairs)"]
    EXP --> EXPOUT["Per-contig FASTA<br/>slice pg.data()+offset directly<br/>NO re-decode, NO FASTA round-trip<br/>[capsule_decode.cpp:903-957]"]

    ARC -->|"mode==coverage"| COV["COVERAGE<br/>EXITS BEFORE pg REBUILD<br/>[capsule_decode.cpp:417-490]"]
    COV --> COVIN["reads only:<br/>pos_abs + read_lengths + orig2uid"]
    COVIN --> COVALGO["sweep-line diff array<br/>diff[a]++, diff[b]--<br/>one pass over PGLEN"]
    COVALGO --> COVOUT["per-base depth,<br/>pg_main / pg_second rows"]
    COVNOTE["real fixed bug:<br/>orig2uid expansion required,<br/>else 20% undercount on duplicates<br/>(measured on E. coli)"] -.-> COVIN

    ARC -->|"mode==query, needs .qidx sidecar"| Q["build sidecar once:<br/>capsule_decode index<br/>2-bit packed pg + placements P,L<br/>+ optional deviations, strand, N-pos"]

    Q --> QDISPATCH{"how is rr[] (start,end)<br/>populated?"}

    QDISPATCH -->|"literal START-END string"| T33["T3.3 — RANGE QUERY<br/>strtoull(modearg,'-')<br/>[line 640-642]<br/>no search, coordinates ARE the input"]

    QDISPATCH -->|"ACGTN probe string"| T34CORE["find_occurrences()<br/>mismatch-tolerant search<br/>over 2-bit packed pg<br/>k_max = floor(P/MINSEED)-1"]

    T34CORE -->|"CAPS_QUERY_CONTAIN=1"| XMI[".xmi completion<br/>k-mer index over<br/>DEVIATION-CARRYING reads only<br/>(~14% of reads)<br/>[line 644-720]"]
    XMI --> T34["T3.4 — EXACT MATCH<br/>UNION(pg-search, xmi-hits)<br/>every candidate verified against<br/>fully-reconstructed read<br/>(strand + N applied)<br/>recall 1.0000, proof not luck"]
    T34CORE -->|"CAPS_QUERY_CONTAIN unset"| T34PLAIN["plain pg-search only<br/>recall 0.86-0.96<br/>(varies by individual)"]

    T33 --> EMIT["shared emit loop<br/>reads placed within rr[] ranges,<br/>deviations + strand + N applied"]
    T34 --> EMIT
    T34PLAIN --> EMIT

    EMIT --> T35ORCH["ORCHESTRATION LAYER<br/>(scripts/t34_realdata/run_window.sh,<br/>score_bilateral.py)<br/>NOT inside capsule_decode.cpp"]
    T35ORCH --> T35A["query, probe UPSTREAM of variant"]
    T35ORCH --> T35B["query, probe DOWNSTREAM of variant"]
    T35A --> T35UNION["UNION in Python<br/>best-overlap read assignment"]
    T35B --> T35UNION
    T35UNION --> T35["T3.5 — LOCUS RETRIEVAL<br/>1,451/1,452 native (99.93%)<br/>1,452/1,452 with .xmi (100%)<br/>fixes the case a single coordinate<br/>structurally cannot: het allele split"]

    style ARC fill:#2b2b2b,color:#fff
    style EXPOUT fill:#1a4d2e,color:#fff
    style COVOUT fill:#1a4d2e,color:#fff
    style T33 fill:#1a3d5c,color:#fff
    style T34 fill:#1a3d5c,color:#fff
    style T35 fill:#5c1a3d,color:#fff
```

---

## Why the shape matters (the point of drawing it this way)

A naive diagram — five boxes side by side, one per table — would imply five
bespoke mechanisms, which is not what the code does and undersells the
actual claim. The real structure is:

1. **Two fast-exit paths** (export, coverage) that skip the assembly layer
   entirely, because they don't need pg *content* — only its length,
   placements, or precomputed spans.
2. **One shared query primitive**, used three different ways depending only
   on what fills `rr[]` — literal coordinates, a pg-search, or a pg-search
   unioned with an auxiliary index.
3. **One orchestration layer** (T3.5) that calls the primitive twice and
   reconciles the results — this is where the paper's actual insight
   (allele-splitting) gets fixed, and it is built entirely on top of
   mechanism #2, not a new one.

That's a stronger, truer picture than "five features": one retained
structure, reused through one primitive, in increasingly demanding ways —
which is the sentence the abstract already makes. The diagram should prove
it, not just illustrate five tables.

## What still needs deciding before this is final

- Visual style (this is Mermaid syntax — renders as a flowchart; may need
  redrawing in whatever tool matches Claim 1/Claim 2's diagrams for visual
  consistency).
- Whether to show the `.qidx` sidecar build step explicitly (currently
  folded into the "build sidecar once" box) or expand it, given `CAPS_XMI`
  and `CAPS_PILEUP` are separate opt-in flags at index-build time, not part
  of every sidecar.
- Level of code-citation detail to keep visible in the final image vs. move
  to a caption — current draft keeps line numbers inline for traceability
  during review; likely too dense for the camera-ready version.

Not wired into `capsul_manuscript.tex`. Not rendered as an image file yet.
