# Architecture, layer by layer: G_CAPSUL vs DiscoSNP++

Third column answers one question only: **is this layer textbook (free to any
implementer), theirs (a DiscoSNP++/GATB design decision we adopted), or ours?**
Written from `~/DiscoSnp/tools/kissnp2/src/*.cpp` and our own source.

## The table

| # | layer | DiscoSNP++ (kissnp2 / GATB) | G_CAPSUL (Method B) | provenance |
|---|---|---|---|---|
| 1 | k-mer counting | GATB: minimizer-partitioned superkmers, disk-partitioned | same shape, but the table is **already built by the compressor** for its own coverage statistic | **textbook** (KMC2/GATB) — our reuse is ours |
| 2 | minimizer ordering | frequency-ranked (partition balance) | lexicographic; frequency ranking implemented, **measured +14% volume, rejected** | **theirs**, and measured not to pay here |
| 3 | memory budget | `-max-memory`, declared | ceiling = 60% of measured MemAvailable, spill self-enables | **textbook** systems practice |
| 4 | graph representation | cascading-Bloom dBG, GATB | flat hash table (`kc`) + `kc_find` as membership oracle | **textbook** |
| 5 | bubble start | branching node, >=2 successors, both strands | same, both orientations | **textbook** |
| 6 | SNV closure | lockstep `expand`, close on `nextNode1==nextNode2` | lockstep pairwise walk | **textbook** — same algorithm |
| 7 | multi-polymorphism | `-P 3` | `MAXPOLY=1` (structural: k=31 vs ~1/1000 het) | **theirs**, ours narrower and we still win |
| 8 | branching policy | `checkBranching`, `b=0` | `branch_walk` stops at first junction — equivalent | **textbook**, arrived at independently |
| 9 | indel proposal | `start_indel_prediction`: BFS per successor pair, **unconditional** | superbubble (Onodera) only — fires at **15.9%** of nodes | **theirs**; our gap, ported opt-in |
| 10 | indel closure shape | guaranteed by lockstep traversal | `LCP+LCS >= |s_short|` test, added after tracing theirs | **theirs** (geometry), **ours** (the explicit test) |
| 11 | indel ambiguity | `checkRepeatSize`, reject `k-2-min(ext) > 20` | ported; **monotone loss on our data**, kept at their default | **theirs**, measured |
| 12 | low complexity | DUST; **off by default** (`l="-l"`) | absent | **theirs**, not a gap |
| 13 | read coherence | **kissreads2 — a separate tool, second pass over all reads** | reads are **already resident**; one indexed sweep, plus a 1-bit quality bitmap | **ours** |
| 14 | quality use | mean phred per path | per-base bitmap (`MINQ=20`), 233 MB vs 2.27 GB of phred strings | **ours** |
| 15 | parallelism | GATB thread pool throughout | traversal was serial under superbubble mode; now parallel, **2.9x**, counts identical | **textbook**, was our deficit |
| 16 | coverage ceiling | none equivalent | `COVCAP = 2*PLOIDY*H`, derived from measured depth | **ours** |
| 17 | ploidy handling | none | HETSCAN pair-fraction gate, declines on haploids (E. coli 0.016) | **ours** |
| 18 | output | FASTA bubbles -> VCF_creator | VCF direct, plus a lossless archive from the same pass | **ours** |

## Why we are different from them, not just similar

**1. The graph is a byproduct, not a build step.** DiscoSNP++ constructs a dBG
in order to call variants. We construct `kc` to compute one scalar the
*compressor* needs, and the graph is what that table already is. The
comparison is not "two callers"; it is "a caller" versus "a compressor that
answers the same question for free". Layers 1, 4 and 13 are all cheaper for
this reason, and it is why the archive and the VCF come out of one pass.

**2. The reads are still in memory, so read coherence is not a second tool.**
kissreads2 exists because a standalone caller has thrown the reads away by the
time it has bubbles. We have not. That converts their separate pass into one
indexed sweep and makes a *per-base* quality test affordable where they use a
per-path mean — a strictly stronger filter (layers 13, 14). This is the layer
where being a compressor is a genuine algorithmic advantage rather than a
packaging difference.

**3. Our parameters are derived from measured depth; theirs are constants.**
`COHC = max(2, H/10)` and `COVCAP = 2*PLOIDY*H` are functions of the sample's
own haploid depth and ploidy. That is why precision holds 0.911-0.950 across a
10-30x sweep and why the tetraploid arm works without retuning. DiscoSNP++
exposes `-b`, `-P`, `-D`, `-max_ambigous_indel` as fixed numbers.

**4. We win where the shared algorithm is shared, and lose where theirs is
richer.** Layers 5-8 are the same textbook bubble calling, and on SNVs we beat
them (0.879 vs 0.847) — that is the fair comparison, and it isolates our
contribution to layers 13-17. On indels their layer 9 attempts a call at every
branching node while ours attempts one at 15.9% of them; that is a real
architectural deficit, honestly a place where their design is better, and it is
why the indel claim is withdrawn rather than argued.

## Genuine architectural headroom

**The one substantial, unexploited thing we own: free read threading.**

McCortex's Linked de Bruijn Graph augments a dBG with long-range connectivity by
**threading reads back through the graph in a separate pass**, and the cost is
real — the published figure is 20 GiB for links on top of 50 GiB for the graph
([Turner et al., *Bioinformatics* 2018](https://academic.oup.com/bioinformatics/article/34/15/2556/4938484)).
Link construction is a whole pipeline stage.

Our pseudogenome construction **already computes, for every read, which contig
it lies on and at which offset** (`ppos`, `read_cid`) — because that is what
mapping reads onto the pseudogenome means. That is precisely the read-to-path
placement McCortex pays a separate pass for, and we currently throw it away
after compression.

Wiring it into the caller would give long-range connectivity at essentially zero
marginal cost, and it attacks the exact failure mode we measured: `find_sb`
exhausts on 960,541 of 1,225,194 branching nodes because it cannot decide which
way to go through a repeat. A read placement says which way an actual molecule
went. This is the strongest remaining lever and it is ours by construction, not
borrowed.

**Prior art to be careful about.** "Searchable compressed archive of reads" is
not new — BEETL-fastq
([Cox et al., *Bioinformatics* 2014](https://academic.oup.com/bioinformatics/article/30/19/2796/2422232))
already does selective read extraction and variant-related queries from a
compressed archive. Claim 3 should be positioned against it explicitly rather
than claiming addressability as novel. The defensible framing is the
*combination*: lossless compression, a retained assembly, and reference-free
calling from one pass — not any one of those alone.

**Two smaller levers, both measured as open rather than assumed:**
- Calling from a *stored* archive (today it is compress-time only) — this is
  what would make "addressable" a demonstrated property rather than a design
  intent.
- Partition balance from layer 2's frequency ranking is still unmeasured; it was
  rejected on volume, which is not what it is for. If the merge is ever made
  memory-bound per partition, re-measure it there.

Sources:
- [Integrating long-range connectivity information into de Bruijn graphs (McCortex)](https://academic.oup.com/bioinformatics/article/34/15/2556/4938484)
- [BEETL-fastq: a searchable compressed archive for DNA reads](https://academic.oup.com/bioinformatics/article/30/19/2796/2422232)
