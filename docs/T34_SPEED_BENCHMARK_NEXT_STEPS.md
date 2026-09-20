# T3.4 speed vs BEETL — investigated, not run. Future work.

**Status: NOT started. This is a plan and a disk-space finding, nothing more.
T3.4's current claim (exact-match recall = 1.0000) is complete and honest as
stated and does NOT depend on this. Do not treat this file as a pending gap in
T3.4 — it is an optional future extension, only needed if a speed comparison
is added on top of the existing recall claim.**

---

## Why this exists

T3.4 currently claims recall only. CLAUDE.md already flags the missing piece
honestly: *"We have not benchmarked against BEETL/CIndex/sFASTQ — that is the
most valuable missing experiment."* That sentence is about **adding a speed
claim**, not about the recall claim being incomplete. The two are separable:

- Recall claim (done, locked): exact-match recall 1.0000 at k=3, own
  measurement, no competitor needed to state it.
- Speed claim (not done, optional): "how fast are we vs BEETL at the same
  exact-match task" — only meaningful if BEETL can do the same job, which it
  can (unlike T3.5/position retrieval, where no competitor exists at all).

## What was checked before deciding to defer

**Disk space, verified, not assumed:**

- Real Windows C: drive free space: measured directly via PowerShell
  (`Get-PSDrive`), **not** via WSL's own `df -h` (which reported a misleading
  "910GB free" — that number is WSL's ext4 filesystem's configured max size
  inside its sparse vhdx container, not real host disk space).
- The actual constraint is the vhdx **container file** on the real C: drive.
  Found at `C:\Users\Lenovo\AppData\Local\wsl\{a6bfff6e-7e13-4353-9bba-b24506722248}\ext4.vhdx`.
- Freed ~14GB of confirmed-dead WSL directories (`bench_5x5_300x`,
  `bench_5x5`, `arcs_speed`, `ggcat_src`, duplicate `miniconda` — verified via
  grep against the whole repo that nothing live references them; `miniconda3`,
  the one WITH the "3", is live — it's on the PATH in every Claim 2 benchmark
  script and must never be deleted).
- That 14GB is freed *inside* WSL's filesystem but is still locked inside the
  vhdx container on the real C: drive — reclaiming it on the Windows side
  needs `diskpart compact vdisk`, which requires an elevated (Administrator)
  PowerShell session this session does not have.
- **Real usable budget right now: 13.9 GB free on C:.**
- Elevation command, to run manually when convenient (safe, no risk):
  ```powershell
  wsl --shutdown
  diskpart
  ```
  then inside diskpart:
  ```
  select vdisk file="C:\Users\Lenovo\AppData\Local\wsl\{a6bfff6e-7e13-4353-9bba-b24506722248}\ext4.vhdx"
  attach vdisk readonly
  compact vdisk
  detach vdisk
  exit
  ```
  This reclaims the ~14GB safely. A non-admin alternative (`wsl --manage
  Ubuntu --set-sparse --allow-unsafe`) exists but was explicitly rejected —
  Microsoft's own docs flag it as carrying a data-corruption risk, and this
  project does not take that kind of shortcut.

**Dataset footprint, verified, not assumed:**

- T3.1/T3.2's 6 locked datasets: `ERR5181310`, `SRR2584863`, `SRR29296997`,
  `SRR37283774`, `DRR976266`, `HG002`.
- Only `SRR2584863` is currently on disk (~358MB raw, both mates). The other
  5 were already cleaned up after their original T3.1/T3.2 run — nothing to
  reclaim, they would need re-downloading regardless of whether this
  extension happens.
- 4 of the 5 missing ones (`ERR5181310`, `SRR29296997`, `SRR37283774`,
  `DRR976266`) are almost certainly small (same order of magnitude as
  `SRR2584863`, going by their SPAdes/coverage runtimes already measured in
  `benchmark/results/claim3/claim3_T3.1_T3.2_T3.3.csv`).
- `HG002` is the outlier — 300x coverage over chr21+22, its own SPAdes run
  took 2,679s versus 50-550s for the others. Its real size must be re-checked
  immediately before running it, not assumed from the others.

## The plan, if this is picked up later

1. Build BEETL once (small footprint, a few hundred MB).
2. Process datasets **one at a time**, smallest/safest first: `SRR2584863` →
   `ERR5181310` → `SRR29296997` → `SRR37283774` → `DRR976266` → `HG002` last.
3. Per dataset: download raw reads if needed -> build BEETL's index -> run the
   same exact-match query workload already used for our own recall number ->
   record wall time for both -> delete that dataset's files -> move to next.
4. Never hold more than one dataset's raw data + index on disk at once. Keeps
   peak usage far under the 13.9GB (or ~28GB post-compaction) ceiling
   regardless of dataset count.
5. Re-check `HG002`'s real size right before its turn, since it is the one
   dataset that could actually threaten the budget.
6. Result feeds **T3.4 only** (both tools can do exact-match retrieval, so a
   speed comparison is fair). It does **not** feed T3.5 — no competitor exists
   for position retrieval, by design, so there is nothing to race there.

## What NOT to do

- Do not present the absence of this benchmark as a gap in T3.4's existing
  recall claim. The recall claim stands on its own.
- Do not run this without re-verifying real disk headroom at the time it is
  picked up — this file's numbers (13.9GB free, 51.52GB vhdx) are a snapshot
  from 2026-09-19 and will have drifted.
- Do not force the unsafe sparse-vhd flag to save the admin-elevation step.
