#!/usr/bin/env python3
"""Run one experiment, recording GNU time plus sampled process-tree PSS.

Usage: python3 scripts/gpt2026_measure.py NEW_OUTPUT_DIR -- COMMAND ARG...
PSS counts shared pages proportionally; summing RSS double-counts forked pages.
The half-second sample is a lower bound on peak PSS, not an exact high-water mark.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def gpt2026_cpu():
    values = list(map(int, Path('/proc/stat').read_text().splitlines()[0].split()[1:]))
    return sum(values[:8]), values[3] + values[4]


def gpt2026_tree(pid):
    todo, seen = [pid], set()
    while todo:
        current = todo.pop()
        if current in seen:
            continue
        seen.add(current)
        try:
            todo.extend(map(int, Path(f'/proc/{current}/task/{current}/children').read_text().split()))
        except (FileNotFoundError, ProcessLookupError):
            pass
    return seen


def gpt2026_memory(pid):
    rss, pss, count = 0, 0, 0
    for current in gpt2026_tree(pid):
        try:
            entries = Path(f'/proc/{current}/smaps_rollup').read_text().splitlines()
        except (FileNotFoundError, ProcessLookupError):
            continue
        except PermissionError:
            return None
        count += 1
        for line in entries:
            if line.startswith('Rss:'):
                rss += int(line.split()[1])
            elif line.startswith('Pss:'):
                pss += int(line.split()[1])
    return {'rss_kib': rss, 'pss_kib': pss, 'processes': count}


def gpt2026_main():
    if len(sys.argv) < 4 or sys.argv[2] != '--':
        raise SystemExit(__doc__)
    output = Path(sys.argv[1]).resolve()
    if output.exists():
        raise SystemExit(f'Refusing to overwrite existing experiment: {output}')
    if shutil.disk_usage(output.parent).free < 40 * 1024**3:
        raise SystemExit('Need at least 40 GiB free before an experiment')
    before = gpt2026_cpu()
    time.sleep(2)
    after = gpt2026_cpu()
    busy = 1 - (after[1]-before[1]) / max(1, after[0]-before[0])
    # This is an operational idle check, not an encoder parameter.
    busy_cores = busy * (os.cpu_count() or 1)
    if busy_cores > 0.6:
        raise SystemExit(f'Machine not idle: {busy_cores:.2f} busy CPU equivalents; no timed job started')
    output.mkdir()
    metadata = {'command': sys.argv[3:], 'cwd': os.getcwd(),
                'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                'cpu_count': os.cpu_count(), 'affinity': sorted(os.sched_getaffinity(0)),
                'pre_run_busy_cores': busy_cores,
                'free_bytes': shutil.disk_usage(output).free,
                'pss_sample_interval_s': 0.5,
                'env': {k: v for k, v in os.environ.items()
                        if k.startswith(('GPT2026_', 'CAPS_', 'OMP_')) or k in ('INPUT','ARCHIVE','BEST')}}
    (output/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    samples = []
    with (output/'stdout.txt').open('wb') as stdout, (output/'stderr.txt').open('wb') as stderr:
        proc = subprocess.Popen(['/usr/bin/time', '-v', '-o', str(output/'time.txt'), *sys.argv[3:]],
                                stdout=stdout, stderr=stderr)
        start = time.monotonic()
        while proc.poll() is None:
            memory = gpt2026_memory(proc.pid)
            if memory is not None:
                samples.append({'elapsed_s': time.monotonic()-start, **memory})
            time.sleep(0.5)
        rc = proc.wait()
    (output/'memory_samples.json').write_text(json.dumps(samples)+'\n')
    summary = {'returncode': rc, 'sampled_tree_peak_pss_kib': max((x['pss_kib'] for x in samples), default=None),
               'sampled_tree_peak_rss_kib': max((x['rss_kib'] for x in samples), default=None),
               'pss_is_sampled_lower_bound': True}
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary))
    return rc


if __name__ == '__main__':
    sys.exit(gpt2026_main())
