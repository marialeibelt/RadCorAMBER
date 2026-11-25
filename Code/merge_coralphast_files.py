#!/usr/bin/env python3
"""Merge hand-picked PHAST ROOT files, even when they live in different folders.

Validate inputs and final tree entry counts; never delete original files.
"""
from __future__ import annotations
import argparse
import csv
import glob
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def read_path_list(file: Path) -> list[str]:
    if not file or not str(file):
        return []
    return [line.strip() for line in file.read_text(encoding='utf-8').splitlines()
            if line.strip() and not line.lstrip().startswith('#')]


def collect_files(args):
    paths = read_path_list(Path(args.input_list))
    if args.extra_list:
        other = Path(args.extra_list)
        if not other.is_file():
            raise RuntimeError(f'Extra list file not found: {other}')
        paths.extend(read_path_list(other))
    for pattern in args.glob_list.splitlines():
        pattern = pattern.strip()
        if not pattern:
            continue
        found = sorted(glob.glob(os.path.expanduser(pattern)))
        if not found:
            raise RuntimeError(f'File glob matched nothing: {pattern}')
        paths.extend(found)
    if not paths:
        raise RuntimeError('No ROOT inputs configured: add paths to INPUT_FILES.')
    normalized = [Path(os.path.expandvars(os.path.expanduser(p))).resolve() for p in paths]
    if len(set(normalized)) != len(normalized):
        raise RuntimeError('Repeated input paths detected; refusing double counting.')
    for p in normalized:
        if not p.is_file() or p.stat().st_size == 0:
            raise RuntimeError(f'Missing or empty ROOT input: {p}')
    return normalized


def read_root_tree(fpath: Path, tree_name: str, root):
    f = root.TFile.Open(str(fpath), 'READ')
    if not f or f.IsZombie():
        raise RuntimeError(f'ROOT file cannot be opened: {fpath}')
    try:
        t = f.Get(tree_name)
        if not t or not t.InheritsFrom('TTree'):
            raise RuntimeError(f'Tree {tree_name!r} missing from {fpath}')
        n = int(t.GetEntries())
        schema = tuple(sorted((str(b.GetName()), str(b.GetClassName()), str(b.GetTitle()))
                              for b in t.GetListOfBranches()))
        return n, schema
    finally:
        f.Close()


def hadd(target: Path, src: list[Path]):
    cmd = ['hadd', '-f', str(target)] + [str(p) for p in src]
    print(f'Running hadd: {len(src)} input files -> {target.name}', flush=True)
    subprocess.run(cmd, check=True)
    if not target.is_file() or target.stat().st_size == 0:
        raise RuntimeError(f'No merged ROOT data written: {target}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-list', required=True)
    p.add_argument('--glob-list', default='')
    p.add_argument('--extra-list', default='')
    p.add_argument('--out', required=True)
    p.add_argument('--tree', default='UserEvent300/PRM')
    p.add_argument('--chunk-size', type=int, default=30)
    p.add_argument('--overwrite', type=int, choices=(0,1), default=0)
    p.add_argument('--keep-chunks', type=int, choices=(0,1), default=0)
    p.add_argument('--expected-files', type=int, default=0)
    p.add_argument('--expected-entries', type=int, default=0)
    a = p.parse_args()
    if a.chunk_size < 2 or a.expected_files < 0 or a.expected_entries < 0:
        p.error('Chunk size must be >= 2, expected numbers must be >= 0')
    inputs = collect_files(a)
    target = Path(a.out).expanduser().resolve()
    if target in inputs:
        raise SystemExit('ERROR: Output file is also in the input list.')
    if target.exists() and not a.overwrite:
        raise SystemExit(f'ERROR: Output already exists: {target} (set OVERWRITE=1 to replace)')
    if a.expected_files and len(inputs) != a.expected_files:
        raise SystemExit(f'ERROR: Expected {a.expected_files} files, found {len(inputs)}')
    if not shutil.which('hadd'):
        raise SystemExit('ERROR: ROOT hadd not found; source ROOT_SETUP.')
    try:
        import ROOT
    except ImportError as e:
        raise SystemExit(f'ERROR: PyROOT unavailable: {e}')
    ROOT.gROOT.SetBatch(True)
    target.parent.mkdir(parents=True, exist_ok=True)
    counts, schema_ref = [], None
    for idx, f in enumerate(inputs, 1):
        n, schema = read_root_tree(f, a.tree, ROOT)
        if schema_ref is not None and schema != schema_ref:
            raise RuntimeError(f'Branch-schema mismatch at input {idx}: {f}')
        schema_ref = schema
        counts.append(n)
        if idx % 100 == 0 or idx == len(inputs):
            print(f'Validated {idx}/{len(inputs)} ROOT files', flush=True)
    expected_total = sum(counts)
    if a.expected_entries and expected_total != a.expected_entries:
        raise RuntimeError(f'Entries sum {expected_total:,} != expected {a.expected_entries:,}')
    print(f'Inputs: {len(inputs):,}; actual combined entries: {expected_total:,}', flush=True)
    tmp = Path(tempfile.mkdtemp(prefix='.custom_phast_merge_', dir=target.parent))
    result_tmp = tmp / 'merged_final.root'
    try:
        level = [Path(x) for x in inputs]
        step = 0
        # Recurse to <= chunk_size inputs; this avoids one hadd with thousands of arguments.
        while len(level) > a.chunk_size:
            out_level = []
            for start in range(0, len(level), a.chunk_size):
                chunk_path = tmp / f'level{step}_chunk{start//a.chunk_size:05d}.root'
                hadd(chunk_path, level[start:start+a.chunk_size])
                out_level.append(chunk_path)
            level = out_level
            step += 1
        hadd(result_tmp, level)
        actual, schema_out = read_root_tree(result_tmp, a.tree, ROOT)
        if actual != expected_total or schema_out != schema_ref:
            raise RuntimeError(f'Merged file validation failed: expected {expected_total}, got {actual}')
        manifest = target.with_suffix('.manifest.tsv')
        stats = target.with_suffix('.merge_summary.txt')
        # validated file is committed only once it is complete
        os.replace(result_tmp, target)
        with manifest.open('w', newline='', encoding='utf-8') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(['source_path','tree_entries'])
            for path, n in zip(inputs, counts):
                w.writerow([str(path),n])
        stats.write_text(
            f'Final ROOT: {target}\nTree: {a.tree}\nInput files: {len(inputs):,}\n'
            f'Summed input entries: {expected_total:,}\nMerged entries: {actual:,}\n'
            'No input files have been deleted.\n', encoding='utf-8')
        print(f'SUCCESS: {target}\nTree entries: {actual:,}\nManifest: {manifest}\nStats: {stats}', flush=True)
        if a.keep_chunks:
            print('Keeping intermediate ROOT chunks:', tmp)
            return
    finally:
        if not a.keep_chunks:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    try:
        main()
    except (OSError, RuntimeError, subprocess.CalledProcessError) as e:
        sys.exit(f'ERROR: {e}')
