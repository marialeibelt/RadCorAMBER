#!/usr/bin/env python3
"""Count expected TGEANT events and safely merge matching PHAST ROOT outputs.

Only Python standard library required. ROOT/hadd are used by the merge action.
Nominal event counts are inferred from '_batchN_COUNT_runM.tgeant.gz' names,
not by decompressing large TGEANT files. Actual PHAST tree counts are accepted
when smaller than nominal and are verified exactly through the final merge.
"""
import datetime
import fnmatch
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

BATCH_RE = re.compile(r'_batch\d+_(\d+)_run\d+\.tgeant\.gz$')


def die(message):
    raise RuntimeError(message)


def env(name):
    value = os.environ.get(name, '').strip()
    if not value:
        die(f'Missing configuration variable: {name}')
    return value


def integer(name, minimum=0):
    value = env(name)
    try:
        n = int(value)
    except ValueError:
        die(f'{name} must be an integer, got {value!r}')
    if n < minimum:
        die(f'{name} must be >= {minimum}')
    return n


def natural_key(filename):
    return [int(x) if x.isdigit() else x.lower()
            for x in re.split(r'(\d+)', str(filename))]


def listed_files(directory, pattern):
    p = Path(directory)
    if not p.is_dir():
        die(f'Directory does not exist: {p}')
    return sorted(
        (f for f in p.iterdir() if f.is_file() and fnmatch.fnmatchcase(f.name, pattern)),
        key=lambda f: natural_key(f.name)
    )


def expected_from_name(name):
    match = BATCH_RE.search(name)
    if not match:
        die(f'Cannot infer event count from TGEANT filename: {name}\n'
            'Expected suffix: _batchNN_NNNN_run000.tgeant.gz')
    count = int(match.group(1))
    if count < 1:
        die(f'Nonpositive event count in filename: {name}')
    return count


def tgeant_index(directory, pattern):
    paths = listed_files(directory, pattern)
    if not paths:
        die(f'No matching TGEANT files in {directory} ({pattern})')
    index = {}
    digest = hashlib.sha256()
    total = 0
    for p in paths:
        if p.name in index:
            die(f'Duplicate TGEANT filename: {p.name}')
        stat = p.stat()
        if stat.st_size <= 0:
            die(f'Empty TGEANT input: {p}')
        nevt = expected_from_name(p.name)
        index[p.name] = nevt
        total += nevt
        # Detect modified/replaced source files between the count and merge steps.
        digest.update(f'{p.name}\t{stat.st_size}\t{stat.st_mtime_ns}\t{nevt}\n'.encode())
    return paths, index, total, digest.hexdigest()


def atomic_write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.tmp.{os.getpid()}')
    try:
        tmp.write_text(text, encoding='utf-8')
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def count():
    directory = env('TGEANT_DIR')
    pattern = env('TGEANT_PATTERN')
    expected_files = integer('EXPECTED_FILES', 1)
    report_path = Path(env('COUNT_REPORT'))
    files, index, total, digest = tgeant_index(directory, pattern)
    if len(files) != expected_files:
        die(f'Expected {expected_files} TGEANT files, found {len(files)}. '
            'Count report NOT written.')
    report = {
        'version': 1,
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'tgeant_dir': str(Path(directory).resolve()),
        'tgeant_pattern': pattern,
        'expected_files': len(files),
        'expected_events': total,
        'source_manifest_sha256': digest,
        'event_count_method': 'batch filename (_batchN_COUNT_runM.tgeant.gz); not file content',
    }
    atomic_write(report_path, json.dumps(report, indent=2) + '\n')
    atomic_write(report_path.with_name('expected_tgeant_filenames.txt'),
                 ''.join(f'{p.name}\t{index[p.name]}\n' for p in files))
    print('==============================================================')
    print('Expected PHAST events from TGEANT batch filenames')
    print('==============================================================')
    print('TGEANT directory:  ', directory)
    print('TGEANT files:      ', len(files))
    print('Expected events:   ', f'{total:,}')
    print('Count report:      ', report_path)
    print('Filename manifest: ', report_path.with_name('expected_tgeant_filenames.txt'))
    print('NOTE: This is the intended generated-event count encoded in filenames.')
    print('      Actual PHAST TTree entries are verified during merging.')


SCAN_MACRO = r'''
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

// Read all PHAST trees, record the ACTUAL entries in each file.
// A smaller entry count than the nominal batch size is allowed.
void scan_phast_roots(const char *inputManifest, const char *treeName,
                      const char *outputManifest) {
    std::ifstream in(inputManifest);
    std::ofstream out(outputManifest, std::ios::trunc);
    if (!in || !out) {
        std::cerr << "ERROR: Cannot open scan manifest(s)" << std::endl;
        gSystem->Exit(30); return;
    }
    std::string line;
    Long64_t totalNominal = 0, totalActual = 0, nFiles = 0, nReduced = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto split = line.rfind('\t');
        if (split == std::string::npos) {
            std::cerr << "ERROR: Invalid input manifest line" << std::endl;
            gSystem->Exit(31); return;
        }
        const std::string filename = line.substr(0, split);
        Long64_t nominal = -1;
        try { nominal = std::stoll(line.substr(split + 1)); }
        catch (...) {
            std::cerr << "ERROR: Invalid nominal count: " << filename << std::endl;
            gSystem->Exit(32); return;
        }
        TFile *file = TFile::Open(filename.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "ERROR: Invalid ROOT file: " << filename << std::endl;
            gSystem->Exit(33); return;
        }
        TTree *tree = dynamic_cast<TTree*>(file->Get(treeName));
        if (!tree) {
            std::cerr << "ERROR: Missing TTree " << treeName << " in " << filename << std::endl;
            gSystem->Exit(34); return;
        }
        const Long64_t actual = tree->GetEntries();
        if (actual <= 0 || actual > nominal) {
            std::cerr << "ERROR: Unexpected entries in " << filename
                      << " nominal=" << nominal << " actual=" << actual << std::endl;
            gSystem->Exit(35); return;
        }
        if (actual < nominal) ++nReduced;
        out << filename << '\t' << actual << '\n';
        if (!out) {
            std::cerr << "ERROR: Cannot write scan results" << std::endl;
            gSystem->Exit(36); return;
        }
        totalNominal += nominal;
        totalActual += actual;
        ++nFiles;
        file->Close();
        delete file;
        if (nFiles % 500 == 0)
            std::cout << "Scanned " << nFiles << " PHAST ROOT files" << std::endl;
    }
    out.close();
    if (!out || nFiles == 0) {
        std::cerr << "ERROR: Incomplete or empty scan output" << std::endl;
        gSystem->Exit(37); return;
    }
    std::cout << "PHAST ROOT files scanned: " << nFiles << std::endl;
    std::cout << "Nominal generated events: " << totalNominal << std::endl;
    std::cout << "Actual PHAST tree entries: " << totalActual << std::endl;
    std::cout << "Events absent before PHAST: " << totalNominal - totalActual << std::endl;
    std::cout << "Files with fewer entries: " << nReduced << std::endl;
    gSystem->Exit(0);
}
'''

VALIDATE_MACRO = r'''
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

void validate_phast_roots(const char *manifest, const char *treeName, Long64_t totalExpected) {
    std::ifstream in(manifest);
    if (!in) {
        std::cerr << "ERROR: Cannot open validation manifest " << manifest << std::endl;
        gSystem->Exit(20);
        return;
    }
    std::string line;
    Long64_t nTotal = 0;
    Long64_t nFiles = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const auto pos = line.rfind('\t');
        if (pos == std::string::npos) {
            std::cerr << "ERROR: Invalid validation manifest line: " << line << std::endl;
            gSystem->Exit(21);
            return;
        }
        std::string filename = line.substr(0, pos);
        Long64_t expected = -1;
        try { expected = std::stoll(line.substr(pos+1)); }
        catch (...) {
            std::cerr << "ERROR: Invalid count: " << line << std::endl;
            gSystem->Exit(22);
            return;
        }
        TFile *f = TFile::Open(filename.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "ERROR: Cannot open ROOT: " << filename << std::endl;
            gSystem->Exit(23);
            return;
        }
        TTree *t = dynamic_cast<TTree*>(f->Get(treeName));
        if (!t) {
            std::cerr << "ERROR: Missing TTree " << treeName << " in " << filename << std::endl;
            gSystem->Exit(24);
            return;
        }
        const Long64_t actual = t->GetEntries();
        if (actual != expected) {
            std::cerr << "ERROR: Entry mismatch in " << filename
                      << " expected=" << expected << " actual=" << actual << std::endl;
            gSystem->Exit(25);
            return;
        }
        nTotal += actual;
        ++nFiles;
        f->Close();
        delete f;
        if (nFiles % 500 == 0) {
            std::cout << "Validated " << nFiles << " input ROOT files" << std::endl;
        }
    }
    std::cout << "Validated ROOT files: " << nFiles << std::endl;
    std::cout << "Total tree entries:   " << nTotal << std::endl;
    std::cout << "Expected tree entries:" << totalExpected << std::endl;
    if (nTotal != totalExpected) {
        std::cerr << "ERROR: Total input entry mismatch" << std::endl;
        gSystem->Exit(26);
        return;
    }
    gSystem->Exit(0);
}
'''



def root_scan(manifest, tree_name, output_manifest, macro_path, expected_files):
    """Obtain observed entry counts for all PHAST inputs (losses are allowed)."""
    temporary = output_manifest.with_name(f'.{output_manifest.name}.scan.{os.getpid()}')
    try:
        temporary.unlink(missing_ok=True)
        args = (f'{macro_path}({json.dumps(str(manifest))},'
                f'{json.dumps(tree_name)},{json.dumps(str(temporary))})')
        subprocess.run(['root', '-l', '-b', '-q', args], check=True)
        if not temporary.is_file():
            die('ROOT scan did not produce an actual-count manifest')
        counts = {}
        with temporary.open(encoding='utf-8') as stream:
            for line in stream:
                path, sep, n = line.rstrip('\n').rpartition('\t')
                if not sep or not path or path in counts:
                    die(f'Invalid or duplicate PHAST scan row: {line!r}')
                val = int(n)
                if val <= 0:
                    die(f'Empty PHAST tree: {path}')
                counts[path] = val
        if len(counts) != expected_files:
            die(f'Incomplete scan: {len(counts)} PHAST trees, expected {expected_files}')
        os.replace(temporary, output_manifest)
        return counts
    finally:
        temporary.unlink(missing_ok=True)

def root_validate(manifest, tree_name, expected_total, macro_path):
    # ROOT calls a function matching the macro's basename.
    args = f'{macro_path}({json.dumps(str(manifest))},{json.dumps(tree_name)},{expected_total}LL)'
    subprocess.run(['root', '-l', '-b', '-q', args], check=True)


def run_hadd(output, inputs):
    print(f'Running hadd: {len(inputs)} inputs -> {output}', flush=True)
    subprocess.run(['hadd', '-f', str(output), *map(str, inputs)], check=True)
    if not output.is_file() or output.stat().st_size <= 0:
        die(f'hadd did not produce a nonempty ROOT file: {output}')


def root_validate_one(path, tree_name, expected, manifest_path, macro_path):
    atomic_write(manifest_path, f'{path}\t{expected}\n')
    root_validate(manifest_path, tree_name, expected, macro_path)


def merge():
    tgeant_dir = env('TGEANT_DIR')
    tgeant_pattern = env('TGEANT_PATTERN')
    input_dir = env('INPUT_DIR')
    input_pattern = env('INPUT_PATTERN')
    report_path = Path(env('COUNT_REPORT'))
    expected_files = integer('EXPECTED_FILES', 1)
    chunksize = integer('MERGE_CHUNK_SIZE', 1)
    use_existing = integer('USE_EXISTING_MERGE')
    overwrite = integer('OVERWRITE')
    delete_tgeant = integer('DELETE_TGEANT_AFTER_SUCCESS')
    if use_existing not in (0, 1) or overwrite not in (0, 1):
        die('USE_EXISTING_MERGE and OVERWRITE must be 0 or 1')
    if delete_tgeant not in (0, 1):
        die('DELETE_TGEANT_AFTER_SUCCESS must be 0 or 1')
    tree_name = env('TREE_NAME')
    outdir = Path(env('MERGED_DIR'))
    target = outdir / env('MERGED_FILE')
    if target.name != env('MERGED_FILE'):
        die('MERGED_FILE must be a filename, not an absolute/relative directory path.')
    if not report_path.is_file():
        die(f'Missing count report: {report_path}\nRun count_phast_expected_events.sh first.')
    report = json.loads(report_path.read_text(encoding='utf-8'))
    source_files, source_index, source_total, digest = tgeant_index(tgeant_dir, tgeant_pattern)
    if report.get('version') != 1:
        die('Unsupported expected-events report version')
    if report.get('tgeant_dir') != str(Path(tgeant_dir).resolve()):
        die('TGEANT_DIR differs between count and merge steps')
    if report.get('tgeant_pattern') != tgeant_pattern:
        die('TGEANT_PATTERN differs between count and merge steps')
    if report.get('source_manifest_sha256') != digest:
        die('TGEANT file list / sizes / timestamps changed since count step. Recount.')
    if (len(source_files) != expected_files or
            report.get('expected_files') != expected_files):
        die(f'Wrong expected TGEANT count: got {len(source_files)}, expected {expected_files}')
    if report.get('expected_events') != source_total:
        die('Expected event total changed since count step')
    phast_files = listed_files(input_dir, input_pattern)
    if len(phast_files) != expected_files:
        die(f'Wrong PHAST file count: found {len(phast_files)}, expected {expected_files}. '
            'Finish or repair the reconstruction before merging.')
    matching = {}
    for p in phast_files:
        if p.stat().st_size <= 0:
            die(f'Empty PHAST output: {p}')
        if not (p.name.startswith('PHAST.') and p.name.endswith('.root')):
            die(f'Unexpected PHAST filename: {p.name}')
        tgeant_name = p.name[len('PHAST.'):-len('.root')] + '.tgeant.gz'
        if tgeant_name not in source_index:
            die(f'PHAST output has no corresponding TGEANT file: {p.name}')
        if tgeant_name in matching:
            die(f'Duplicate PHAST output for TGEANT source: {tgeant_name}')
        matching[tgeant_name] = p
    missing = sorted(set(source_index) - set(matching), key=natural_key)
    if missing:
        die(f'{len(missing)} TGEANT files lack corresponding PHAST output; e.g. {missing[:4]}')
    # Stable, natural order by original TGEANT batch naming.
    nominal_pairs = [(matching[p.name], source_index[p.name]) for p in source_files]
    outdir.mkdir(parents=True, exist_ok=True)
    chunksdir = outdir / '.merge_chunks'
    chunksdir.mkdir(parents=True, exist_ok=True)
    root_macro = chunksdir / 'validate_phast_roots.C'
    scan_macro = chunksdir / 'scan_phast_roots.C'
    atomic_write(root_macro, VALIDATE_MACRO)
    atomic_write(scan_macro, SCAN_MACRO)
    sources_manifest = outdir / 'phast_input_manifest.tsv'
    atomic_write(sources_manifest, ''.join(f'{p}\t{n}\n' for p, n in nominal_pairs))

    print('==============================================================', flush=True)
    print('PHAST merge preflight', flush=True)
    print('Input files:     ', len(nominal_pairs), flush=True)
    print('Nominal events:  ', f'{source_total:,}', flush=True)
    print('Output:          ', target, flush=True)
    print('Chunk size:      ', chunksize, flush=True)
    print('Delete TGEANT:   ', 'YES after merge' if delete_tgeant else 'NO', flush=True)
    print('==============================================================', flush=True)

    # CORAL may process 999 of 1000 events (as confirmed for batch 109).
    # Require one readable, nonempty PHAST tree per TGEANT file, but allow losses.
    print('Counting actual entries in all PHAST ROOT trees ...', flush=True)
    actual_manifest = outdir / 'phast_actual_entries.tsv'
    observed = root_scan(sources_manifest, tree_name, actual_manifest,
                         scan_macro, expected_files)
    if set(observed) != {str(p) for p, _ in nominal_pairs}:
        die('PHAST scan contained missing/unexpected file paths')
    pairs = [(p, observed[str(p)]) for p, _ in nominal_pairs]
    actual_total = sum(n for _, n in pairs)
    lost_total = source_total - actual_total
    if lost_total < 0:
        die('Actual PHAST entries exceed nominal generated events')
    loss_stats = {
        'nominal_events': source_total,
        'actual_phast_entries': actual_total,
        'missing_from_phast': lost_total,
        'files_below_nominal': sum(a < nom for (p, nom), (_, a)
                                   in zip(nominal_pairs, pairs)),
    }
    atomic_write(outdir / 'phast_entry_counts.json',
                 json.dumps(loss_stats, indent=2) + '\n')
    print('Actual PHAST entries:', f'{actual_total:,}', flush=True)
    print('Difference to nominal:', f'{lost_total:,}', flush=True)

    # Now validate EXACTLY the observed count for every ROOT file.
    root_validate(actual_manifest, tree_name, actual_total, root_macro)
    if target.exists() and use_existing == 1:
        print('Existing merged ROOT file detected: validating only.', flush=True)
        root_validate_one(target, tree_name, actual_total,
                          chunksdir / 'validate_final.tsv', root_macro)
    else:
        if target.exists() and not overwrite:
            die(f'Output already exists: {target}\n'
                'Set USE_EXISTING_MERGE=1 to validate or OVERWRITE=1 to replace.')
        chunk_outputs = []
        for i in range(0, len(pairs), chunksize):
            k = i // chunksize
            section = pairs[i:i + chunksize]
            fragment = ''.join(
                f'{p}\t{p.stat().st_size}\t{p.stat().st_mtime_ns}\t{nevt}\n'
                for p, nevt in section
            )
            sha = hashlib.sha256(fragment.encode()).hexdigest()[:16]
            chunk = chunksdir / f'chunk_{k:04d}_{sha}.root'
            count_expected = sum(n for _, n in section)
            chunk_manifest = chunksdir / f'chunk_{k:04d}_{sha}.validated'
            if chunk.exists() and chunk_manifest.exists():
                try:
                    root_validate_one(chunk, tree_name, count_expected,
                                      chunksdir / 'validate_chunk.tsv', root_macro)
                    print(f'Reusing verified chunk {k+1}', flush=True)
                except subprocess.CalledProcessError:
                    print(f'Invalid cached chunk {k+1}, rebuilding.', flush=True)
                    chunk.unlink(missing_ok=True)
                    chunk_manifest.unlink(missing_ok=True)
            if not (chunk.exists() and chunk_manifest.exists()):
                print(f'Chunk {k+1}/{(len(pairs)+chunksize-1)//chunksize}', flush=True)
                run_hadd(chunk, [p for p, _ in section])
                root_validate_one(chunk, tree_name, count_expected,
                                  chunksdir / 'validate_chunk.tsv', root_macro)
                atomic_write(chunk_manifest, f'entries={count_expected}\nsha={sha}\n')
            chunk_outputs.append(chunk)
        tmp = outdir / f'.{target.stem}.incomplete.{os.getpid()}.root'
        try:
            run_hadd(tmp, chunk_outputs)
            root_validate_one(tmp, tree_name, actual_total,
                              chunksdir / 'validate_final.tsv', root_macro)
            os.replace(tmp, target)
        finally:
            tmp.unlink(missing_ok=True)
        if os.environ.get('KEEP_CHUNKS', '0') == '0':
            for p in chunk_outputs:
                p.unlink(missing_ok=True)
            for p in chunksdir.glob('chunk_*.validated'):
                p.unlink(missing_ok=True)
    summary = {
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'phast_dir': str(Path(input_dir).resolve()),
        'tgeant_dir': str(Path(tgeant_dir).resolve()),
        'phast_files': len(pairs),
        'nominal_generated_events': source_total,
        'actual_phast_entries': actual_total,
        'difference_nominal_minus_phast': lost_total,
        'phast_tree_count_method': 'GetEntries() for each PHAST ROOT file',
        'tree_name': tree_name,
        'merged_root': str(target),
        'merged_size_bytes': target.stat().st_size,
        'source_manifest_sha256': digest,
        'tgeant_deleted': False,
        'tgeant_cleanup_state': 'not_requested',
        'note': 'Validated PHAST outputs only; does not archive the complete detector simulation or MC truth.',
    }
    summary_path = outdir / 'merge_summary.json'
    atomic_write(summary_path, json.dumps(summary, indent=2) + '\n')

    if delete_tgeant:
        # Deliberate opt-in. Do not touch original TGEANT data before ALL
        # 7,000 source trees and the final merged tree have been validated.
        root_validate_one(target, tree_name, actual_total,
                          chunksdir / 'validate_final_before_cleanup.tsv', root_macro)
        refreshed_paths, refreshed_index, refreshed_total, refreshed_digest = tgeant_index(
            tgeant_dir, tgeant_pattern)
        if (refreshed_digest != digest or refreshed_total != source_total or
                set(refreshed_index) != set(source_index)):
            die('TGEANT inputs changed during merge: refusing deletion')
        tgeant_base = Path(tgeant_dir).resolve()
        if target.resolve().is_relative_to(tgeant_base):
            die('Merged output lies within TGEANT source directory: refusing deletion')
        for p in refreshed_paths:
            if p.is_symlink() or not p.is_file() or p.parent.resolve() != tgeant_base:
                die(f'Unexpected TGEANT path (not a regular direct child): {p}')

        # A durable manifest is written BEFORE the first source file is removed.
        cleanup_manifest = outdir / 'tgeant_deletion_manifest.tsv'
        atomic_write(cleanup_manifest,
                     'filename\tbytes\tmtime_ns\tnominal_events\n' +
                     ''.join(f'{p.name}\t{p.stat().st_size}\t{p.stat().st_mtime_ns}\t'
                             f'{source_index[p.name]}\n' for p in refreshed_paths))
        summary['tgeant_cleanup_state'] = 'in_progress'
        summary['tgeant_cleanup_manifest'] = str(cleanup_manifest)
        atomic_write(summary_path, json.dumps(summary, indent=2) + '\n')
        print('==============================================================', flush=True)
        print('MERGE VALIDATED. Deleting exactly the matching TGEANT inputs:', flush=True)
        print('Source directory:', tgeant_base, flush=True)
        print('Number of files:', len(refreshed_paths), flush=True)
        print('Deletion manifest:', cleanup_manifest, flush=True)
        print('This deletion is irreversible.', flush=True)
        print('==============================================================', flush=True)

        # Record each completed removal; useful if an I/O failure interrupts cleanup.
        journal_path = outdir / 'tgeant_deletion_completed.txt'
        if journal_path.exists():
            die(f'Existing deletion journal found: {journal_path}; refusing automatic restart')
        deleted_bytes = 0
        with journal_path.open('x', encoding='utf-8') as journal:
            for i, p in enumerate(refreshed_paths, 1):
                size = p.stat().st_size
                p.unlink()
                deleted_bytes += size
                journal.write(p.name + '\n')
                if i % 100 == 0 or i == len(refreshed_paths):
                    journal.flush()
                    os.fsync(journal.fileno())
                    print(f'Deleted {i}/{len(refreshed_paths)} TGEANT files', flush=True)
        summary['tgeant_deleted'] = True
        summary['tgeant_cleanup_state'] = 'complete'
        summary['tgeant_deleted_files'] = len(refreshed_paths)
        summary['tgeant_deleted_bytes'] = deleted_bytes
        atomic_write(summary_path, json.dumps(summary, indent=2) + '\n')
        print('MERGE SUCCESSFUL; matching TGEANT source files deleted.', flush=True)
    else:
        print('MERGE SUCCESSFUL; source TGEANT files were NOT deleted.', flush=True)

    print('==============================================================', flush=True)
    print('Merged ROOT: ', target, flush=True)
    print('Tree:        ', tree_name, flush=True)
    print('Entries:     ', f'{actual_total:,}', flush=True)
    print('Size:        ', f'{target.stat().st_size / 1024**2:.1f} MiB', flush=True)
    print('Input list:  ', sources_manifest, flush=True)
    print('Summary:     ', outdir / 'merge_summary.json', flush=True)
    print('==============================================================', flush=True)


if __name__ == '__main__':
    try:
        if len(sys.argv) != 2 or sys.argv[1] not in ('count', 'merge'):
            die('Usage: coralphast_merge_tools.py count|merge')
        {'count': count, 'merge': merge}[sys.argv[1]]()
    except (RuntimeError, OSError, subprocess.CalledProcessError, json.JSONDecodeError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(1)
