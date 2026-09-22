#!/usr/bin/env python3
"""Compare the DATA time series of an uninterrupted run and of its restart.

    compare_series.py <uninterrupted run's DATA dir> <restart's DATA dir>

Every record the restart wrote must equal the uninterrupted run's record with
the same time label.  Line-oriented files are keyed by their first field, block
files by their 'time' header line.  Prints one line per file, and exits
non-zero if any common record differs or a file has no records in common.

This covers what restart.dat does not carry, such as the copy numbers that
init_counterCopyNums() recounts at a restart: the restart writes its first
record right after setting up, so that record is compared against the
uninterrupted run's record for the same step.
"""
import os
import re
import sys

LINE_FILES = ['copy_numbers_time.dat', 'mono_dimer_time.dat', 'bound_pair_time.dat']
BLOCK_FILES = {
    'histogram_complexes_time.dat': re.compile(r'^Time \(s\): (\S+)'),
    'event_counters_time.dat': re.compile(r'^time \(s\): (\S+)'),
    'transition_matrix_time.dat': re.compile(r'^time: (\S+)'),
}


def key(label):
    # Both runs print the time with the stream's default precision, but they
    # compute it by different expressions; round so 1e-4 and 0.0001 agree.
    return round(float(label), 12)


def line_records(path):
    recs = {}
    with open(path, errors='replace') as f:
        lines = f.read().splitlines()
    for line in lines[1:]:
        fields = re.split(r'[,\t ]+', line.strip())
        if not fields or not fields[0]:
            continue
        try:
            k = key(fields[0])
        except ValueError:
            continue
        recs.setdefault(k, []).append(fields[1:])
    return recs


def block_records(path, header):
    recs = {}
    cur = None
    with open(path, errors='replace') as f:
        for line in f.read().splitlines():
            m = header.match(line)
            if m:
                cur = key(m.group(1))
                recs.setdefault(cur, []).append([])
            elif cur is not None:
                recs[cur][-1].append(line.rstrip())
    return recs


def compare(name, full, cont):
    common = sorted(set(full) & set(cont))
    bad = [k for k in common if full[k][-1] != cont[k][-1]]
    status = 'ok' if common and not bad else ('NO COMMON RECORDS' if not common else 'DIFF')
    msg = '%-30s %-17s %4d common, %d differ' % (name, status, len(common), len(bad))
    if bad:
        msg += '; first at t=%g' % bad[0]
    print(msg)
    return status == 'ok'


def main():
    fdir, cdir = sys.argv[1], sys.argv[2]
    good = True
    for name in LINE_FILES:
        p, q = os.path.join(fdir, name), os.path.join(cdir, name)
        if os.path.exists(p) and os.path.exists(q):
            good &= compare(name, line_records(p), line_records(q))
    for name, header in BLOCK_FILES.items():
        p, q = os.path.join(fdir, name), os.path.join(cdir, name)
        if os.path.exists(p) and os.path.exists(q):
            good &= compare(name, block_records(p, header), block_records(q, header))
    sys.exit(0 if good else 1)


main()
