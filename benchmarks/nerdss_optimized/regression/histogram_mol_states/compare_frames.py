#!/usr/bin/env python3
"""Check that two histogram_complexes_time.dat files agree at every time they share.

Usage: compare_frames.py FULL CONTINUATION

Lines within a frame are an unordered set (the histogram is written from an
unordered_map), so frames are compared as multisets of (composition, count).
Exit 1 if any shared time differs or no time is shared.
"""
import sys
from check_histogram import read_histogram, same_time

full = read_histogram(sys.argv[1])
cont = read_histogram(sys.argv[2])
shared = 0
bad = 0
for t_cont, frame_cont in cont:
    frame_full = next((f for t, f in full if same_time(t, t_cont)), None)
    if frame_full is None:
        continue
    shared += 1
    if frame_full != frame_cont:
        bad += 1
        print(f'FAIL: t={t_cont}: continuation frame differs from the full run', file=sys.stderr)
if shared == 0:
    sys.exit('FAIL: the two files share no frame time')
print(f'shared frames: {shared}, differing: {bad}')
sys.exit(1 if bad else 0)
