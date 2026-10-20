#!/usr/bin/env python3
"""Run a command under a pty and prefix every output line with seconds since launch.
Usage: tsrun.py <logfile> <cmd> [args...]   (prints the child's pid on stderr)"""
import os, pty, sys, time, select, signal
log = open(sys.argv[1], 'w', buffering=1)
cmd = sys.argv[2:]
t0 = time.monotonic()
log.write(f'# t0_monotonic {t0:.6f}\n')
pid, fd = pty.fork()
if pid == 0:
    os.execvp(cmd[0], cmd)
sys.stderr.write(f"child pid {pid}\n"); sys.stderr.flush()
open(sys.argv[1] + '.pid', 'w').write(str(pid))
buf = b''
while True:
    try:
        r, _, _ = select.select([fd], [], [], 1.0)
    except InterruptedError:
        continue
    if fd in r:
        try:
            data = os.read(fd, 65536)
        except OSError:
            data = b''
        if not data:
            break
        buf += data
        while b'\n' in buf:
            line, buf = buf.split(b'\n', 1)
            log.write(f"{time.monotonic()-t0:12.6f} {line.decode(errors='replace').rstrip()}\n")
if buf:
    log.write(f"{time.monotonic()-t0:12.6f} {buf.decode(errors='replace')}\n")
_, status = os.waitpid(pid, 0)
log.write(f"{time.monotonic()-t0:12.6f} [exit status {status}]\n")
