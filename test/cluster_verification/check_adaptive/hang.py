#!/usr/bin/env python3
"""A binary that will not stop, for the timeout path (plan stage 1, Task 6).

  python hang.py ignore-term      # ignores SIGTERM: only SIGKILL ends it
  python hang.py obey             # stops on SIGTERM
  python hang.py child            # spawns a child and waits, so the whole bundle
                                  # has to be signalled, not just this process

Prints a line first so run.log is never empty, then waits. Nothing here is part
of the port; it exists so the timeout can be tested without a real run.
"""
import os
import signal
import subprocess
import sys
import time


def main(argv):
    mode = argv[0] if argv else 'obey'
    if mode == 'ignore-term':
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
    print(f'hanging ({mode})', flush=True)
    if mode == 'child':
        subprocess.Popen([sys.executable, __file__, 'ignore-term'])
    while True:
        time.sleep(0.05)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
