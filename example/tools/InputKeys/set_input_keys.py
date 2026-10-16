#!/usr/bin/env python3
"""Write a copy of a CCEX input with some keys set (A63: the values the removed command-line options gave).

    python3 set_input_keys.py IN OUT KEY=VALUE ...

KEY is a key or a dotted path into a block (configuration.defect_axis_file); the block has to exist in IN, the key itself may be new.
Keys are matched without regard to case, as CCEX reads them. VALUE is read as JSON when it is JSON (50, true,
["bath_1"]) and as a string otherwise (Random3). OUT is plain JSON (the comments of IN are not kept) and may not be IN.
File paths in an input are relative to the working directory of the run, not to the input, so OUT may be anywhere.
Any error: a message and exit 2.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from convert_input_keys import strip_comments   # noqa: E402


def fail(message):
    print(f"set_input_keys.py: {message}", file=sys.stderr)
    sys.exit(2)


def same_key(block, name):
    """The key of block that is name without regard to case, or name itself when there is none."""
    return next((k for k in block if k.lower() == name.lower()), name)


def set_key(root, path, value):
    parts = path.split('.')
    block = root
    for i, part in enumerate(parts[:-1]):
        key = same_key(block, part)
        if not isinstance(block.get(key), dict):
            fail(f"{'.'.join(parts[:i + 1])} is not a block of the input (needed for {path})")
        block = block[key]
    block[same_key(block, parts[-1])] = value


def main(argv):
    if len(argv) < 2:
        fail("usage: set_input_keys.py IN OUT KEY=VALUE ...")
    src, out, settings = argv[0], argv[1], argv[2:]
    if os.path.abspath(out) == os.path.abspath(src):
        fail(f"OUT is IN ({src}); write the copy elsewhere")
    try:
        with open(src) as f:
            root = json.loads(strip_comments(f.read()))
    except (OSError, ValueError) as e:
        fail(f"cannot read {src}: {e}")
    if not isinstance(root, dict):
        fail(f"{src} is not a JSON object")
    for setting in settings:
        key, eq, text = setting.partition('=')
        if not eq or not key:
            fail(f"{setting!r} is not KEY=VALUE")
        try:
            value = json.loads(text)
        except ValueError:
            value = text
        set_key(root, key, value)
    try:
        with open(out, 'w') as f:
            json.dump(root, f, indent=4)
            f.write('\n')
    except OSError as e:
        fail(f"cannot write {out}: {e}")


if __name__ == '__main__':
    main(sys.argv[1:])
