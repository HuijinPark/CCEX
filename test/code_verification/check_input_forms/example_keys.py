#!/usr/bin/env python3
"""example_keys.py EXAMPLE REPO: the example input `main.out -e` writes names every input key (A63).

The readers of the input (src/input/option.cpp, and src/main.cpp for the adaptive and auto_convergence blocks) look
their keys up by name, many of them only when another key asks for them, so no single run shows them all. The contract
is checked on the source instead, and is closed: every identifier-like string literal of option.cpp has to appear in
the example, so a key added to a reader and not to the example fails here.

1. The example without its comments is valid JSON (its active part is the input that runs).
2. A literal looked up as a key - cJSON_GetObjectItem(x, "key") or readBlock(x, "key") in option.cpp and main.cpp - is
   written in the example as "key": (on an active or a commented line).
3. Every other identifier-like literal of option.cpp (keys read through helpers, key tables, allowed values, units)
   appears in the example at least as a word, except the refused old names (taken from option.cpp's own table of
   renamed keys) and the internal names listed in NOT_INPUT below.
4. The keys of the two blocks read by other parts: every option of the refinement (options_table.h kOptions, '-' as
   '_', without the ones the run decides itself) and every key of the automatic convergence (autoconv/runner.py KEYS,
   TOP_KEYS, SECTIONS, SPLIT_KEYS) is written as "key":.
5. Conversely, every "key": in the example is a name one of those sources knows (a typo in the example fails).
Exit 0 when all hold, 1 with the missing names otherwise.
"""
import ast, json, re, sys

example_path, repo = sys.argv[1:3]
sys.path.insert(0, f"{repo}/example/tools/InputKeys")
from convert_input_keys import strip_comments   # the CCEX comment rules ('#', '!' outside strings)

text = open(example_path).read()
problems = []
try:
    json.loads(strip_comments(text))
except ValueError as e:
    problems.append(f"the example without comments is not valid JSON: {e}")
written_keys = set(re.findall(r'"([A-Za-z_][A-Za-z0-9_]*)"\s*:', text))
words = set(re.findall(r'[A-Za-z_][A-Za-z0-9_]*', text))

def code_of(path):
    """The source without // comment lines (a literal there is no code)."""
    return '\n'.join(l for l in open(path).read().splitlines() if not l.lstrip().startswith('//'))

option = code_of(f"{repo}/src/input/option.cpp")
main = code_of(f"{repo}/src/main.cpp")
IDENT = r'"([A-Za-z_][A-Za-z0-9_]*)"'
looked_up = set()
for src in (option, main):
    looked_up |= set(re.findall(r'(?:cJSON_GetObjectItem|readBlock)\(\s*[^,()]+,\s*(?:\(char\s*\*\)\s*)?' + IDENT, src))
literals = set(re.findall(IDENT, option))
# the refused old names: option.cpp's table of renamed keys, {"block", "key", "is now ..." / "was removed ..."}
old = set(re.findall(r'\{\s*"[^"]*",\s*' + IDENT + r',\s*"(?:is now|was removed)', option))
# names in option.cpp that are not input keys or values: the internal method names the readers translate to (and the
# lower-case forms they compare), labels of the defect readers' messages, a section label of a message, the fopen mode,
# the prefix of the default qubit names, the refused names outside the table above (the flat gpu_* keys, the old
# "Refine" block, save_nodiv, magnetic_field_axis), and rdip_ladder - read by CCEX, but an internal key the automatic
# convergence writes into its runs' inputs (A65), not shown in the example
NOT_INPUT = {
    "pCCE", "pgCCE", "pmeCCE", "pmegCCE", "ccce", "cce", "gcce", "mecce", "megcce", "pgcce", "pmecce", "pmegcce",
    "cJSON_ReadDefectInfo_IntCharDouble", "cJSON_ReadDefectInfo_IntCharMatrixXcd1d", "cJSON_ReadDefectInfo_IntCharMatrixXcd2d",
    "r", "q", "rdip_ladder",
    "gpu_min_dim", "gpu_batch_clusters", "gpu_batch_steps", "gpu_cache_mb", "gpu_cluster_batch_mb", "refine", "save_nodiv",
    "magnetic_field_axis",
}
stale = sorted(n for n in NOT_INPUT if n not in literals or n in words)
if stale:
    problems.append("NOT_INPUT names that option.cpp no longer has, or that the example shows now (drop them here): "
                    + ", ".join(stale))

missing_keys = sorted(k for k in looked_up if k not in written_keys and k not in old and k not in NOT_INPUT)
if missing_keys:
    problems.append('looked up by the readers, not written as "key": in the example: ' + ', '.join(missing_keys))
missing_words = sorted(l for l in literals - looked_up if l not in words and l not in old and l not in NOT_INPUT)
if missing_words:
    problems.append('names of option.cpp missing from the example: ' + ', '.join(missing_words))

# the refinement's options (the input block's keys are the option names with '_'); input, binary and ranks are the run's
table = open(f"{repo}/include/clustering/algo_adaptive/options_table.h").read()
adaptive = {n.replace('-', '_') for n in re.findall(r'\{"([a-z-]+)",\s*"[a-z_]+",\s*Kind::', table)}
fixed = re.search(r'static const char\* fixed\[\] = \{([^}]*)\}', main)   # adaptiveArgsFromInput: decided by the run
decided = set(re.findall(r'"([a-z_]+)"', fixed.group(1))) if fixed else set()
if not fixed: problems.append("src/main.cpp: the list of refinement options the run decides (fixed[]) was not found")
adaptive_keys = (adaptive - decided) | {"work_file"}
# the automatic convergence's keys (the block is read by its runner, not by CCEX): the top-level keys, the sections and
# the keys of each section as the input writes them (runner.py TOP_KEYS, SECTIONS, SPLIT_KEYS; its FILE_KEYS are the CCEX
# configuration block's keys since A65, read by option.cpp, not keys of the block)
tree = ast.parse(open(f"{repo}/example/tools/AutoConvergence/autoconv/runner.py").read())
tables = {n.targets[0].id: n.value for n in tree.body
          if isinstance(n, ast.Assign) and len(n.targets) == 1 and isinstance(n.targets[0], ast.Name)}
def strings(node): return {c.value for c in ast.walk(node) if isinstance(c, ast.Constant) and isinstance(c.value, str)}
def call_keywords(node): return {k.arg for k in node.keywords} if isinstance(node, ast.Call) else set()
auto = {k.value for k in tables["SECTIONS"].keys}
for section in tables["SECTIONS"].values: auto |= call_keywords(section)
auto |= call_keywords(tables["SPLIT_KEYS"]) | strings(tables["TOP_KEYS"])
for group, keys in (("adaptive", adaptive_keys), ("auto_convergence", auto)):
    miss = sorted(k for k in keys if k not in written_keys)
    if miss: problems.append(f'{group} keys not written as "key": in the example: ' + ', '.join(miss))

known = literals | set(re.findall(IDENT, main)) | adaptive_keys | auto
unknown = sorted(k for k in written_keys if k not in known)
if unknown:
    problems.append('"key": in the example that no reader knows: ' + ', '.join(unknown))

print(f"example keys: {len(written_keys)} written; {len(looked_up)} looked up by the readers, {len(literals)} names in "
      f"option.cpp ({len(old)} refused old names), {len(adaptive_keys)} adaptive and {len(auto)} auto_convergence keys")
for p in problems: print("  " + p)
sys.exit(1 if problems else 0)
