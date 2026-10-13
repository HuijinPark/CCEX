#!/bin/sh
# What a change to the refinement code (src/clustering/algo_adaptive) has to pass.
#   verify.sh unit                 # build the driver and its unit tests, run them
#   verify.sh regress              # the policy against its recorded answers, run for real
#   verify.sh inprocess            # the driver alone, inside bin/main.out, and from an "adaptive" block; the probe
#                                  # step on 1 and 3 ranks with its memory budget (probe_rotation.py)
#   verify.sh sanitize             # unit and regress once more, built with ASan and UBSan
#   verify.sh fuzz [SEED [COUNT]]  # random JSON documents through the reader and the writer
# regress, inprocess and sanitize need the built CCEX (bin/main.out) and mpirun.
#
# here is absolute, and every path handed to a program is absolute: building them
# relative to the test directory broke as soon as the script was called from
# somewhere else, which is exactly how the servers call it.
set -e
here=$(cd "$(dirname "$0")" && pwd)
# One run at a time per folder: a second run would delete the first one's work files (D29).
# A script this run starts (check_refusals runs legacy.sh) inherits the lock and does not take it again.
if [ "${CCEX_RUN_LOCK:-}" != "$here" ]; then
    exec 9>"$here/.run.lock"
    if command -v flock >/dev/null 2>&1 && ! flock -n 9; then
        echo "ERROR: another check is running in $here (it would delete that run's work files); wait for it or use another checkout" >&2; exit 1
    fi
    CCEX_RUN_LOCK="$here"; export CCEX_RUN_LOCK
fi
py=${PYTHON:-python3}

# The unit tests run from cpp/ and leave what they read and wrote under cpp/out.
unit() {    # unit BIN [make arguments]
    bin=$1; shift
    rm -rf "$here/cpp/out" "$here/cpp/$bin"
    make -C "$here/cpp" BIN="$bin" "$@" all
    # Every power of two and both its neighbours: a random corpus never lands on
    # these, and this is where the float text once went wrong.
    (cd "$here/cpp" && "./$bin/test_float_text" "$here/reference/semantics/boundaries.tsv")
    for f in "$here"/regression/inputs/*.json; do
        (cd "$here/cpp" && "./$bin/test_json_config" --file "$f" > /dev/null)
    done
    # What the writer wrote and what the reader read, against the stored expectations.
    "$py" "$here/json_roundtrip.py" "$here/cpp/out/written.json"
    "$py" "$here/json_expect.py" cases "$here/cpp/out/read" "$here/reference/json-cases"
    (cd "$here" && "$py" -m unittest test_adaptive_tools)
}

case "$1" in
unit)
    unit bin
    ;;
regress)
    "$py" "$here/regress.py" "$here/cpp/out/regress"
    ;;
inprocess)
    "$py" "$here/inprocess_compare.py" "$here/cpp/out/inprocess"
    # the probe step across ranks and states, with its memory budget (A1)
    "$py" "$here/probe_rotation.py" "$here/cpp/out/probe-rotation"
    ;;
sanitize)
    # A sanitizer stops the program it finds something in, so a finding in the
    # driver shows up as an exit status the golden files do not hold.
    unit bin-asan CXXFLAGS="-I../../../../external/sqlite -std=c++17 -O1 -g -Wall -Wextra -MMD -MP -I../../../../include/clustering/algo_adaptive -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=undefined" \
        SQLITE_SANITIZE="-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=undefined"
    CCEX_ADAPTIVE_BIN="$here/cpp/bin-asan/ccex-adaptive" "$py" "$here/regress.py" "$here/cpp/out/regress-asan"
    ;;
fuzz)
    # The defects found in the JSON reader were all inputs nobody had listed, so
    # new ones are generated for each seed and python3's json is the judge, of
    # reading, of writing and of which configs are refused.
    seed=${2:-1}
    count=${3:-20000}
    corpus=$(mktemp -d)
    trap 'rm -rf "$corpus"' EXIT
    "$py" "$here/fuzz_json.py" "$corpus" "$count" "$seed"
    rm -rf "$here/cpp/out/all" "$here/cpp/out/back"
    (cd "$here/cpp" && ./bin/test_json_read --all "$corpus")
    "$py" "$here/json_expect.py" tree "$here/cpp/out/all" "$corpus"
    (cd "$here/cpp" && ./bin/test_json_read --write-back "$corpus")
    "$py" "$here/json_expect.py" back "$here/cpp/out/back" "$corpus"
    (cd "$here/cpp" && ./bin/test_json_config --verdicts "$corpus" > "$corpus/verdicts.tsv")
    "$py" "$here/config_verdicts.py" "$corpus" "$corpus/verdicts.tsv"
    ;;
*)
    echo "usage: verify.sh unit | regress | inprocess | sanitize | fuzz [SEED [COUNT]]" >&2
    exit 2
    ;;
esac
