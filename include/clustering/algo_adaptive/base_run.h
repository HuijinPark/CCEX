// Running CCEX once and keeping what came back: the launcher, the clock, SHA-256,
// the input manifest that identifies a run, and the on-disk cache of raw probe curves.
#pragma once
#include <functional>
#include <map>
#include <string>
#include <vector>
#include "base_json.h"

// ----------------------------------------------------------------------------
// SHA-256 (plan 0B-4). The cache keys, the source manifests and the binary
// fingerprint are all this hash, so it has to agree with Python's hashlib.
// ----------------------------------------------------------------------------
namespace adaptive {
std::string sha256(const std::string& data);
std::string sha256_file(const std::string& path);   // raises like Python's open() would
// The CCEX binary's SHA-256, kept for the life of the process and computed again only when the file's metadata (size,
// mtime and ctime in ns, inode, device; stat on the path each time, so a replaced file is seen) differ from the last
// time (review 2026-10-08: hashing the 123 MB binary whole was 99 % of the driver's own time). A same-size overwrite
// that also keeps mtime and ctime is not seen - nor one inside the same tick of a coarse file-system clock (a few ms on
// WSL, measured), so a binary must not be rebuilt within milliseconds of a check; input files keep sha256_file. A file that changes while it is hashed is
// a RuntimeError. binary_hash_reads() counts the whole-file reads, for the tests.
std::string binary_sha256(const std::string& path);
long long binary_hash_reads();
}   // namespace adaptive

// ----------------------------------------------------------------------------
// The clock the decision code reads.
//
// Time is part of the decision, not just a measurement: a budget's remaining
// seconds gate whether another order is tried, and the timeout handed to each
// native call is the validation budget minus what has been spent. So the
// decision code never calls the clock directly - it is handed one, which is what
// lets a test hand it a clock of its own.
// ----------------------------------------------------------------------------
namespace adaptive {
class Clock {
 public:
    virtual ~Clock() = default;
    virtual double monotonic() = 0;
    virtual double perf_counter() = 0;
};

class SystemClock : public Clock {
 public:
    double monotonic() override;
    double perf_counter() override;
};
}   // namespace adaptive

// ----------------------------------------------------------------------------
// Running the CCEX binary (plan stage 1, Task 6; adaptive_refine.py).
//
// The process is put in a session of its own so that a timeout can end the whole
// bundle rather than just the launcher: mpiexec's children would otherwise outlive
// it. SIGTERM first, five seconds, then SIGKILL - the same two steps Python takes.
// ----------------------------------------------------------------------------
namespace adaptive {
struct NativeResult {
    int returncode = 0;
    bool timed_out = false;
};

// Runs `command` with `cwd` as its directory, everything it prints going to
// `log_path`. `timeout` is in seconds; a timeout ends the process bundle.
NativeResult run_native(const std::vector<std::string>& command, const std::string& cwd,
                        const std::string& log_path, double timeout);

// Run the calculation in this process instead of spawning one (plan stage 8).
//
// Set by bin/main.out, which has already read the bath once and can answer a
// step without paying for that again; bin/ccex-adaptive leaves it unset and keeps
// spawning. The runner is handed the same command run_native would have run, so
// the step directory and its input.json reach it unchanged - the file protocol
// the recordings are written in does not move.
using LocalRunner = std::function<NativeResult(const std::vector<std::string>& command,
                                               const std::string& cwd,
                                               const std::string& log_path)>;

// Install a runner for the steps, or clear it by passing nothing. Used by
// main.out --adaptive, which already is the binary a step would launch.
void set_local_runner(LocalRunner runner);

// The environment a run gets: the parent's, with the thread counts pinned to one
// and I_MPI_FABRICS defaulted - defaulted, not forced, so a machine that needs
// something else keeps it.
std::vector<std::string> native_environment();
}   // namespace adaptive

// The refinement tool's body, callable from another main() (plan stage 8).
// Defined in src/clustering/algo_adaptive/adaptive_main.cpp, at global scope as its own main() was.
int adaptive_main(const std::vector<std::string>& arguments);

// ----------------------------------------------------------------------------
// The input files a run depends on, and the context that identifies it
// (plan stage 1, Task 5; adaptive_probe_cache.py).
//
// 0B-4 ported the cache itself but not these: the tests wrote the context by
// hand. The cache key is the SHA-256 of this context's text, so what is in it and
// how its paths are spelled decides which runs share probes.
// ----------------------------------------------------------------------------
namespace adaptive {
// str(pathlib.Path(cwd, name).resolve()): absolute, symlinks followed, "." and
// ".." taken out - and NOT the same as the plain normalisation the CLI does to a
// path option, which leaves links alone. A component that does not exist is kept
// as written, the way resolve(strict=False) keeps it.
std::string resolved(const std::string& cwd, const std::string& name);

// Every input file the config names, with its hash, or null where the file is
// not there.
Value input_manifest(const Value& cfg, const std::string& cwd, const std::vector<int>& states);

// The manifest as one comparable string. Two drivers check the cache against it,
// and they must compare it the same way.
std::string manifest_text(const Value& files);

// What makes two runs share a cache: the config minus the options that only
// choose outputs or bundleing, the directory, the states, the ranks, the binary,
// and the input files.
Value identity(const Value& cfg, const std::string& cwd, const std::vector<int>& states,
               const std::string& binary_sha256, long long ranks);
}   // namespace adaptive

// ----------------------------------------------------------------------------
// The raw probe cache, a file shared with the Python implementation (plan 0B-4).
//
// adaptive_probe_cache.py writes it and reads it; so does this. That is
// only possible if both sides produce the same text, because the context key and
// every payload's integrity check are SHA-256 hashes of the JSON itself.
// ----------------------------------------------------------------------------
struct sqlite3;

namespace adaptive {

using Sites = std::vector<int>;                       // one probe's site tuple
using Curve = std::vector<std::pair<double, double>>; // one state's curve, as complex pairs
using Curves = std::vector<Curve>;                    // one per state

class ProbeCache {
 public:
    ProbeCache(const std::string& path, const Value& context);
    ~ProbeCache();

    std::map<Sites, Curves> get(const std::vector<Sites>& wanted);
    void put(const std::map<Sites, Curves>& values);
    const std::string& key() const { return key_; }   // the context hash
    const std::string& path() const { return path_; }
    void close();

 private:
    sqlite3* db_ = nullptr;
    std::string path_;
    std::string key_;
    std::size_t states_ = 0;
    long long steps_ = 0;
};
}   // namespace adaptive
