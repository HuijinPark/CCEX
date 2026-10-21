/**
 * @file options_table.h
 * @brief Adaptive refinement: the table of every option, its type and default.
 */
// Every option the refinement takes, on the command line (--name) or as a key of
// the input file's cluster.adaptive block (name with _ for -; "workdir" for output). This table is the one
// place an option is declared: add a row here to add an option. 27 options (20 that
// only the removed policies read went with them on 2026-10-10; the three regroup ones
// came on 2026-10-11, and regroup-factors and regroup-whole later that day).
#pragma once
#include <cstddef>

namespace adaptive {
namespace options {
enum class Kind { Path, Text, Int, Float, Flag };

struct Option {
    const char* name;        // without the leading --
    const char* dest;        // the name the policy code uses
    Kind kind;
    bool required;
    const char* default_text;
    bool no_default;         // Python's None: absent rather than empty
    bool int_default;        // the default is an int even where the option is a float
    const char* choices;     // "" when anything goes, else the allowed values, '|' apart
};

inline const Option kOptions[] = {
    {"input", "input", Kind::Path, true, "", true, false, ""},
    {"binary", "binary", Kind::Path, true, "", true, false, ""},
    {"output", "output", Kind::Path, true, "", true, false, ""},
    {"cwd", "cwd", Kind::Path, false, "<cwd>", false, false, ""},
    {"mpiexec", "mpiexec", Kind::Text, false, "mpirun", false, false, ""},
    {"ranks", "ranks", Kind::Int, false, "1", false, true, ""},
    {"atol", "atol", Kind::Float, true, "", true, false, ""},
    {"rtol", "rtol", Kind::Float, true, "", true, false, ""},
    {"max-rounds", "max_rounds", Kind::Int, false, "8", false, true, ""},
    {"max-patches", "max_patches", Kind::Int, false, "64", false, true, ""},
    {"max-probes", "max_probes", Kind::Int, false, "4096", false, true, ""},
    {"max-candidates", "max_candidates", Kind::Int, false, "", true, false, ""},
    {"neighbors", "neighbors", Kind::Int, false, "8", false, true, ""},
    {"cache-mb", "cache_mb", Kind::Int, false, "256", false, true, ""},
    {"timeout", "timeout", Kind::Float, false, "600", false, true, ""},
    {"raw-cache", "raw_cache", Kind::Path, false, "", true, false, ""},
    {"probe-batch-clusters", "probe_batch_clusters", Kind::Int, false, "0", false, true, ""},
    {"merge-priority", "merge_priority", Kind::Text, false, "gain-per-cost", false, false, "gain-per-cost|worst-error"},
    // Only merge since 2026-10-10; adaptive_main says why a removed name is refused.
    {"policy", "policy", Kind::Text, false, "merge", false, false, "merge"},
    {"max-total-work", "max_total_work", Kind::Float, false, "8000000000000.0", false, false, ""},
    {"prefer-physical", "prefer_physical", Kind::Flag, false, "", false, false, ""},
    {"hub-candidates", "hub_candidates", Kind::Flag, false, "", false, false, ""},
    // Method D ((internal notes)). Off, these three are left out
    // of the report's settings, so a run without them records what it always did.
    {"regroup", "regroup", Kind::Flag, false, "", false, false, ""},
    {"regroup-max-regions", "regroup_max_regions", Kind::Int, false, "64", false, true, ""},
    {"regroup-max-candidates", "regroup_max_candidates", Kind::Int, false, "256", false, true, ""},
    // Where D's culprit factors and its whole curve for condition (c) come from
    // (adaptive-regroup-stage3.md, -stage3b.md). At "probes" they are left out of the
    // report's settings like the three above, so every run that does not ask records
    // what it did before they existed.
    {"regroup-factors", "regroup_factors", Kind::Text, false, "probes", false, false, "probes|factors"},
    {"regroup-whole", "regroup_whole", Kind::Text, false, "probes", false, false, "probes|except"},
};

inline std::size_t count() { return sizeof kOptions/sizeof *kOptions; }
}   // namespace options
}   // namespace adaptive
