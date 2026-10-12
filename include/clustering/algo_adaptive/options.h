// The command line, from the 22 options of options_table.h (42 until the four removed policies took 20 with them, 2026-10-10).
//
// The table is generated from adaptive_refine.argument_parser(), so the two stay
// in step. What is reproduced is the set of options, their types, their defaults
// and which are required - not argparse's usage text, which ends the run before
// any decision is made (declared exception).
#pragma once
#include <map>
#include <string>
#include <vector>
#include "base_json.h"

namespace adaptive {

// Parsed options, by their Python name (dest). Values carry their type: a Path and
// a Text are strings, an Int a big integer, a Float a double, a Flag a bool.
Value parse_options(const std::vector<std::string>& argv, const std::string& working_directory);
}   // namespace adaptive
