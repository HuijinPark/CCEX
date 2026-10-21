// Make a directory path, one segment at a time (plan 0B-1, Task 6). Test only.
// The saving function deliberately does not do this: Python's save() raises
// FileNotFoundError when the parent is missing, and the port reproduces that.
#pragma once
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

namespace testing {
inline void mkdirs(const std::string& path) {
    for (std::size_t i = 1; i <= path.size(); ++i)
        if (i == path.size() || path[i] == '/')
            mkdir(path.substr(0, i).c_str(), 0777);     // EEXIST is the normal case
}
}   // namespace testing
