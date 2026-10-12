// Running a real process, and ending one that will not end (plan stage 1, Task 6).
#include "check.h"
#include "mkdirs.h"
#include "base_run.h"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {
std::string slurp(const std::string& path) {
    std::ifstream file(path);
    std::stringstream text;
    text << file.rdbuf();
    return text.str();
}

double seconds_of(void (*body)()) {
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    body();
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}
}   // namespace

int main(int argc, char** argv) {
    const std::string python = argc > 1 ? argv[1] : "python3";
    testing::mkdirs("out/native");

    // A process that finishes on its own: the log holds what it printed.
    adaptive::NativeResult ok = adaptive::run_native(
        {python, "-c", "print('hello'); import sys; sys.stderr.write('and stderr\\n')"},
        ".", "out/native/plain.log", 60.);
    CHECK_EQ(ok.returncode, 0);
    CHECK(!ok.timed_out);
    CHECK(slurp("out/native/plain.log").find("hello") != std::string::npos);
    CHECK(slurp("out/native/plain.log").find("and stderr") != std::string::npos);

    // A non-zero exit is reported, not thrown.
    adaptive::NativeResult failed = adaptive::run_native({python, "-c", "raise SystemExit(3)"},
                                                     ".", "out/native/fail.log", 60.);
    CHECK_EQ(failed.returncode, 3);
    CHECK(!failed.timed_out);

    // The environment the run gets.
    bool pinned = false, fabrics = false;
    for (const std::string& entry : adaptive::native_environment()) {
        if (entry == "OMP_NUM_THREADS=1") pinned = true;
        if (entry.compare(0, 16, "I_MPI_FABRICS=sh") == 0) fabrics = true;
    }
    CHECK(pinned);
    CHECK(fabrics);

    // One that ignores SIGTERM has to be killed, and the whole bundle with it.
    const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    adaptive::NativeResult hung = adaptive::run_native(
        {python, "../hang.py", "ignore-term"}, ".", "out/native/hang.log", 1.);
    const double took = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    CHECK(hung.timed_out);
    CHECK(hung.returncode == -9);                 // SIGKILL, the way Python reports it
    CHECK(took >= 6.);                            // one second, then five of grace
    CHECK(took < 20.);
    CHECK(slurp("out/native/hang.log").find("hanging") != std::string::npos);

    // One that obeys SIGTERM stops in the grace period rather than after it.
    const std::chrono::steady_clock::time_point began = std::chrono::steady_clock::now();
    adaptive::NativeResult polite = adaptive::run_native(
        {python, "../hang.py", "obey"}, ".", "out/native/obey.log", 1.);
    const double waited = std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
    CHECK(polite.timed_out);
    CHECK_EQ(polite.returncode, -15);             // SIGTERM
    CHECK(waited < 5.);

    (void)seconds_of;
    std::cout << "ran four processes, two of them ended by the timeout\n";
    return check::report("test_native");
}
