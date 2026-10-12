// binary_sha256 (review 2026-10-08): the binary's hash is read once per process and again only when its metadata change.
#include "check.h"
#include "base_run.h"
#include <cstdio>
#include <fcntl.h>
#include <string>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

namespace {
void write(const std::string& path, const std::string& text) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    std::fwrite(text.data(), 1, text.size(), f);
    std::fclose(f);
}

void set_mtime(const std::string& path, long seconds) {
    struct timespec times[2];
    times[0].tv_sec = seconds; times[0].tv_nsec = 0;
    times[1].tv_sec = seconds; times[1].tv_nsec = 0;
    utimensat(AT_FDCWD, path.c_str(), times, 0);
}
}  // namespace

int main() {
    char dir[] = "/tmp/ccex-binary-hash-XXXXXX";
    CHECK(mkdtemp(dir) != nullptr);
    const std::string bin = std::string(dir) + "/main.out";
    write(bin, "first binary");
    const long long start = adaptive::binary_hash_reads();

    // unchanged: one whole read, then the kept hash
    const std::string first = adaptive::binary_sha256(bin);
    CHECK_EQ(first, adaptive::sha256_file(bin));
    for (int i = 0; i < 5; ++i) CHECK_EQ(adaptive::binary_sha256(bin), first);
    CHECK_EQ(adaptive::binary_hash_reads() - start, 1LL);

    // touched, same content: read again, same hash
    set_mtime(bin, 1000000000);
    CHECK_EQ(adaptive::binary_sha256(bin), first);
    CHECK_EQ(adaptive::binary_hash_reads() - start, 2LL);

    // same size, other content, mtime put back: ctime still moves, so it is read again and the hash differs. The ctime
    // a file system keeps can be coarse (a few ms here, measured 2026-10-08): two changes inside one tick look alike,
    // so the test waits longer than that; a rebuild takes seconds.
    usleep(50000);
    write(bin, "other binary");
    set_mtime(bin, 1000000000);
    const std::string other = adaptive::binary_sha256(bin);
    CHECK(other != first);
    CHECK_EQ(other, adaptive::sha256_file(bin));

    // replaced by rename (another inode): read again
    const std::string fresh = std::string(dir) + "/fresh";
    write(fresh, "third binary content");
    CHECK_EQ(std::rename(fresh.c_str(), bin.c_str()), 0);
    CHECK_EQ(adaptive::binary_sha256(bin), adaptive::sha256_file(bin));

    // gone: an OSError, not the old hash
    std::remove(bin.c_str());
    CHECK_RAISES(adaptive::binary_sha256(bin), adaptive::Error);
    rmdir(dir);
    return check::report("test_binary_hash");
}
