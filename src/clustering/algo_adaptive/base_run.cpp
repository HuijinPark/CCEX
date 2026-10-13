// Running CCEX once and keeping what came back: the launcher, the clock, SHA-256,
// the input manifest that identifies a run, and the on-disk cache of raw probe curves.
#include "base_run.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <chrono>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <limits.h>
#include <map>
#include <set>
#include <sqlite3.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include <vector>
#include "base_numeric.h"
#include "base_json.h"

// ----------------------------------------------------------------------------
// SHA-256 of a text or a file
// ----------------------------------------------------------------------------
namespace adaptive {
namespace {

const std::uint32_t kRound[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
    0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
    0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
    0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
    0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
    0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};

std::uint32_t rotate(std::uint32_t value, int bits) {
    return (value >> bits) | (value << (32-bits));
}

class Digest {
 public:
    void add(const unsigned char* data, std::size_t size) {
        length_ += size;
        while (size > 0) {
            const std::size_t room = 64-filled_;
            const std::size_t take = size < room ? size : room;
            for (std::size_t i = 0; i < take; ++i)
                block_[filled_+i] = data[i];
            filled_ += take;
            data += take;
            size -= take;
            if (filled_ == 64) {
                compress();
                filled_ = 0;
            }
        }
    }

    std::string hex() {
        const std::uint64_t bits = length_*8;
        unsigned char padding[72] = {0x80};
        std::size_t pad = filled_ < 56 ? 56-filled_ : 120-filled_;
        add(padding, pad);
        unsigned char tail[8];
        for (int i = 0; i < 8; ++i)
            tail[i] = static_cast<unsigned char>(bits >> (56-8*i));
        length_ -= pad+8;                      // the padding is not message length
        add(tail, 8);
        static const char* kHex = "0123456789abcdef";
        std::string out;
        for (int i = 0; i < 8; ++i)
            for (int shift = 28; shift >= 0; shift -= 4)
                out.push_back(kHex[(state_[i] >> shift) & 0xf]);
        return out;
    }

 private:
    void compress() {
        std::uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = (static_cast<std::uint32_t>(block_[i*4]) << 24) |
                   (static_cast<std::uint32_t>(block_[i*4+1]) << 16) |
                   (static_cast<std::uint32_t>(block_[i*4+2]) << 8) |
                   static_cast<std::uint32_t>(block_[i*4+3]);
        for (int i = 16; i < 64; ++i) {
            const std::uint32_t s0 = rotate(w[i-15], 7)^rotate(w[i-15], 18)^(w[i-15] >> 3);
            const std::uint32_t s1 = rotate(w[i-2], 17)^rotate(w[i-2], 19)^(w[i-2] >> 10);
            w[i] = w[i-16]+s0+w[i-7]+s1;
        }
        std::uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        std::uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
        for (int i = 0; i < 64; ++i) {
            const std::uint32_t s1 = rotate(e, 6)^rotate(e, 11)^rotate(e, 25);
            const std::uint32_t choose = (e & f)^((~e) & g);
            const std::uint32_t t1 = h+s1+choose+kRound[i]+w[i];
            const std::uint32_t s0 = rotate(a, 2)^rotate(a, 13)^rotate(a, 22);
            const std::uint32_t majority = (a & b)^(a & c)^(b & c);
            const std::uint32_t t2 = s0+majority;
            h = g; g = f; f = e; e = d+t1;
            d = c; c = b; b = a; a = t1+t2;
        }
        state_[0] += a; state_[1] += b; state_[2] += c; state_[3] += d;
        state_[4] += e; state_[5] += f; state_[6] += g; state_[7] += h;
    }

    std::uint32_t state_[8] = {0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
                               0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
    unsigned char block_[64] = {};
    std::size_t filled_ = 0;
    std::uint64_t length_ = 0;
};
}   // namespace

std::string sha256(const std::string& data) {
    Digest digest;
    digest.add(reinterpret_cast<const unsigned char*>(data.data()), data.size());
    return digest.hex();
}

namespace {
struct FileStamp {
    long long size = -1, mtime_s = 0, mtime_ns = 0, ctime_s = 0, ctime_ns = 0;
    unsigned long long inode = 0, device = 0;
    bool operator==(const FileStamp& o) const {
        return size == o.size && mtime_s == o.mtime_s && mtime_ns == o.mtime_ns && ctime_s == o.ctime_s &&
               ctime_ns == o.ctime_ns && inode == o.inode && device == o.device;
    }
};

FileStamp stamp_of(const std::string& path) {
    struct stat info;
    errno = 0;
    if (::stat(path.c_str(), &info) != 0)
        throw os_error(errno, path);
    FileStamp out;
    out.size = static_cast<long long>(info.st_size);
#if defined(__APPLE__)
    out.mtime_s = info.st_mtimespec.tv_sec; out.mtime_ns = info.st_mtimespec.tv_nsec;
    out.ctime_s = info.st_ctimespec.tv_sec; out.ctime_ns = info.st_ctimespec.tv_nsec;
#else
    out.mtime_s = info.st_mtim.tv_sec; out.mtime_ns = info.st_mtim.tv_nsec;
    out.ctime_s = info.st_ctim.tv_sec; out.ctime_ns = info.st_ctim.tv_nsec;
#endif
    out.inode = static_cast<unsigned long long>(info.st_ino);
    out.device = static_cast<unsigned long long>(info.st_dev);
    return out;
}

struct BinaryHash { FileStamp stamp; std::string sha; };
std::map<std::string, BinaryHash>& binary_hashes() { static std::map<std::string, BinaryHash> held; return held; }
long long whole_reads = 0;
}   // namespace

std::string binary_sha256(const std::string& path) {
    const FileStamp now = stamp_of(path);
    std::map<std::string, BinaryHash>::const_iterator held = binary_hashes().find(path);
    if (held != binary_hashes().end() && held->second.stamp == now)
        return held->second.sha;
    const std::string sha = sha256_file(path);
    ++whole_reads;
    if (!(stamp_of(path) == now))         // never keep a hash of a file that moved under the read
        throw RuntimeError("the binary changed while it was being hashed: "+path);
    binary_hashes()[path] = BinaryHash{now, sha};
    return sha;
}

long long binary_hash_reads() { return whole_reads; }

std::string sha256_file(const std::string& path) {
    errno = 0;
    std::FILE* file = std::fopen(path.c_str(), "rb");
    if (file == nullptr)
        throw os_error(errno, path);
    Digest digest;
    std::vector<unsigned char> buffer(65536);
    for (;;) {
        errno = 0;
        const std::size_t read = std::fread(buffer.data(), 1, buffer.size(), file);
        digest.add(buffer.data(), read);
        if (read < buffer.size()) {
            const int failure = std::ferror(file) ? (errno != 0 ? errno : EIO) : 0;
            std::fclose(file);
            if (failure != 0)
                throw os_error(failure, path);
            return digest.hex();
        }
    }
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// the clock the budgets read
// ----------------------------------------------------------------------------
namespace adaptive {
namespace {
double seconds(std::chrono::steady_clock::time_point at) {
    return std::chrono::duration<double>(at.time_since_epoch()).count();
}
}   // namespace

// Both are monotonic, and both are only ever used as differences here.
double SystemClock::monotonic() { return seconds(std::chrono::steady_clock::now()); }
double SystemClock::perf_counter() { return seconds(std::chrono::steady_clock::now()); }
}   // namespace adaptive

// ----------------------------------------------------------------------------
// launching CCEX for one step
// ----------------------------------------------------------------------------
extern char** environ;

namespace adaptive {
namespace {
void set(std::vector<std::string>* env, const std::string& name, const std::string& value,
         bool only_if_absent) {
    const std::string prefix = name+"=";
    for (std::string& entry : *env)
        if (entry.compare(0, prefix.size(), prefix) == 0) {
            if (!only_if_absent)
                entry = prefix+value;
            return;
        }
    env->push_back(prefix+value);
}

// waitpid without blocking, so the deadline can be watched. A hundredth of a
// second is far below any real run and far above the cost of asking.
bool finished(pid_t child, int* status, bool blocking) {
    const pid_t done = waitpid(child, status, blocking ? 0 : WNOHANG);
    if (done == child)
        return true;
    if (done < 0 && errno != EINTR)
        throw os_error(errno, "waitpid");
    return false;
}

void pause_briefly() {
    struct timespec how_long = {0, 10*1000*1000};
    nanosleep(&how_long, nullptr);
}

double now_seconds() {
    struct timespec at;
    clock_gettime(CLOCK_MONOTONIC, &at);
    return static_cast<double>(at.tv_sec)+static_cast<double>(at.tv_nsec)/1e9;
}
}   // namespace

std::vector<std::string> native_environment() {
    std::vector<std::string> env;
    for (char** at = environ; at != nullptr && *at != nullptr; ++at)
        env.push_back(*at);
    set(&env, "OMP_NUM_THREADS", "1", false);
    set(&env, "OPENBLAS_NUM_THREADS", "1", false);
    set(&env, "MKL_NUM_THREADS", "1", false);
    set(&env, "I_MPI_FABRICS", "shm", true);        // setdefault, not an override
    return env;
}

namespace { LocalRunner local_runner; }

void set_local_runner(LocalRunner runner) { local_runner = std::move(runner); }

NativeResult run_native(const std::vector<std::string>& command, const std::string& cwd,
                        const std::string& log_path, double timeout) {
    // In-process: no fork, and therefore no timeout to enforce - the calculation
    // this process runs cannot be killed without killing the decision with it.
    if (local_runner)
        return local_runner(command, cwd, log_path);
    const std::vector<std::string> env = native_environment();
    std::vector<char*> argv, envp;
    for (const std::string& word : command)
        argv.push_back(const_cast<char*>(word.c_str()));
    argv.push_back(nullptr);
    for (const std::string& entry : env)
        envp.push_back(const_cast<char*>(entry.c_str()));
    envp.push_back(nullptr);
    // execvpe is GNU only (macOS has none: review 20261012-233058-452278 R2), so the PATH
    // search is done here as glibc's execvpe does it: a name with a slash is run as given,
    // otherwise each PATH entry in turn (an empty entry is the working directory; no PATH:
    // /bin:/usr/bin), going past an entry that does not hold it or may not run it. A file
    // the kernel will not execute (ENOEXEC) is handed to /bin/sh. Everything is built
    // before fork; the child only calls execve, so relative entries resolve after its
    // chdir, as before.
    std::vector<std::string> candidates;
    if (command.empty() || command[0].find('/') != std::string::npos) {
        candidates.push_back(command.empty() ? std::string() : command[0]);
    } else {
        const char* path = std::getenv("PATH");
        const std::string dirs = path != nullptr ? path : "/bin:/usr/bin";
        for (std::size_t at = 0; at <= dirs.size(); ) {
            const std::size_t colon = std::min(dirs.find(':', at), dirs.size());
            const std::string dir = dirs.substr(at, colon-at);
            candidates.push_back(dir.empty() ? command[0] : dir+"/"+command[0]);
            at = colon+1;
        }
    }
    std::vector<char*> shell_argv;
    shell_argv.push_back(const_cast<char*>("/bin/sh"));
    shell_argv.push_back(nullptr);                       // the file, set in the child
    for (std::size_t k = 1; k < argv.size(); ++k)
        shell_argv.push_back(argv[k]);                   // ends with argv's nullptr

    const pid_t child = fork();
    if (child < 0)
        throw os_error(errno, "fork");
    if (child == 0) {
        // A session of its own, so the timeout can reach mpiexec's children.
        setsid();
        if (!cwd.empty() && chdir(cwd.c_str()) != 0)
            _exit(127);
        const int log = open(log_path.c_str(), O_WRONLY|O_CREAT|O_TRUNC, 0644);
        if (log < 0)
            _exit(127);
        dup2(log, STDOUT_FILENO);
        dup2(log, STDERR_FILENO);          // stderr joins stdout, as Python asks
        if (log > STDERR_FILENO)
            close(log);
        for (const std::string& file : candidates) {
            execve(file.c_str(), argv.data(), envp.data());
            if (errno == ENOEXEC) {
                shell_argv[1] = const_cast<char*>(file.c_str());
                execve("/bin/sh", shell_argv.data(), envp.data());
                break;
            }
            if (errno != EACCES && errno != ENOENT && errno != ENOTDIR && errno != ESTALE &&
                errno != ENODEV && errno != ETIMEDOUT)
                break;
        }
        _exit(127);
    }

    NativeResult out;
    int status = 0;
    const double deadline = now_seconds()+timeout;
    bool done = false;
    while (!(done = finished(child, &status, false))) {
        if (now_seconds() >= deadline) {
            out.timed_out = true;
            kill(-child, SIGTERM);
            const double grace = now_seconds()+5.;
            while (!(done = finished(child, &status, false)) && now_seconds() < grace)
                pause_briefly();
            if (!done) {
                kill(-child, SIGKILL);
                while (!finished(child, &status, true)) {}
            }
            break;
        }
        pause_briefly();
    }
    // Python reads Popen.returncode, which is -N when a signal ended it.
    if (WIFEXITED(status))
        out.returncode = WEXITSTATUS(status);
    else if (WIFSIGNALED(status))
        out.returncode = -WTERMSIG(status);
    return out;
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// the input manifest: which files a run depends on, and their hashes
// ----------------------------------------------------------------------------
namespace adaptive {
namespace {
// the input keys that name files, as paths into the input's blocks (A47); the order is the manifest's
const std::vector<std::vector<const char*>> kInputFiles = {
    {"gyrofile"}, {"Qubit", "file"}, {"bathfile"}, {"bath_state_file"}, {"Defect", "state_file"},
    {"Defect", "axis_file"}, {"hamiltonian", "common", "hyperfine", "file"},
    {"hamiltonian", "common", "quadrupole", "file"},
    {"hamiltonian", "common", "quadrupole", "file_without_defect"}};

std::string parent_of(const std::string& path) {
    const std::size_t slash = path.rfind('/');
    if (slash == std::string::npos || slash == 0)
        return "/";
    return path.substr(0, slash);
}

std::vector<std::string> split(const std::string& text) {
    std::vector<std::string> parts;
    for (std::size_t at = 0; at <= text.size(); ) {
        const std::size_t slash = text.find('/', at);
        const std::size_t end = slash == std::string::npos ? text.size() : slash;
        if (end > at)
            parts.push_back(text.substr(at, end-at));
        if (slash == std::string::npos)
            break;
        at = slash+1;
    }
    return parts;
}

// One component at a time, following a symlink the moment it is met. That order
// matters: "link/.." is the parent of what the link points at, not of the link.
//
// A link that leads back to itself is left exactly as written rather than raising:
// resolve(strict=False) gives .../loop for a link called loop that points at loop,
// and .../a for a pair a -> b -> a (measured against CPython 3.13.10, 2026-09-16).
void walk(const std::string& text, std::string* out, std::set<std::string>* resolving) {
    for (const std::string& part : split(text)) {
        if (part == ".")
            continue;
        if (part == "..") {
            *out = parent_of(*out);
            continue;
        }
        const std::string here = (*out == "/" ? "" : *out)+"/"+part;
        struct stat info;
        if (lstat(here.c_str(), &info) == 0 && S_ISLNK(info.st_mode) &&
            resolving->count(here) == 0) {
            std::vector<char> target(static_cast<std::size_t>(info.st_size)+1, '\0');
            const ssize_t got = readlink(here.c_str(), target.data(), target.size()-1);
            if (got < 0)
                throw os_error(errno, here);
            const std::string link(target.data(), static_cast<std::size_t>(got));
            resolving->insert(here);
            // A relative target resolves against the directory holding the link,
            // which is exactly where *out already stands.
            if (!link.empty() && link[0] == '/')
                *out = "/";
            walk(link, out, resolving);
            resolving->erase(here);
            continue;
        }
        *out = here;
    }
}
}   // namespace

std::string resolved(const std::string& cwd, const std::string& name) {
    std::string out = "/";
    std::set<std::string> resolving;
    if (!name.empty() && name[0] == '/') {
        walk(name, &out, &resolving);
        return out;
    }
    // A relative base is relative to where the process is, the way Path(cwd) is
    // before resolve() is called on it. --cwd defaults to an absolute path but
    // nothing stops a relative one being given.
    if (cwd.empty() || cwd[0] != '/') {
        std::vector<char> here(4096, '\0');
        if (getcwd(here.data(), here.size()) == nullptr)
            throw os_error(errno, ".");
        walk(here.data(), &out, &resolving);
    }
    walk(cwd, &out, &resolving);
    walk(name, &out, &resolving);
    return out;
}

Value input_manifest(const Value& cfg, const std::string& cwd, const std::vector<int>& states) {
    Value files = Value::object();
    for (const std::vector<const char*>& key : kInputFiles) {
        std::vector<std::string> paths;
        const Value* at = &cfg;
        for (const char* step : key)
            at = at && at->kind() == Value::Kind::Object ? at->find(step) : nullptr;
        if (at && at->kind() != Value::Kind::Null) {
            const Value& given = *at;
            if (given.kind() == Value::Kind::String) {
                paths.push_back(given.as_string());
            } else {
                for (std::size_t i = 0; i < given.size(); ++i)
                    paths.push_back(given.at(i).as_string());
            }
        }
        const std::string name(key.back());
        if (name == "bath_state_file" || name == "state_file") {
            std::vector<std::string> spread;
            for (const std::string& path : paths)
                for (int state : states)
                    if (state != 0)                 // Python's "for s in states if s"
                        spread.push_back(path+std::to_string(state));
            paths = spread;
        }
        for (const std::string& path : paths) {
            const std::string here = resolved(cwd, path);
            struct stat info;
            const bool is_file = stat(here.c_str(), &info) == 0 && S_ISREG(info.st_mode);
            files.set(here, is_file ? Value::text(sha256_file(here)) : Value());
        }
    }
    return files;
}

Value identity(const Value& cfg, const std::string& cwd, const std::vector<int>& states,
               const std::string& binary_sha256, long long ranks) {
    // These choose where output goes or how bundles are drawn, not the dynamics of
    // an explicit site list, so two runs that differ only here share probes: in the
    // bundle block, all but its expansion limits (max_clusters, max_work).
    static const std::set<std::string> ignored{"outfile", "adaptive", "profile", "savemode"};
    static const std::set<std::string> ignoredBundle{"max_dim", "merge_cutoff", "connections", "connection_budget"};
    Value config = Value::object();
    for (const std::string& key : cfg.keys()) {
        if (ignored.count(key) != 0)
            continue;
        if (key == "bundle" && cfg.at(key).kind() == Value::Kind::Object) {
            Value kept = Value::object();
            for (const std::string& inner : cfg.at(key).keys())
                if (ignoredBundle.count(inner) == 0)
                    kept.set(inner, cfg.at(key).at(inner));
            if (!kept.keys().empty())
                config.set(key, kept);
            continue;
        }
        config.set(key, cfg.at(key));
    }
    Value out = Value::object();
    out.set("schema", Value::integer(static_cast<long long>(1)));
    out.set("config", config);
    out.set("cwd", Value::text(resolved(cwd, ".")));
    Value listed = Value::array();
    for (int state : states)
        listed.push(Value::integer(static_cast<long long>(state)));
    out.set("states", listed);
    out.set("ranks", Value::integer(static_cast<long long>(ranks)));
    out.set("binary_sha256", Value::text(binary_sha256));
    out.set("files", input_manifest(cfg, cwd, states));
    return out;
}

// Python compares the two manifests as dicts. They hold nothing but strings and
// nulls, so the sorted, compact text says the same thing and is easy to read in a
// message.
std::string manifest_text(const Value& files) {
    return json_dumps(files, true, true, false);
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// the on-disk cache of raw probe curves
// ----------------------------------------------------------------------------
namespace adaptive {
namespace {

void run(sqlite3* db, const char* statement) {
    char* message = nullptr;
    if (sqlite3_exec(db, statement, nullptr, nullptr, &message) != SQLITE_OK) {
        const std::string text = message ? message : "sqlite error";
        sqlite3_free(message);
        throw Error("OperationalError", text);
    }
}

// Python's sqlite3 raises on every result that is neither a row nor the end.
// Dropping the result codes read a locked, full or damaged cache as "nothing
// cached", and let put() commit a response it had only half written.
sqlite3_stmt* prepared(sqlite3* db, const char* statement) {
    sqlite3_stmt* out = nullptr;
    if (sqlite3_prepare_v2(db, statement, -1, &out, nullptr) != SQLITE_OK)
        throw Error("OperationalError", sqlite3_errmsg(db));
    return out;
}

// A row, or the end - and the statement is finalized before anything is thrown,
// so a caller's ROLLBACK is not held up by it.
bool stepped(sqlite3* db, sqlite3_stmt* statement) {
    const int code = sqlite3_step(statement);
    if (code == SQLITE_ROW || code == SQLITE_DONE)
        return code == SQLITE_ROW;
    const std::string message = sqlite3_errmsg(db);
    sqlite3_finalize(statement);
    throw Error("OperationalError", message);
}

void make_parents(const std::string& path) {
    // Python's ProbeCache makes its parent directories; save_json does not. Both
    // are reproduced as they are.
    for (std::size_t i = 1; i < path.size(); ++i)
        if (path[i] == '/')
            mkdir(path.substr(0, i).c_str(), 0777);
}

// Python resolves the cache path before recording it in the report.
std::string resolved_path(const std::string& path) {
    char* real = realpath(path.c_str(), nullptr);
    if (real == nullptr)
        return path;
    const std::string out(real);
    std::free(real);
    return out;
}

std::string sites_text(const Sites& sites) {
    // json.dumps(tuple) - the DEFAULT separators, unlike the context payload.
    Value list = Value::array();
    for (int site : sites)
        list.push(Value::integer(static_cast<long long>(site)));
    return json_dumps(list, false, false, true);
}

std::string curves_text(const Curves& curves) {
    Value outer = Value::array();
    for (const Curve& curve : curves) {
        Value points = Value::array();
        for (const std::pair<double, double>& point : curve) {
            Value pair = Value::array();
            pair.push(Value::text(float_text(point.first)));    // repr keeps every double exactly,
            pair.push(Value::text(float_text(point.second)));   // NaN and infinity included
            points.push(pair);
        }
        outer.push(points);
    }
    return json_dumps(outer, false, false, true);
}
}   // namespace

ProbeCache::ProbeCache(const std::string& path, const Value& context) {
    path_ = resolved_path(path);
    const std::string payload = json_dumps(context, true, true, false);
    key_ = sha256(payload);
    states_ = context.at("states").size();
    steps_ = static_cast<long long>(context.at("config").at("nstep").as_int());
    make_parents(path);
    if (sqlite3_open(path.c_str(), &db_) != SQLITE_OK)
        throw Error("OperationalError", std::string("unable to open ")+path);
    sqlite3_busy_timeout(db_, 60000);                  // Python's timeout=60
    // A constructor that throws has no destructor run for it: closing here is what
    // releases the file, and with it any transaction left open.
    try {
        run(db_, "BEGIN");
        run(db_, "CREATE TABLE IF NOT EXISTS contexts (id TEXT PRIMARY KEY, payload TEXT NOT NULL)");
        run(db_, "CREATE TABLE IF NOT EXISTS probes (context TEXT, sites TEXT, payload TEXT NOT NULL,"
                 " sha256 TEXT NOT NULL, PRIMARY KEY (context, sites))");
        sqlite3_stmt* insert = prepared(db_, "INSERT OR IGNORE INTO contexts VALUES (?,?)");
        sqlite3_bind_text(insert, 1, key_.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(insert, 2, payload.c_str(), -1, SQLITE_TRANSIENT);
        stepped(db_, insert);
        sqlite3_finalize(insert);
        run(db_, "COMMIT");

        sqlite3_stmt* select = prepared(db_, "SELECT payload FROM contexts WHERE id=?");
        sqlite3_bind_text(select, 1, key_.c_str(), -1, SQLITE_TRANSIENT);
        std::string stored;
        if (stepped(db_, select))
            stored = reinterpret_cast<const char*>(sqlite3_column_text(select, 0));
        sqlite3_finalize(select);
        if (stored != payload)
            throw RuntimeError("raw cache context mismatch");
    } catch (...) {
        close();
        throw;
    }
}

ProbeCache::~ProbeCache() {
    close();
}

void ProbeCache::close() {
    if (db_ != nullptr) {
        sqlite3_close(db_);
        db_ = nullptr;
    }
}

std::map<Sites, Curves> ProbeCache::get(const std::vector<Sites>& wanted) {
    std::map<Sites, Curves> found;
    sqlite3_stmt* select = prepared(db_, "SELECT payload,sha256 FROM probes WHERE context=? AND sites=?");
    for (const Sites& sites : wanted) {
        const std::string text = sites_text(sites);
        sqlite3_reset(select);
        sqlite3_bind_text(select, 1, key_.c_str(), -1, SQLITE_TRANSIENT);
        sqlite3_bind_text(select, 2, text.c_str(), -1, SQLITE_TRANSIENT);
        if (!stepped(db_, select))
            continue;
        const std::string payload = reinterpret_cast<const char*>(sqlite3_column_text(select, 0));
        const std::string digest = reinterpret_cast<const char*>(sqlite3_column_text(select, 1));
        if (sha256(payload) != digest) {
            sqlite3_finalize(select);
            throw RuntimeError("raw cache checksum mismatch");
        }
        const Value data = read_json(payload);
        bool shaped = data.size() == states_;
        for (std::size_t i = 0; shaped && i < data.size(); ++i)
            shaped = static_cast<long long>(data.at(i).size()) == steps_;
        if (!shaped) {
            sqlite3_finalize(select);
            throw RuntimeError("raw cache curve shape mismatch");
        }
        Curves curves;
        for (std::size_t i = 0; i < data.size(); ++i) {
            Curve curve;
            for (std::size_t k = 0; k < data.at(i).size(); ++k) {
                const Value& point = data.at(i).at(k);
                curve.emplace_back(parse_float(point.at(std::size_t(0)).as_string()),
                                   parse_float(point.at(std::size_t(1)).as_string()));
            }
            curves.push_back(curve);
        }
        found[sites] = curves;
    }
    sqlite3_finalize(select);
    return found;
}

void ProbeCache::put(const std::map<Sites, Curves>& values) {
    // One transaction for the whole response: an interrupted call leaves nothing.
    run(db_, "BEGIN");
    try {
        for (const std::pair<const Sites, Curves>& entry : values) {
            const std::string text = sites_text(entry.first);
            const std::string payload = curves_text(entry.second);
            const std::string digest = sha256(payload);
            sqlite3_stmt* select = prepared(db_, "SELECT sha256 FROM probes WHERE context=? AND sites=?");
            sqlite3_bind_text(select, 1, key_.c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_text(select, 2, text.c_str(), -1, SQLITE_TRANSIENT);
            std::string existing;
            if (stepped(db_, select))
                existing = reinterpret_cast<const char*>(sqlite3_column_text(select, 0));
            sqlite3_finalize(select);
            if (!existing.empty() && existing != digest)
                throw RuntimeError("raw cache has conflicting results for identical conditions");
            sqlite3_stmt* insert = prepared(db_, "INSERT OR IGNORE INTO probes VALUES (?,?,?,?)");
            sqlite3_bind_text(insert, 1, key_.c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_text(insert, 2, text.c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_text(insert, 3, payload.c_str(), -1, SQLITE_TRANSIENT);
            sqlite3_bind_text(insert, 4, digest.c_str(), -1, SQLITE_TRANSIENT);
            stepped(db_, insert);
            sqlite3_finalize(insert);
        }
    } catch (...) {
        run(db_, "ROLLBACK");
        throw;
    }
    run(db_, "COMMIT");
}

}   // namespace adaptive
