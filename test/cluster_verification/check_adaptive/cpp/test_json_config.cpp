// The one place we deliberately differ from Python: a config integer past int64 is refused.
#include "check.h"
#include "base_json.h"
#include <cmath>
#include <string>

int main(int argc, char* argv[]) {
    if (argc > 2 && std::string(argv[1]) == "--file") {
        adaptive::Value config = adaptive::read_config(adaptive::read_file(argv[2]));
        std::cout << argv[2] << ": " << config.keys().size() << " keys\n";
        return check::report("test_json_config");
    }
    if (argc > 2 && std::string(argv[1]) == "--verdicts") {
        // One line per file: accepted, or the exception it raised. Python's json has
        // a parse_int hook that fires on every integer token, so the same corpus can
        // be judged there and the two verdicts compared.
        for (const std::string& path : adaptive::list_json(argv[2])) {
            std::string verdict = "ok";
            try {
                adaptive::read_config(adaptive::read_file(std::string(argv[2])+"/"+path));
            } catch (const adaptive::Error& error) {
                verdict = error.type()+": "+error.what();
            }
            for (char& c : verdict)
                if (c == '\n' || c == '\t' || c == '\r') c = ' ';   // one verdict, one line
            std::cout << path << "\t" << verdict << "\n";
        }
        return check::report("test_json_config");
    }
    adaptive::Value ok = adaptive::read_config("{\"seed\": 19, \"order\": 2, \"deltat\": 0.01}");
    CHECK(ok.at("seed").kind() == adaptive::Value::Kind::Int);
    CHECK(ok.at("deltat").kind() == adaptive::Value::Kind::Double);

    adaptive::Value edge = adaptive::read_config("{\"a\": 9223372036854775807, \"b\": -9223372036854775808}");
    CHECK_EQ(adaptive::decimal(edge.at("a").as_int()), std::string("9223372036854775807"));
    CHECK_EQ(adaptive::decimal(edge.at("b").as_int()), std::string("-9223372036854775808"));

    CHECK_THROWS(adaptive::read_config("{\"a\": 9223372036854775808}"), adaptive::ValueError,
                 "integer 9223372036854775808 at a does not fit in 64 bits");
    CHECK_THROWS(adaptive::read_config("{\"a\": [{\"b\": -9223372036854775809}]}"), adaptive::ValueError,
                 "integer -9223372036854775809 at a[0].b does not fit in 64 bits");
    // An overwritten token is still a token: checking only the final value would let this pass.
    CHECK_THROWS(adaptive::read_config("{\"a\": 9223372036854775808, \"a\": 1}"), adaptive::ValueError,
                 "integer 9223372036854775808 at a does not fit in 64 bits");
    // A float of any size is fine: Python reads 1e400 as inf and so do we.
    CHECK_EQ(adaptive::read_config("{\"a\": 1e400}").at("a").as_double(), HUGE_VAL);
    // The path shown in this message escapes control characters, lone surrogates and
    // the backslash itself. A declared compatibility exception, and it goes no
    // further: the key kept in the value is untouched and the exception type is the
    // same (11th review). Without it a NUL truncated the message at what(),
    // and a WTF-8 key made the log line invalid UTF-8.
    CHECK_THROWS(adaptive::read_config("{\"a\\u0000b\": 9223372036854775808}"), adaptive::ValueError,
                 "integer 9223372036854775808 at a\\x00b does not fit in 64 bits");
    CHECK_THROWS(adaptive::read_config("{\"a\\ud800\": 9223372036854775808}"), adaptive::ValueError,
                 "integer 9223372036854775808 at a\\ud800 does not fit in 64 bits");
    CHECK_THROWS(adaptive::read_config("{\"a\\\\x00\": 9223372036854775808}"), adaptive::ValueError,
                 "integer 9223372036854775808 at a\\\\x00 does not fit in 64 bits");

    // Asking an array for a key is a TypeError, as in Python, not a KeyError:
    // it decides which catch the caller lands in.
    adaptive::Value list = adaptive::read_json("[1, 2]");
    CHECK_THROWS(list.at("a"), adaptive::TypeError, "expected dict, found list");
    CHECK_THROWS(list.has("a"), adaptive::TypeError, "expected dict, found list");
    adaptive::Value object = adaptive::read_json("{\"a\": 1}");
    CHECK_THROWS(object.at(std::size_t(0)), adaptive::TypeError,
                 "dict indices must be integers or slices");
    CHECK_THROWS(object.at("b"), adaptive::KeyError, "'b'");
    CHECK_THROWS(list.at(std::size_t(5)), adaptive::IndexError, "list index out of range");

    // Reading a directory came back as an empty string, because fread's failure was
    // never checked: a silent success is worse than an error (10th review).
    CHECK_THROWS(adaptive::read_file("."), adaptive::Error, "[Errno 21] Is a directory: '.'");
    CHECK_THROWS(adaptive::read_file("no-such-file.json"), adaptive::Error,
                 "[Errno 2] No such file or directory: 'no-such-file.json'");
    CHECK_THROWS(adaptive::read_file("check.h/inside.json"), adaptive::Error,
                 "[Errno 20] Not a directory: 'check.h/inside.json'");
    return check::report("test_json_config");
}
