// Writing what Python reads back unchanged, in value, type and sign (plan 0B-1, Task 6).
#include "check.h"
#include "base_json.h"
#include "mkdirs.h"
#include <cmath>

int main() {
    adaptive::Value root = adaptive::Value::object();
    root.set("work", adaptive::Value::integer(static_cast<long long>(1) << 62));
    root.set("neg", adaptive::Value::integer(-(static_cast<long long>(1) << 62)));
    root.set("small", adaptive::Value::integer(static_cast<long long>(0)));
    root.set("one", adaptive::Value::number(1.));          // must not come back as int 1
    root.set("two", adaptive::Value::number(2.));
    root.set("pos_zero", adaptive::Value::number(0.));     // must not come back as int 0
    root.set("neg_zero", adaptive::Value::number(-0.));    // must keep its sign
    root.set("tenth", adaptive::Value::number(.1));
    root.set("big", adaptive::Value::number(1e308));
    root.set("tiny", adaptive::Value::number(5e-324));
    root.set("inf", adaptive::Value::number(HUGE_VAL));
    root.set("ninf", adaptive::Value::number(-HUGE_VAL));
    root.set("nan", adaptive::Value::number(std::nan("")));
    adaptive::Value list = adaptive::Value::array();
    list.push(adaptive::Value::text("112\\3\n4"));
    list.push(adaptive::Value::boolean(true));
    list.push(adaptive::Value());
    root.set("mixed", list);
    // Strings built here, never through the reader: a reader and a writer that
    // make the same mistake cancel out when every string arrives via the reader.
    adaptive::Value texts = adaptive::Value::object();
    texts.set("nul", adaptive::Value::text(std::string("a\0b", 3)));
    texts.set("accented", adaptive::Value::text("\xc3\xa9"));          // U+00E9
    texts.set("emoji", adaptive::Value::text("\xf0\x9f\x98\x80"));   // U+1F600
    texts.set("lone", adaptive::Value::text("\xed\xa0\x80"));         // WTF-8 U+D800
    texts.set("literal", adaptive::Value::text("\\u0041"));            // six characters, not 'A'
    root.set("texts", texts);
    root.set("dup", adaptive::Value::integer(static_cast<long long>(1)));
    root.set("dup", adaptive::Value::integer(static_cast<long long>(2)));   // last wins, first position
    adaptive::save_json("out/written.json", root, 2);

    const std::string text = adaptive::write_json(root, 2);
    CHECK_EQ(text.back(), '\n');
    CHECK(text.find("\"one\": 1.0") != std::string::npos);
    CHECK(text.find("\"neg_zero\": -0.0") != std::string::npos);
    CHECK(text.find("\"inf\": \"inf\"") != std::string::npos);
    CHECK(text.find("\"literal\": \"\\\\u0041\"") != std::string::npos);   // the backslash is escaped
    CHECK(text.find("\"lone\": \"\\ud800\"") != std::string::npos);         // WTF-8 goes back out as an escape
    CHECK_EQ(root.keys().size(), std::size_t(16));
    CHECK_EQ(root.keys().at(0), std::string("work"));
    CHECK_EQ(adaptive::decimal(root.at("dup").as_int()), std::string("2"));
    // Nesting on the way out: building a deep Value in memory and writing it used
    // to overflow the stack at 10,000 levels (10th review).
    {
        adaptive::Value node = adaptive::Value::array();
        for (int i = 0; i < 499; ++i) {
            adaptive::Value outer = adaptive::Value::array();
            outer.push(node);
            node = outer;
        }
        CHECK(!adaptive::write_json(node, 0).empty());
        adaptive::Value deeper = adaptive::Value::array();
        deeper.push(node);
        deeper.push(adaptive::Value::integer(static_cast<long long>(1)));
        adaptive::Value deepest = adaptive::Value::array();
        deepest.push(deeper);
        CHECK_THROWS(adaptive::write_json(deepest, 0), adaptive::RecursionError,
                     "maximum JSON nesting depth (500) exceeded");
    }

    // The file error carries Python's own class, not one blanket FileNotFoundError.
    testing::mkdirs("out/errors");
    adaptive::save_json("out/errors/plain.json", adaptive::Value::object(), 2);
    CHECK_THROWS(adaptive::save_json("out/errors/plain.json/deeper.json", adaptive::Value::object(), 2),
                 adaptive::Error, "[Errno 20] Not a directory: 'out/errors/plain.json/deeper.json.tmp'");
    testing::mkdirs("out/errors/directory.json.tmp");
    CHECK_THROWS(adaptive::save_json("out/errors/directory.json", adaptive::Value::object(), 2),
                 adaptive::Error, "[Errno 21] Is a directory: 'out/errors/directory.json.tmp'");
    CHECK_THROWS(adaptive::save_json("out/errors/missing/x.json", adaptive::Value::object(), 2),
                 adaptive::Error, "[Errno 2] No such file or directory: 'out/errors/missing/x.json.tmp'");
    return check::report("test_json_write");
}
