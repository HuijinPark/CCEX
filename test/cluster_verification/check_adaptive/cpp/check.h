// Minimal assertion harness for the refinement port (plan 0B-1). Test only.
// No MPI, no MKL: these tests must run on a laptop and inside a job alike.
#pragma once
#include <iostream>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

namespace check {
inline std::vector<std::string>& failures() {
    static std::vector<std::string> found;
    return found;
}

inline std::string show(const std::string& value) { return "\""+value+"\""; }

template <class T> std::string show(const T& value) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out.precision(17);
    out << value;
    return out.str();
}

// Containers are compared often enough (key order, file lists) that printing
// them belongs to the harness, not to each test.
template <class T> std::string show(const std::vector<T>& items) {
    std::string text = "[";
    for (std::size_t i = 0; i < items.size(); ++i)
        text += (i ? ", " : "")+show(items[i]);
    return text+"]";
}

inline void fail(const char* file, int line, const std::string& text) {
    failures().push_back(std::string(file)+":"+std::to_string(line)+" "+text);
}

inline int report(const char* name) {
    for (const std::string& text : failures())
        std::cout << "FAIL " << text << "\n";
    std::cout << name << ": " << failures().size() << " failures\n";
    return failures().empty() ? 0 : 1;
}
}   // namespace check

#define CHECK(cond) \
    do { if (!(cond)) check::fail(__FILE__, __LINE__, "not true: " #cond); } while (0)

#define CHECK_EQ(left, right) \
    do { auto l_ = (left); auto r_ = (right); \
         if (!(l_ == r_)) check::fail(__FILE__, __LINE__, \
             std::string(#left " == " #right ": ")+check::show(l_)+" vs "+check::show(r_)); } while (0)

// The exception type AND its message are part of the decision (design 3.4).
// Where the message is long and carries a file diff, only the kind is checked.
#define CHECK_RAISES(expr, Type) \
    do { bool raised_ = false; \
         try { (void)(expr); } \
         catch (const Type&) { raised_ = true; } \
         catch (...) { check::fail(__FILE__, __LINE__, #expr ": wrong exception type"); raised_ = true; } \
         if (!raised_) check::fail(__FILE__, __LINE__, #expr ": nothing raised"); } while (0)

#define CHECK_THROWS(expr, Type, message) \
    do { bool raised_ = false; \
         try { (void)(expr); } \
         catch (const Type& error_) { raised_ = true; \
             if (std::string(error_.what()) != std::string(message)) \
                 check::fail(__FILE__, __LINE__, std::string(#expr ": message ")+error_.what()+ \
                             " vs "+(message)); } \
         catch (...) { check::fail(__FILE__, __LINE__, #expr ": wrong exception type"); raised_ = true; } \
         if (!raised_) check::fail(__FILE__, __LINE__, #expr ": nothing raised"); } while (0)
