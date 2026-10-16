/**
 * @file options.cpp
 * @brief Adaptive refinement: reading and checking the options.
 */
#include "options.h"

#include "options_table.h"
#include <cctype>
#include <vector>
#include "base_json.h"

namespace adaptive {
namespace {

const options::Option* find(const std::string& name) {
    for (std::size_t i = 0; i < options::count(); ++i)
        if (name == options::kOptions[i].name)
            return &options::kOptions[i];
    return nullptr;
}

// str(pathlib.Path(text)). argparse hands the option's type the raw word, and for
// these options that type is Path, which normalises: repeated slashes collapse, "."
// components go, a trailing slash goes, and an empty path becomes ".". What it does
// NOT do is resolve ".." or "~". POSIX keeps exactly two leading slashes, and drops
// three or more to one (measured against CPython 3.13.10, 2026-09-16).
std::string path_text(const std::string& text) {
    std::size_t leading = 0;
    while (leading < text.size() && text[leading] == '/')
        ++leading;
    const std::string root = leading == 0 ? "" : leading == 2 ? "//" : "/";
    std::vector<std::string> parts;
    for (std::size_t at = leading; at <= text.size(); ) {
        const std::size_t slash = text.find('/', at);
        const std::size_t end = slash == std::string::npos ? text.size() : slash;
        const std::string part = text.substr(at, end-at);
        if (!part.empty() && part != ".")
            parts.push_back(part);
        if (slash == std::string::npos)
            break;
        at = slash+1;
    }
    std::string out = root;
    for (std::size_t i = 0; i < parts.size(); ++i)
        out += (i == 0 ? "" : "/")+parts[i];
    if (out.empty())
        return ".";
    return out;
}

Value convert(const options::Option& option, const std::string& text, const std::string& cwd) {
    switch (option.kind) {
        case options::Kind::Int:
            try {
                return Value::integer(parse_int(text));
            } catch (const ValueError&) {
                throw ValueError(std::string("argument --")+option.name+
                                   ": invalid int value: '"+text+"'");
            }
        case options::Kind::Float:
            try {
                return Value::number(parse_float(text));
            } catch (const ValueError&) {
                throw ValueError(std::string("argument --")+option.name+
                                   ": invalid float value: '"+text+"'");
            }
        case options::Kind::Path:
            return Value::text(text == "<cwd>" ? cwd : path_text(text));
        default:
            return Value::text(text);
    }
}
}   // namespace

namespace {

// argparse decides whether a word is a value or an option before it converts
// anything, and a word starting with '-' is an option unless it is a plain negative
// number. So "--atol -1" is a tolerance of -1 but "--atol -1e-8" is an error: the
// exponent takes it out of argparse's _negative_number_matcher, ^-\d+$|^-\d*\.\d+$
// (measured against CPython 3.13.10, 2026-09-16).
bool looks_like_option(const std::string& word) {
    if (word.size() < 2 || word[0] != '-')
        return false;
    std::size_t at = 1;
    bool digit = false, dot = false;
    for (; at < word.size(); ++at) {
        if (std::isdigit(static_cast<unsigned char>(word[at]))) {
            digit = true;
        } else if (word[at] == '.' && !dot) {
            dot = true;
        } else {
            break;
        }
    }
    const bool negative_number = at == word.size() && digit && (!dot || word.size() > 2);
    return !negative_number;
}

bool allowed(const std::string& choices, const std::string& value) {
    std::size_t at = 0;
    while (at <= choices.size()) {
        const std::size_t bar = choices.find('|', at);
        const std::size_t end = bar == std::string::npos ? choices.size() : bar;
        if (choices.substr(at, end-at) == value)
            return true;
        if (bar == std::string::npos)
            return false;
        at = bar+1;
    }
    return false;
}

std::string listed(const std::string& choices) {
    std::string out = choices;
    for (std::size_t at = out.find('|'); at != std::string::npos; at = out.find('|', at+2))
        out.replace(at, 1, ", ");
    return out;
}

}   // namespace

Value parse_options(const std::vector<std::string>& argv, const std::string& working_directory) {
    Value parsed = Value::object();
    std::map<std::string, bool> seen;
    // argparse lays every default down first, in the order the options were
    // declared, and parsing then overwrites in place. That is why vars(args) is in
    // declaration order and not command-line order - and the report copies
    // vars(args) straight into summary.json, so the order is part of the file.
    for (std::size_t i = 0; i < options::count(); ++i) {
        const options::Option& option = options::kOptions[i];
        if (option.kind == options::Kind::Flag)
            parsed.set(option.dest, Value::boolean(option.default_text[0] != '\0'));
        else if (option.no_default || option.required)
            parsed.set(option.dest, Value());
        else if (option.int_default && option.kind == options::Kind::Float)
            parsed.set(option.dest, Value::integer(parse_int(option.default_text)));
        else
            parsed.set(option.dest, convert(option, option.default_text, working_directory));
    }
    for (std::size_t i = 0; i < argv.size(); ++i) {
        const std::string& word = argv[i];
        if (word.rfind("--", 0) != 0)
            throw ValueError("unrecognized arguments: "+word);
        std::string name = word.substr(2), value;
        bool inline_value = false;
        const std::size_t equals = name.find('=');
        if (equals != std::string::npos) {
            value = name.substr(equals+1);
            name = name.substr(0, equals);
            inline_value = true;
        }
        const options::Option* option = find(name);
        if (option == nullptr)
            throw ValueError("unrecognized arguments: "+word);
        if (option->kind == options::Kind::Flag) {
            if (inline_value)
                throw ValueError(std::string("argument --")+name+": ignored explicit argument");
            parsed.set(option->dest, Value::boolean(true));
        } else {
            if (!inline_value) {
                if (i+1 >= argv.size() || looks_like_option(argv[i+1]))
                    throw ValueError(std::string("argument --")+name+": expected one argument");
                value = argv[++i];
            }
            if (option->choices[0] != '\0' && !allowed(option->choices, value))
                throw ValueError(std::string("argument --")+name+": invalid choice: '"+
                                   value+"' (choose from "+listed(option->choices)+")");
            parsed.set(option->dest, convert(*option, value, working_directory));
        }
        seen[option->dest] = true;
    }
    std::string missing;
    for (std::size_t i = 0; i < options::count(); ++i) {
        const options::Option& option = options::kOptions[i];
        if (seen.count(option.dest) != 0)
            continue;
        // The defaults are already in place; what is left is to notice a required
        // option nobody gave.
        if (option.required)
            missing += (missing.empty() ? "" : ", ")+std::string("--")+option.name;
    }
    if (!missing.empty())
        throw ValueError("the following arguments are required: "+missing);
    return parsed;
}
}   // namespace adaptive
