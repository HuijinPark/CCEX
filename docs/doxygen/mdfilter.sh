#!/bin/sh
# Doxygen input filter for the Markdown pages (docs/user and docs/methods).
# Only POSIX sed and awk features: the CI action runs it with busybox.
#  - A link to a file that is not part of the site (docs/validation, docs/development)
#    becomes its text followed by the repository path in code, so the site has no dead
#    link while GitHub keeps the real one.
#  - A link ../../<path> into the repository becomes a link to that path in the public
#    repository on GitHub.
#  - A link written as [text](Page.md#anchor) or [text](../user/Page.md#anchor) works on
#    GitHub, but Doxygen cannot turn it into a page link; it can resolve the bare anchor,
#    which is global. So the path and page part are dropped before Doxygen reads the file.
#  - An image written as ![text](../images/<dir>/file.png) is found by GitHub through
#    that path, but Doxygen copies an image into the site only when it is named bare
#    and found in IMAGE_PATH. So the directory part is dropped too.
#  - Math is written for GitHub, $...$ inline and $$...$$ displayed (on one line or over
#    several). Doxygen reads \f$...\f$ and \f[...\f], so the delimiters are rewritten,
#    outside code blocks and code spans only. A $ that does not open math (no closing $
#    on the line, a space after it) is left alone, as GitHub does.
sed -E \
  -e 's#\[([^]]*)\]\(\.\./(validation|development)/([^)]*)\)#\1 (`docs/\2/\3`)#g' \
  -e 's#\]\(\.\./\.\./([^)]*)\)#](https://github.com/CQML-SKKU/CCEX/tree/main/\1)#g' \
  -e 's#\]\((\.\./(user|methods)/)?[A-Za-z0-9_-]+\.md\##](\##g' \
  -e 's#\]\(\.\./images/([A-Za-z0-9_]+/)*#](#g' \
  "$1" | awk '
BEGIN { fence = ""; disp = 0 }
function conv(s,    out, i, n, c, j, k, run, rl, body) {
    out = ""; i = 1; n = length(s)
    while (i <= n) {
        c = substr(s, i, 1)
        if (c == "`") {                          # code span: copy through the matching run
            j = i; while (j <= n && substr(s, j, 1) == "`") j++
            run = substr(s, i, j - i); rl = j - i
            k = index(substr(s, j), run)
            if (k > 0) { out = out substr(s, i, rl + k - 1 + rl); i = j + k - 1 + rl }
            else { out = out run; i = j }
            continue
        }
        if (disp) {                              # inside a displayed formula
            if (substr(s, i, 2) == "$$") { out = out "\\f]"; disp = 0; i += 2 }
            else { out = out c; i++ }
            continue
        }
        if (c == "\\" && substr(s, i + 1, 1) == "$") { out = out "\\$"; i += 2; continue }
        if (substr(s, i, 2) == "$$") {
            k = index(substr(s, i + 2), "$$")
            if (k > 0) { out = out "\\f[" substr(s, i + 2, k - 1) "\\f]"; i = i + 2 + k - 1 + 2 }
            else { out = out "\\f["; disp = 1; i += 2 }
            continue
        }
        if (c == "$" && i < n && substr(s, i + 1, 1) != " ") {
            j = i + 1; k = 0
            while (j <= n) {
                if (substr(s, j, 1) == "\\") { j += 2; continue }
                if (substr(s, j, 1) == "$") {
                    if (substr(s, j - 1, 1) != " " && substr(s, j + 1, 1) !~ /[0-9]/) k = j
                    break
                }
                j++
            }
            if (k > 0) { body = substr(s, i + 1, k - i - 1); out = out "\\f$" body "\\f$"; i = k + 1; continue }
        }
        out = out c; i++
    }
    return out
}
{
    line = $0
    if (fence == "" && (line ~ /^[ ]*```[^`]*$/ || line ~ /^[ ]*~~~/)) { fence = (line ~ /^[ ]*```/) ? "```" : "~~~"; print; next }
    if (fence != "") { if (index(line, fence) > 0 && line ~ /^[ ]*(```|~~~)[ ]*$/) fence = ""; print; next }
    print conv(line)
}'
