#!/bin/sh
# Doxygen input filter for the Markdown pages (docs/user and docs/methods).
#  - A link to a file that is not part of the site (docs/validation, docs/development,
#    or ../../ into the repository) becomes its text followed by the repository path in
#    code, so the site has no dead link while GitHub keeps the real one.
#  - A link written as [text](Page.md#anchor) or [text](../user/Page.md#anchor) works on
#    GitHub, but Doxygen cannot turn it into a page link; it can resolve the bare anchor,
#    which is global. So the path and page part are dropped before Doxygen reads the file.
#  - A page link [text](../user/Page.md) or [text](../methods/page.md) loses its
#    directory: Doxygen finds a page by its file name.
#  - An image written as ![text](../images/<dir>/file.png) is found by GitHub through
#    that path, but Doxygen copies an image into the site only when it is named bare
#    and found in IMAGE_PATH. So the directory part is dropped too.
sed -E \
  -e 's#\[([^]]*)\]\(\.\./(validation|development)/([^)]*)\)#\1 (`docs/\2/\3`)#g' \
  -e 's#\[([^]]*)\]\(\.\./\.\./([^)]*)\)#\1 (`\2`)#g' \
  -e 's#\]\((\.\./(user|methods)/)?[A-Za-z0-9_-]+\.md\##](\##g' \
  -e 's#\]\(\.\./images/([A-Za-z0-9_]+/)*#](#g' \
  "$1"
