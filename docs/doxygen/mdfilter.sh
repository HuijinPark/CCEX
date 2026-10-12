#!/bin/sh
# Doxygen input filter for the Markdown pages.
#  - A link written as [text](Page.md#anchor) works on GitHub, but Doxygen cannot turn
#    it into a page link; it can resolve the bare anchor, which is global. So the page
#    part is dropped before Doxygen reads the file.
#  - An image written as ![text](../images/<dir>/file.png) is found by GitHub through
#    that path, but Doxygen copies an image into the site only when it is named bare
#    and found in IMAGE_PATH. So the directory part is dropped too.
sed -E -e 's/\]\([A-Za-z0-9_]+\.md#/](#/g' -e 's#\]\(\.\./images/([A-Za-z0-9_]+/)*#](#g' "$1"
