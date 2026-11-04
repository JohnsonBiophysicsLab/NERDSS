#!/usr/bin/env bash
# Rebuild the online NERDSS user guide from its Markdown source.
#
#   1. refreshes the table of contents inside nerdss_user_guide.md
#      (between the <!-- toc --> and <!-- tocstop --> markers);
#   2. renders nerdss_user_guide.html with template.html + filters.lua.
#
# Needs pandoc (tested with 3.1) and python3. Run from anywhere:
#   ./user_guide/build.sh
set -euo pipefail
cd "$(dirname "$0")"

SRC=nerdss_user_guide.md
OUT=nerdss_user_guide.html
MATHJAX=https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js

strip_toc() { sed '/^<!-- toc -->$/,/^<!-- tocstop -->$/d' "$SRC"; }

# 1. Table of contents (levels 2-3) written back into the Markdown file.
toc=$(strip_toc | pandoc -f gfm -t gfm --lua-filter=filters.lua \
        --toc --toc-depth=3 --template=toc.tpl --wrap=none)
TOC="$toc" python3 - "$SRC" <<'PY'
import os, re, sys
path = sys.argv[1]
text = open(path, encoding='utf-8').read()
block = '<!-- toc -->\n## Contents\n\n' + os.environ['TOC'].strip() + '\n<!-- tocstop -->'
new, n = re.subn(r'<!-- toc -->\n.*?<!-- tocstop -->', lambda m: block, text, flags=re.S)
if n != 1:
    sys.exit('build.sh: expected exactly one <!-- toc --> ... <!-- tocstop --> block in ' + path)
if new != text:
    open(path, 'w', encoding='utf-8').write(new)
PY

# 2. Standalone HTML page in the website's look.
strip_toc | pandoc -f gfm -t html5 \
  --template=template.html \
  --lua-filter=filters.lua \
  --toc --toc-depth=3 \
  --mathjax="$MATHJAX" -V mathjax-url="$MATHJAX" \
  --no-highlight \
  --wrap=none \
  -o "$OUT"

echo "Built $(pwd)/$OUT"
