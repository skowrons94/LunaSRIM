#!/bin/sh
# Download the JavaScript libraries of the web application into
# lunasrim/web/static/vendor so that it works without internet access
# (otherwise the browser loads them from a CDN):
#   plotly.js 2.35.2  (plots)
#   marked 12.0.2     (Markdown of the Physics page)
#   KaTeX 0.16.11     (formulas of the Physics page: js, css, auto-render, woff2 fonts)
# The KaTeX fonts are put next to katex.min.css (vendor/, no fonts/
# sub-directory, so that the package data pattern static/vendor/* covers
# them); the url(fonts/...) references of the css are rewritten accordingly.
# usage: tools/fetch_web_vendor.sh [--force]
set -e
dest="$(cd "$(dirname "$0")/.." && pwd)/lunasrim/web/static/vendor"
mkdir -p "$dest"
force="$1"
CDN=https://cdn.jsdelivr.net/npm

download() {   # url file
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL "$1" -o "$2"
    else
        python3 -c "import urllib.request,sys; urllib.request.urlretrieve(sys.argv[1], sys.argv[2])" "$1" "$2"
    fi
}

# npm_unpack PACKAGE@VERSION -> prints the directory of the unpacked package
npm_unpack() {
    command -v npm >/dev/null 2>&1 || return 1
    t=$(mktemp -d)
    (cd "$t" && npm pack "$1" >/dev/null 2>&1 && tar xzf ./*.tgz) || { rm -rf "$t"; return 1; }
    echo "$t/package"
}

have() { [ -s "$dest/$1" ] && [ "$force" != "--force" ]; }

# ------------------------------------------------------------------ plotly
if have plotly.min.js; then
    echo "plotly.js already present in $dest"
elif p=$(npm_unpack plotly.js-dist-min@2.35.2); then
    cp "$p/plotly.min.js" "$dest/plotly.min.js"; rm -rf "$(dirname "$p")"
    echo "plotly.js saved to $dest/plotly.min.js"
else
    download $CDN/plotly.js-dist-min@2.35.2/plotly.min.js "$dest/plotly.min.js"
    echo "plotly.js saved to $dest/plotly.min.js"
fi

# ------------------------------------------------------------------ marked
if have marked.min.js; then
    echo "marked already present in $dest"
elif p=$(npm_unpack marked@12.0.2); then
    cp "$p/marked.min.js" "$dest/marked.min.js"; rm -rf "$(dirname "$p")"
    echo "marked saved to $dest/marked.min.js"
else
    download $CDN/marked@12.0.2/marked.min.js "$dest/marked.min.js"
    echo "marked saved to $dest/marked.min.js"
fi

# ------------------------------------------------------------------- KaTeX
if have katex.min.js && have katex.min.css && have auto-render.min.js; then
    echo "KaTeX already present in $dest"
else
    if p=$(npm_unpack katex@0.16.11); then
        cp "$p/dist/katex.min.js" "$p/dist/katex.min.css" "$p/dist/contrib/auto-render.min.js" "$dest/"
        cp "$p/dist/fonts/"*.woff2 "$dest/"
        rm -rf "$(dirname "$p")"
    else
        k=$CDN/katex@0.16.11/dist
        download $k/katex.min.js "$dest/katex.min.js"
        download $k/katex.min.css "$dest/katex.min.css"
        download $k/contrib/auto-render.min.js "$dest/auto-render.min.js"
        for f in $(grep -o 'fonts/KaTeX_[A-Za-z0-9_-]*\.woff2' "$dest/katex.min.css" | sort -u); do
            download "$k/$f" "$dest/${f#fonts/}"
        done
    fi
    # fonts next to the css; woff2 only (every browser that runs the app reads it)
    python3 - "$dest/katex.min.css" <<'EOF'
import re, sys
p = sys.argv[1]
css = open(p, encoding='utf-8').read()
css = re.sub(r'src:url\(fonts/([^)]+\.woff2)\) format\("woff2"\)[^;}]*', r'src:url(\1) format("woff2")', css)
css = css.replace('url(fonts/', 'url(')
open(p, 'w', encoding='utf-8').write(css)
EOF
    echo "KaTeX saved to $dest (katex.min.js, katex.min.css, auto-render.min.js, fonts)"
fi
