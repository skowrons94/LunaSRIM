#!/bin/sh
# Download plotly.js into lunasrim/web/static/vendor so the web application
# works without internet access (otherwise the browser loads it from a CDN).
set -e
dest="$(dirname "$0")/../lunasrim/web/static/vendor"
mkdir -p "$dest"
if [ -s "$dest/plotly.min.js" ] && [ "$1" != "--force" ]; then
    echo "plotly.js already present in $dest"
    exit 0
fi
if command -v npm >/dev/null 2>&1; then
    tmp=$(mktemp -d)
    (cd "$tmp" && npm pack plotly.js-dist-min@2.35.2 >/dev/null && tar xzf plotly.js-dist-min-2.35.2.tgz)
    cp "$tmp/package/plotly.min.js" "$dest/plotly.min.js"
    rm -rf "$tmp"
elif command -v curl >/dev/null 2>&1; then
    curl -fsSL https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js -o "$dest/plotly.min.js"
else
    python3 -c "import urllib.request,sys; urllib.request.urlretrieve(sys.argv[1], sys.argv[2])" \
        https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js "$dest/plotly.min.js"
fi
echo "plotly.js saved to $dest/plotly.min.js"
