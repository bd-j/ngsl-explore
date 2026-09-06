#!/bin/bash
# Resumable fetch of the XSL DR3 tarball. The server drops long connections,
# so a single curl truncates; it supports range requests, so -C - resumes.
cd /Users/bjohnson/Projects/ngsl/data/xsl
URL=http://xsl.u-strasbg.fr/tarball/XSL_DR3_release.tar
WANT=772163584
for i in $(seq 1 60); do
  have=$(stat -f%z XSL_DR3_release.tar 2>/dev/null || echo 0)
  [ "$have" -ge "$WANT" ] && { echo "complete: $have bytes after $((i-1)) resume(s)"; exit 0; }
  echo "attempt $i: have $have / $WANT"
  curl -sSL -C - -A "Mozilla/5.0" --max-time 900 -o XSL_DR3_release.tar "$URL" || true
done
echo "gave up at $(stat -f%z XSL_DR3_release.tar) bytes"
