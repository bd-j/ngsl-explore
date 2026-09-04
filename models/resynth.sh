#!/bin/bash
# Re-synthesize existing converged atmospheres over the wide range that
# carries both the Balmer (3646 A) and Paschen (8206 A) breaks.
cd /Users/bjohnson/Projects/ngsl/models/work
for s in HD194453 HD040573; do
  echo "=== $s ==="
  $ATLAS12/bin/synthe.exe "$s.atm" wlbeg=320 wlend=950 resolu=300000 \
      > "$s.synthe.log" 2>&1 && echo "  ok: $(wc -l < $s.spec) points" \
      || { echo "  FAILED"; tail -3 "$s.synthe.log"; }
done
