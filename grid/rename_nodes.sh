#!/usr/bin/env bash
#
# Rename grid node files from the old naming to the C3K v2.3 convention.
#
#     t10000g4.00m-0.10.spec  ->  at12_feh-0.10_afe+0.0_t10000g4.00.spec
#
# Everything after the node name is preserved verbatim, so .atm, .spec,
# .flux, .iter, .atlas.log, .synthe.log, .FAILED and _stepN.atm all come
# along without being enumerated here.  Files that do not match the old
# pattern -- including any C3K starting atmospheres sitting in the same
# directory, which are already in the new convention -- are left alone.
#
# Nothing is moved until the whole plan is known to be clean: if any
# destination already exists the script reports every conflict and exits
# without touching a file.  Re-run with -f to overwrite.
#
# Usage: rename_nodes.sh [-n] [-f] [DIR]
#   -n   dry run: print the plan, move nothing
#   -f   force: overwrite existing destinations
#   DIR  directory to convert (default: the current directory)

set -euo pipefail

dry=0
force=0
while getopts ':nfh' opt; do
    case "$opt" in
        n) dry=1 ;;
        f) force=1 ;;
        h) sed -n '2,25p' "$0"; exit 0 ;;
        *) echo "unknown option -$OPTARG (try -h)" >&2; exit 2 ;;
    esac
done
shift $((OPTIND - 1))
dir="${1:-.}"

[ -d "$dir" ] || { echo "ERROR: not a directory: $dir" >&2; exit 2; }

# t<TTTTT>g<G.GG>m<+M.MM><anything>
re='^t([0-9]{5})g(-?[0-9]\.[0-9]{2})m([+-][0-9]\.[0-9]{2})(.*)$'

srcs=()
dsts=()
conflicts=()

shopt -s nullglob
for path in "$dir"/*; do
    [ -f "$path" ] || continue
    base="${path##*/}"
    [[ "$base" =~ $re ]] || continue
    teff="${BASH_REMATCH[1]}"
    logg="${BASH_REMATCH[2]}"
    feh="${BASH_REMATCH[3]}"
    rest="${BASH_REMATCH[4]}"
    new="at12_feh${feh}_afe+0.0_t${teff}g${logg}${rest}"
    srcs+=("$base")
    dsts+=("$new")
    if [ -e "$dir/$new" ]; then
        conflicts+=("$base -> $new")
    fi
done
shopt -u nullglob

n=${#srcs[@]}
if [ "$n" -eq 0 ]; then
    echo "no files in $dir match the old node naming; nothing to do"
    exit 0
fi

if [ "${#conflicts[@]}" -gt 0 ] && [ "$force" -eq 0 ] && [ "$dry" -eq 0 ]; then
    echo "ERROR: ${#conflicts[@]} destination(s) already exist; nothing renamed." >&2
    printf '  %s\n' "${conflicts[@]}" >&2
    echo >&2
    echo "These are usually the copied C3K starting atmospheres, which already" >&2
    echo "carry the new name.  The node's own output is the file worth keeping," >&2
    echo "so re-run with -f to overwrite them -- or move the .atm starts aside" >&2
    echo "first if you want to keep them." >&2
    exit 1
fi

echo "$dir: $n file(s) to rename"
if [ "$dry" -eq 1 ]; then
    for ((i = 0; i < n; i++)); do
        printf '  %s -> %s\n' "${srcs[$i]}" "${dsts[$i]}"
    done
    if [ "${#conflicts[@]}" -gt 0 ]; then
        echo "  ${#conflicts[@]} of these would overwrite an existing file:"
        printf '    %s\n' "${conflicts[@]}"
        echo "  the real run refuses unless given -f"
    fi
    echo "(dry run; nothing moved)"
    exit 0
fi

moved=0
for ((i = 0; i < n; i++)); do
    mv -f -- "$dir/${srcs[$i]}" "$dir/${dsts[$i]}"
    moved=$((moved + 1))
done

# Confirm by counting what is actually on disk, not by the loop finishing.
shopt -s nullglob
left=("$dir"/t[0-9][0-9][0-9][0-9][0-9]g*m[+-]*)
now=("$dir"/at12_feh*_afe+0.0_t*)
shopt -u nullglob
echo "renamed $moved; ${#left[@]} old-style file(s) left, ${#now[@]} new-style file(s) present"
[ "${#left[@]}" -eq 0 ] || { echo "ERROR: old-style files remain" >&2; exit 1; }
