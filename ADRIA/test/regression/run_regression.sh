#!/usr/bin/env bash
# Opt-in baseline regression runner (see README.md). Not part of the test suite.
#
#   run_regression.sh ref [REV]                 build + cache the reference for REV (default HEAD)
#   run_regression.sh compare [CHECKOUT] [REV] [--tol X]   (--tol may come first)
#                                               run the driver on CHECKOUT (default: this repo)
#                                               and compare against the cached reference of REV
#                                               (default: the one recorded by the last `ref`)
#   run_regression.sh selftest                  reference-vs-itself passes; perturbations fail
#
# Env: ADRIA_REGRESSION_REF (cache dir, default ~/.cache/adria-regression),
#      ADRIA_REGRESSION_DRIVER (driver script, default driver.jl next to this file),
#      ADRIA_REGRESSION_THREADS (default 4), TMPDIR (scratch),
#      JULIA (Julia command, default julia; split on spaces, so a juliaup channel works:
#      JULIA="julia +1.12", or a binary path without spaces).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(git -C "$HERE" rev-parse --show-toplevel)"
CACHE="${ADRIA_REGRESSION_REF:-$HOME/.cache/adria-regression}"
DRIVER="${ADRIA_REGRESSION_DRIVER:-$HERE/driver.jl}"
THREADS="${ADRIA_REGRESSION_THREADS:-4}"
# Word-split on purpose, so JULIA="julia +1.12" selects a juliaup channel
read -r -a JULIA_CMD <<< "${JULIA:-julia}"
die() { echo "ERROR: $*" >&2; exit 2; }
[ ${#JULIA_CMD[@]} -gt 0 ] || die "JULIA is empty or whitespace only; unset it or set e.g. JULIA=\"julia +1.12\""
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/adria-regression.XXXXXX")"
WT=""
cleanup() {
    [ -n "$WT" ] && git -C "$REPO" worktree remove --force "$WT" 2>/dev/null || true
    rm -rf "$SCRATCH"
}
trap cleanup EXIT

# Julia is always started from this empty dir, so no config.toml / Outputs is picked up
# from (or written to) the cwd, and a directory-specific juliaup override applies equally
# to the driver and to `julia_version`
mkdir -p "$SCRATCH/cwd"
run_driver() { # checkout outdir source
    (cd "$SCRATCH/cwd" && "${JULIA_CMD[@]}" --project="$1/ADRIA" -t "$THREADS" "$DRIVER" "$1/ADRIA" "$2")
    write_meta "$2" "$3"
}
driver_sha() {
    local sha
    if command -v sha256sum >/dev/null; then sha="$(sha256sum "$DRIVER" | cut -d' ' -f1)"
    elif command -v shasum >/dev/null; then sha="$(shasum -a 256 "$DRIVER" | cut -d' ' -f1)"
    else die "neither sha256sum nor shasum found"; fi
    [ -n "$sha" ] || die "could not checksum $DRIVER"
    echo "$sha"
}
julia_version() {
    local v
    v="$(cd "$SCRATCH/cwd" && "${JULIA_CMD[@]}" --startup-file=no -e 'print(VERSION)')" || true
    [ -n "$v" ] || die "could not get the Julia version from: ${JULIA_CMD[*]}"
    echo "$v"
}
# outdir source: complete the META the driver wrote (its own Julia VERSION and thread
# count) with the driver and the source; every value must be non-empty
write_meta() {
    [ -f "$1/META" ] || die "the driver wrote no $1/META"
    local sha k
    sha="$(driver_sha)"
    {
        echo "driver=$(basename "$DRIVER")"
        echo "driver_sha256=$sha"
        echo "source=$2"
    } >> "$1/META"
    for k in julia_version threads driver driver_sha256 source; do
        [ -n "$(meta_get "$1" "$k")" ] || die "empty or missing $k in $1/META"
    done
    [ "$(meta_get "$1" threads)" = "$THREADS" ] ||
        die "the driver ran with $(meta_get "$1" threads) threads, not $THREADS"
}
meta_get() { sed -n "s/^$2=//p" "$1/META" 2>/dev/null; } # dir key
compare() { # ref cand [extra args]; uses the candidate's environment for Zarr/JSON
    local ref="$1" cand="$2"; shift 2
    "${JULIA_CMD[@]}" --project="$COMPARE_ENV" "$HERE/compare.jl" "$ref" "$cand" "$@"
}
COMPARE_ENV="${COMPARE_ENV:-$REPO/ADRIA}"
resolve_rev() { # rev-ish -> full commit sha, resolved in this repository
    git -C "$REPO" rev-parse --verify --quiet "$1^{commit}" ||
        die "cannot resolve '$1' to a commit in $REPO"
}
current_rev() { # the rev recorded by the last `ref`
    [ -f "$CACHE/CURRENT" ] ||
        die "no reference recorded in $CACHE (no CURRENT file); run: $0 ref <rev>"
    local rev
    rev="$(tr -d '[:space:]' < "$CACHE/CURRENT")"
    [ -n "$rev" ] || die "$CACHE/CURRENT is empty; run: $0 ref <rev>"
    echo "$rev"
}
require_ref() { # rev: the cached reference folder must exist
    [ -d "$CACHE/$1" ] || die "no cached reference for $1 in $CACHE; run: $0 ref $1"
}

cmd="${1:-}"; shift || true
case "$cmd" in
ref)
    rev="$(resolve_rev "${1:-HEAD}")"
    dest="$CACHE/$rev"
    if [ -d "$dest" ]; then
        echo "reference $rev already cached at $dest"
        if [ ! -f "$dest/META" ]; then
            echo "ERROR: the cached reference has no META (built before META was recorded)," \
                 "so its Julia version, threads and driver cannot be checked." \
                 "Remove $dest to rebuild it, or use another ADRIA_REGRESSION_REF." >&2
            exit 2
        fi
        jv="$(julia_version)"; sha="$(driver_sha)"
        if [ "$(meta_get "$dest" julia_version)" != "$jv" ] ||
           [ "$(meta_get "$dest" driver_sha256)" != "$sha" ] ||
           [ "$(meta_get "$dest" threads)" != "$THREADS" ]; then
            echo "ERROR: the cached reference was built with Julia $(meta_get "$dest" julia_version)," \
                 "$(meta_get "$dest" threads) threads, $(meta_get "$dest" driver);" \
                 "now Julia $jv, $THREADS threads, $(basename "$DRIVER") (or the driver changed)." \
                 "Remove $dest to rebuild it, or use another ADRIA_REGRESSION_REF." >&2
            exit 2
        fi
    else
        WT="$SCRATCH/wt"
        git -C "$REPO" worktree add --detach "$WT" "$rev" >/dev/null
        mkdir -p "$CACHE"
        run_driver "$WT" "$SCRATCH/out" "$rev"
        mv "$SCRATCH/out" "$dest"
        echo "cached reference $rev at $dest"
    fi
    echo "$rev" > "$CACHE/CURRENT"
    ;;
compare)
    # Positional CHECKOUT then REV, each optional; options (--tol X) anywhere
    cand=""; rev=""; extra=()
    while [ $# -gt 0 ]; do
        case "$1" in
        --tol) [ $# -ge 2 ] || { echo "--tol needs a value" >&2; exit 2; }
               extra+=("$1" "$2"); shift 2 ;;
        --*) echo "unknown option $1" >&2; exit 2 ;;
        *) if [ -z "$cand" ]; then cand="$1"; elif [ -z "$rev" ]; then rev="$1"; else
               echo "too many arguments: $1" >&2; exit 2; fi
           shift ;;
        esac
    done
    cand="${cand:-$REPO}"
    if [ -n "$rev" ]; then rev="$(resolve_rev "$rev")"; else rev="$(current_rev)"; fi
    require_ref "$rev"
    [ -d "$cand/ADRIA" ] || die "$cand is not an ADRIA checkout (no ADRIA/ folder)"
    cand_abs="$(cd "$cand" && pwd)"
    cand_rev="$(git -C "$cand_abs" rev-parse HEAD 2>/dev/null || echo unknown)"
    [ -n "$(git -C "$cand_abs" status --porcelain --untracked-files=no 2>/dev/null)" ] && cand_rev="$cand_rev+dirty"
    COMPARE_ENV="$cand_abs/ADRIA"
    run_driver "$cand_abs" "$SCRATCH/cand" "$cand_abs at $cand_rev"
    compare "$CACHE/$rev" "$SCRATCH/cand" ${extra[@]+"${extra[@]}"}
    ;;
selftest)
    rev="$(current_rev)"
    require_ref "$rev"
    ref="$CACHE/$rev"
    echo "== reference vs itself (must pass)"
    compare "$ref" "$ref"
    "${JULIA_CMD[@]}" --project="$COMPARE_ENV" "$HERE/selftest.jl" "$ref" "$HERE/compare.jl"
    ;;
*) sed -n 2,15p "${BASH_SOURCE[0]}"; exit 2 ;;
esac
