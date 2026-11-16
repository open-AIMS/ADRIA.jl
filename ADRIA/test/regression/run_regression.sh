#!/usr/bin/env bash
# Opt-in baseline regression runner (see README.md). Not part of the test suite.
#
#   run_regression.sh ref [REV]                 build + cache the reference for REV (default HEAD)
#   run_regression.sh compare [CHECKOUT] [REV] [--tol X]
#                                               run the driver on CHECKOUT (default: this repo)
#                                               and compare against the cached reference of REV
#                                               (default: the one recorded by the last `ref`)
#   run_regression.sh selftest                  reference-vs-itself passes; perturbations fail
#
# Env: ADRIA_REGRESSION_REF (cache dir, default ~/.cache/adria-regression),
#      ADRIA_REGRESSION_DRIVER (driver script, default driver.jl next to this file),
#      ADRIA_REGRESSION_THREADS (default 4), JULIA (default julia), TMPDIR (scratch).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(git -C "$HERE" rev-parse --show-toplevel)"
CACHE="${ADRIA_REGRESSION_REF:-$HOME/.cache/adria-regression}"
DRIVER="${ADRIA_REGRESSION_DRIVER:-$HERE/driver.jl}"
THREADS="${ADRIA_REGRESSION_THREADS:-4}"
JULIA="${JULIA:-julia}"
SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/adria-regression.XXXXXX")"
WT=""
cleanup() {
    [ -n "$WT" ] && git -C "$REPO" worktree remove --force "$WT" 2>/dev/null || true
    rm -rf "$SCRATCH"
}
trap cleanup EXIT

run_driver() { # checkout outdir
    # Run from an empty dir so no config.toml / Outputs is picked up from (or written to) cwd
    mkdir -p "$SCRATCH/cwd"
    (cd "$SCRATCH/cwd" && "$JULIA" --project="$1/ADRIA" -t "$THREADS" "$DRIVER" "$1/ADRIA" "$2")
}
compare() { # ref cand [extra args]; uses the candidate's environment for Zarr/JSON
    local ref="$1" cand="$2"; shift 2
    "$JULIA" --project="$COMPARE_ENV" "$HERE/compare.jl" "$ref" "$cand" "$@"
}
COMPARE_ENV="${COMPARE_ENV:-$REPO/ADRIA}"

cmd="${1:-}"; shift || true
case "$cmd" in
ref)
    rev="$(git -C "$REPO" rev-parse "${1:-HEAD}")"
    dest="$CACHE/$rev"
    if [ -d "$dest" ]; then echo "reference $rev already cached at $dest"; else
        WT="$SCRATCH/wt"
        git -C "$REPO" worktree add --detach "$WT" "$rev" >/dev/null
        mkdir -p "$CACHE"
        run_driver "$WT" "$SCRATCH/out"
        mv "$SCRATCH/out" "$dest"
        echo "cached reference $rev at $dest"
    fi
    echo "$rev" > "$CACHE/CURRENT"
    ;;
compare)
    cand="${1:-$REPO}"; [ $# -gt 0 ] && shift || true
    rev=""; extra=()
    if [ $# -gt 0 ] && [ "$1" != "--tol" ]; then rev="$1"; shift; fi
    extra=("$@")
    [ -z "$rev" ] && rev="$(cat "$CACHE/CURRENT")"
    [ -d "$CACHE/$rev" ] || { echo "no cached reference for $rev; run: $0 ref $rev" >&2; exit 2; }
    COMPARE_ENV="$(cd "$cand" && pwd)/ADRIA"
    run_driver "$(cd "$cand" && pwd)" "$SCRATCH/cand"
    compare "$CACHE/$rev" "$SCRATCH/cand" "${extra[@]}"
    ;;
selftest)
    rev="$(cat "$CACHE/CURRENT")"
    ref="$CACHE/$rev"
    echo "== reference vs itself (must pass)"
    compare "$ref" "$ref"
    "$JULIA" --project="$COMPARE_ENV" "$HERE/selftest.jl" "$ref" "$HERE/compare.jl"
    ;;
*) sed -n 2,13p "${BASH_SOURCE[0]}"; exit 2 ;;
esac
