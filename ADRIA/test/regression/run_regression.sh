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
    (cd "$SCRATCH/cwd" && "${JULIA_CMD[@]}" --project="$1/ADRIA" -t "$THREADS" "$DRIVER" "$1/ADRIA" "$2")
    write_meta "$2" "$3"
}
driver_sha() { sha256sum "$DRIVER" | cut -d' ' -f1; }
julia_version() { "${JULIA_CMD[@]}" --startup-file=no -e 'print(VERSION)'; }
write_meta() { # outdir source: what compare.jl checks (Julia version, threads, driver)
    {
        echo "julia_version=$(julia_version)"
        echo "threads=$THREADS"
        echo "driver=$(basename "$DRIVER")"
        echo "driver_sha256=$(driver_sha)"
        echo "source=$2"
    } > "$1/META"
}
meta_get() { sed -n "s/^$2=//p" "$1/META" 2>/dev/null; } # dir key
compare() { # ref cand [extra args]; uses the candidate's environment for Zarr/JSON
    local ref="$1" cand="$2"; shift 2
    "${JULIA_CMD[@]}" --project="$COMPARE_ENV" "$HERE/compare.jl" "$ref" "$cand" "$@"
}
COMPARE_ENV="${COMPARE_ENV:-$REPO/ADRIA}"

cmd="${1:-}"; shift || true
case "$cmd" in
ref)
    rev="$(git -C "$REPO" rev-parse "${1:-HEAD}")"
    dest="$CACHE/$rev"
    if [ -d "$dest" ]; then
        echo "reference $rev already cached at $dest"
        if [ ! -f "$dest/META" ]; then
            echo "WARNING: it has no META (built before META was recorded); rebuild it to check versions" >&2
        elif [ "$(meta_get "$dest" julia_version)" != "$(julia_version)" ] ||
             [ "$(meta_get "$dest" driver_sha256)" != "$(driver_sha)" ] ||
             [ "$(meta_get "$dest" threads)" != "$THREADS" ]; then
            echo "ERROR: the cached reference was built with Julia $(meta_get "$dest" julia_version)," \
                 "$(meta_get "$dest" threads) threads, $(meta_get "$dest" driver);" \
                 "now Julia $(julia_version), $THREADS threads, $(basename "$DRIVER") (or the driver changed)." \
                 "Remove $dest or use another ADRIA_REGRESSION_REF." >&2
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
    [ -z "$rev" ] && rev="$(cat "$CACHE/CURRENT")"
    [ -d "$CACHE/$rev" ] || { echo "no cached reference for $rev; run: $0 ref $rev" >&2; exit 2; }
    cand_abs="$(cd "$cand" && pwd)"
    cand_rev="$(git -C "$cand_abs" rev-parse HEAD 2>/dev/null || echo unknown)"
    [ -n "$(git -C "$cand_abs" status --porcelain --untracked-files=no 2>/dev/null)" ] && cand_rev="$cand_rev+dirty"
    COMPARE_ENV="$cand_abs/ADRIA"
    run_driver "$cand_abs" "$SCRATCH/cand" "$cand_abs at $cand_rev"
    compare "$CACHE/$rev" "$SCRATCH/cand" ${extra[@]+"${extra[@]}"}
    ;;
selftest)
    rev="$(cat "$CACHE/CURRENT")"
    ref="$CACHE/$rev"
    echo "== reference vs itself (must pass)"
    compare "$ref" "$ref"
    "${JULIA_CMD[@]}" --project="$COMPARE_ENV" "$HERE/selftest.jl" "$ref" "$HERE/compare.jl"
    ;;
*) sed -n 2,15p "${BASH_SOURCE[0]}"; exit 2 ;;
esac
