#!/usr/bin/env bash
# Throwaway diagnostic: print cache/env state and probe package loading.
set +e
echo "##### DIAG label=$1 host=$(hostname)"
T=$GITHUB_WORKSPACE/ADRIA/test/.testenv
sha256sum $T/Project.toml $T/Manifest.toml; stat -c '%y %n' $T/Project.toml $T/Manifest.toml
lscpu | grep -i 'model name'
env | grep -i JULIA
C=~/.julia/compiled/v1.12
ls -la $C | head -8
echo "files: $(find ~/.julia/compiled -type f | wc -l)"; du -sh ~/.julia/compiled
echo "ji/so counts: $(find $C -name '*.ji' | wc -l) / $(find $C -name '*.so' | wc -l)"
find $C -type f | sort | sed "s#$HOME##" | sha256sum
ls -la ~/.julia/dev ~/.julia/environments 2>&1 | head
cat > /tmp/uuids.jl <<'JL'
using TOML
m = TOML.parsefile(joinpath(ARGS[1], "Manifest.toml"))["deps"]
for n in ["RangeArrays","ExprTools","Reexport"]
    println(n, " ", m[n][1]["uuid"])
end
JL
UU=$(julia --startup-file=no /tmp/uuids.jl $T)
echo "$UU"
for n in RangeArrays ExprTools Reexport; do echo "--- $C/$n"; ls -la --time-style=full-iso $C/$n; done
cd $GITHUB_WORKSPACE
BASE="--project=$T --check-bounds=yes --startup-file=no"
COV="--code-coverage=@$GITHUB_WORKSPACE/ADRIA"
cat > /tmp/probe.jl <<'JL'
using Base: PkgId, UUID
n = "RangeArrays"
u = UUID(split(first(filter(l->startswith(l,n), readlines(`cat /tmp/uuids.txt`))))[2])
@time Base.require(PkgId(u, n))
JL
echo "$UU" > /tmp/uuids.txt
probe() {
  name=$1; shift
  echo "=== PROBE $name: $*"
  start=$(date +%s)
  out=$(env "${ENVX[@]}" JULIA_DEBUG=loading julia "$@" /tmp/probe.jl 2>&1)
  echo "elapsed $(( $(date +%s)-start ))s"
  echo "$out" | grep -E "Rejecting|Loading object cache|Precompiling|Skipping|stale|seconds" | head -15
  echo "rejections total: $(echo "$out" | grep -c Rejecting)"
  echo "$out" | grep -q Precompiling && { echo "PROBE $name TRIGGERED PRECOMPILE"; echo "$out" | grep -A3 -E "Rejecting cache file.*RangeArrays|Rejecting.*RangeArrays" | head -10; return 1; }
  return 0
}
ENVX=(JULIA_NUM_THREADS=auto)
probe A $BASE $COV --warn-overwrite=yes --depwarn=yes || exit 0
ENVX=(X=1)
probe B $BASE $COV --warn-overwrite=yes --depwarn=yes || exit 0
probe C $BASE --warn-overwrite=yes --depwarn=yes || exit 0
probe D $BASE $COV || exit 0
probe E --project=$T --startup-file=no || exit 0
