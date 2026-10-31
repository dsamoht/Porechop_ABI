#!/bin/bash
# Parity runs for the Rust branch against upstream Porechop_ABI 0.5.1 (checked out in ./upstream).
# Every Rust output is compared with upstream on this runner and with the hashes from the macOS reference run.
set -u
W=$PWD; DATA=$W/data; RUNS=$W/runs; mkdir -p $DATA $RUNS
RES=$RUNS/results.tsv; : > $RES
NT=$(python -c 'import os; print(os.cpu_count())')
PLAT="$(uname -s)-$(uname -m)"
now() { python -c 'import time; print(time.time())'; }
sha() { shasum "$1" 2>/dev/null | cut -c1-40; }
log() { echo -e "$*" | tee -a $RES; }

log "platform\t$PLAT\tcpus=$NT\t$(python --version 2>&1)\t$(rustc --version)\t$( (c++ --version || g++ --version) 2>/dev/null | head -1)"

# ---- data (ENA) ----
while read -r acc url; do
  f=$DATA/$acc.fastq.gz
  for i in 1 2 3 4 5 6 7 8 9 10; do
    gzip -t $f 2>/dev/null && break
    curl -sS -C - --retry 5 --retry-all-errors -o $f "$url" || true
  done
  gzip -t $f || { log "data\t$acc\tDOWNLOAD FAILED"; exit 1; }
done < $W/rust/.github/parity/datasets.tsv

# ---- builds ----
python -m venv $W/venv_up && $W/venv_up/bin/pip install -q networkx setuptools
if [ "$(uname -s)" = Darwin ]; then
  # SeqAn 2.4 uses std::unary_function, removed from libc++ in C++17 mode
  export CXXFLAGS="-D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION"
fi
( cd $W/upstream && $W/venv_up/bin/pip install -v --no-build-isolation . ) > $RUNS/build_up.log 2>&1
unset CXXFLAGS
UP_OK=0; $W/venv_up/bin/porechop_abi --version >/dev/null 2>&1 && UP_OK=1
log "build\tupstream\t$([ $UP_OK = 1 ] && echo OK || echo FAILED)"
[ $UP_OK = 1 ] || grep -E "error" $RUNS/build_up.log | sort -u | head -20 | tee -a $RES

python -m venv $W/venv_rs && $W/venv_rs/bin/pip install -q networkx maturin
( cd $W/rust && $W/venv_rs/bin/maturin build --release -i $W/venv_rs/bin/python -o $W/wheels && $W/venv_rs/bin/pip install -q $W/wheels/*.whl ) > $RUNS/build_rs.log 2>&1
$W/venv_rs/bin/porechop_abi --version >/dev/null 2>&1 && log "build\trust\tOK" || { log "build\trust\tFAILED"; tail -30 $RUNS/build_rs.log; exit 1; }

# ---- GIL check: upstream thread scaling on 1,000 reads ----
UP_T=$NT
if [ $UP_OK = 1 ]; then
  d=$RUNS/gil; mkdir -p $d; gzip -dc $DATA/DRR253169.fastq.gz | head -4000 > $d/s.fastq
  for v in up rs; do for t in 1 $NT; do
    t0=$(now); $W/venv_$v/bin/porechop_abi -i $d/s.fastq -o $d/$v$t.fastq -t $t > /dev/null 2>&1; t1=$(now)
    eval "T_${v}_$t=$(python -c "print(round($t1-$t0,1))")"
  done; done
  log "gil\tup\tt1=${T_up_1}s\tt$NT=$(eval echo \$T_up_$NT)s\trs\tt1=${T_rs_1}s\tt$NT=$(eval echo \$T_rs_$NT)s\tsame=$( [ $(cat $d/*.fastq | wc -c) -gt 0 ] && [ $(for f in $d/up*.fastq $d/rs*.fastq; do sha $f; done | sort -u | wc -l) = 1 ] && echo yes || echo no)"
  python -c "import sys; sys.exit(0 if $T_up_1 < $(eval echo \$T_up_$NT) else 1)" && UP_T=1
fi
log "threads\tupstream=$UP_T\trust=$NT"

# ---- parity runs ----
run() {  # run MODE NAME EXPECTED_SHA ARGS...
  local mode=$1 name=$2 exp=$3; shift 3
  local d=$RUNS/$mode/$name; mkdir -p $d/up $d/rs
  local tu=- hu=-
  if [ $UP_OK = 1 ] && [ "$mode" != seeded ]; then
    t0=$(now); ( cd $d/up && $W/venv_up/bin/porechop_abi "$@" -o out.fastq -t $UP_T -tmp tmp > log.txt 2>&1 < /dev/null ); t1=$(now)
    tu=$(python -c "print(round($t1-$t0,1))"); hu=$(sha $d/up/out.fastq)
  fi
  t0=$(now); ( cd $d/rs && $W/venv_rs/bin/porechop_abi "$@" -o out.fastq -t $NT -tmp tmp > log.txt 2>&1 < /dev/null ); t1=$(now)
  local tr=$(python -c "print(round($t1-$t0,1))") hr=$(sha $d/rs/out.fastq)
  local vs_up=$([ "$hu" = - ] && echo n/a || ([ "$hu" = "$hr" ] && echo IDENTICAL || echo DIFFERENT))
  local vs_mac=$([ "$hr" = "$exp" ] && echo IDENTICAL || echo DIFFERENT)
  log "$mode\t$name\tup=${tu}s\trs=${tr}s\trust_vs_upstream=$vs_up\trust_vs_macos=$vs_mac\t${hr:0:12}"
  rm -f $d/up/out.fastq $d/rs/out.fastq
}

printf 'sn=2000\nk=16\nsl=100\nlim=500\nlc=1.0\n' > $RUNS/sampled.config
while read -r acc url; do
  run db $acc $(awk -v a=$acc '$1=="db" && $2==a {print $3}' $W/rust/.github/parity/expected.tsv) -i $DATA/$acc.fastq.gz
done < $W/rust/.github/parity/datasets.tsv
while read -r acc url; do
  run abi $acc $(awk -v a=$acc '$1=="abi" && $2==a {print $3}' $W/rust/.github/parity/expected.tsv) -abi -i $DATA/$acc.fastq.gz
done < $W/rust/.github/parity/datasets.tsv
for s in 1 3; do
  run seeded ERR6059982_seed$s $(awk -v s=seed$s '$1=="seeded" && $2==s {print $3}' $W/rust/.github/parity/expected.tsv) \
    -abi -abc $RUNS/sampled.config --seed $s -i $DATA/ERR6059982.fastq.gz
done

# ---- repo test suite against each installed build ----
for v in up rs; do
  [ $v = up ] && [ $UP_OK = 0 ] && continue
  d=$RUNS/unittest/$v; mkdir -p $d; cp -r $W/rust/test $d/test; cp $W/rust/porechop-abi-runner.py $d/
  ( cd $d && PATH=$W/venv_$v/bin:$PATH python -m unittest > result.txt 2>&1 ); rc=$?
  log "unittest\t$v\texit=$rc\t$(grep -E '^(Ran|OK|FAILED)' $d/result.txt | paste -sd' ' -)"
done
