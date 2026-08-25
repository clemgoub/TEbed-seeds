#!/usr/bin/env bash
# =============================================================================
# overassembly_seqtest.sh -- sequence proof/disproof that a library consensus
#                            is an OVER-ASSEMBLY of k tandem units.
#
# Companion to tools/detect_overassembly.py, which flags candidates from BED16
# coordinates alone.  The coordinate detector cannot distinguish "k tandem
# units" from "a chimera cut at a commensurate position" (measured decoy-null
# FDR ceiling ~0.80 for the integer-k claim).  Only the sequence can.
# Background: RESUME.md section 5, PIPELINE_FINDINGS.md F4.
#
# WHAT IT RUNS
#   1. SELF-ALIGNMENT   blastn -word_size 7 -dust no -strand both, query=subject.
#      Every off-diagonal HSP is reported with its OFFSET (qstart - sstart) and
#      that offset expressed as a fraction of the consensus length L and as
#      offset*k/L for k = 2..5.
#        * k tandem units  => a LONG SAME-STRAND off-diagonal HSP at offset
#          ~ L/k (and weaker ones at 2L/k ... ), ideally covering ~(k-1)*L/k.
#        * a palindrome / MITE with terminal inverted repeats => a SHORT
#          OPPOSITE-STRAND HSP pinned to the two ENDS (qstart ~ 1 and
#          send ~ L), and nothing else.  This is the alternative hypothesis
#          F4 tested and rejected on coordinates; the self-alignment settles it.
#        * an internal satellite => same-strand off-diagonals at a period that
#          is NOT L/k for small k (measured: rm2 rnd-3_family-344, L = 632,
#          period ~45-66 bp, i.e. k ~ 10-14).
#   2. TIR MEASUREMENT  the 5' end vs the reverse complement of the 3' end,
#      ungapped, reporting the length at which mismatches start accumulating.
#   3. CROSS-ALIGNMENT  the family against every other library given with -c.
#      Over-assembly predicts a SHORT consensus (~L/k) from another tool
#      matching the long one k times, at the block boundaries.
#   4. VERDICT per k, from the measured HSPs.
#
# USAGE
#   tools/overassembly_seqtest.sh -l LIB.fa -f FAMILY [options]
#
#     -l LIB.fa      library FASTA containing the family                (required)
#     -f FAMILY      family name; matched against the FASTA id, with or
#                    without a '#Class/Subclass' suffix                 (required)
#     -c OTHER.fa    cross-align library; repeatable.  Give the libraries of
#                    the OTHER tools in the same cluster.
#     -k 2,3,4,5     candidate unit counts to score            (default 2,3,4,5)
#     -o OUTDIR      output directory        (default overassembly_seqtest/FAMILY)
#     -m MINLEN      minimum HSP length to report                    (default 12)
#     -t TOL         offset tolerance as a fraction of L/k          (default 0.15)
#     -h             this help
#
# TOMORROW, when the REPET library lands, the whole cluster-62 test is:
#
#   L=libs/repet_consensi.fa                      # the file the user provides
#   for F in Gnig_TEdenovoGr-B-G1303-Map20 Gnig_TEdenovoGr-B-G1473-Map8 \
#            Gnig_TEdenovoGr-B-G1039-Map3  Gnig_TEdenovoGr-B-G1568-Map3 ; do
#     tools/overassembly_seqtest.sh -l "$L" -f "$F" \
#       -c libs/GCA_951799975.1-families.fa \
#       -c libs/fGobNig_all_rounds_classified.fasta \
#       -o runs/GCA_951799975.1/seqtest/"$F"
#   done
#
# All four are flagged k=3 with unit_len 254.7-257.0 bp by the coordinate
# detector, and their rm2/pantera cluster-mates are 264-270 bp.  The prediction
# is a same-strand self-alignment HSP of ~510 bp at offset ~255.
#
# REQUIRES: blastn (NCBI BLAST+), awk, python3 (for the TIR scan only).
# =============================================================================
set -euo pipefail

LIB=""; FAM=""; KS="2,3,4,5"; OUT=""; MINLEN=12; TOL=0.15
CROSS=()

usage() { sed -n '2,63p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while getopts ":l:f:c:k:o:m:t:h" opt; do
  case "$opt" in
    l) LIB="$OPTARG" ;;
    f) FAM="$OPTARG" ;;
    c) CROSS+=("$OPTARG") ;;
    k) KS="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    m) MINLEN="$OPTARG" ;;
    t) TOL="$OPTARG" ;;
    h) usage 0 ;;
    \?) echo "unknown option -$OPTARG" >&2; usage 1 ;;
    :)  echo "option -$OPTARG needs an argument" >&2; usage 1 ;;
  esac
done

[ -n "$LIB" ] && [ -n "$FAM" ] || { echo "ERROR: -l and -f are required" >&2; usage 1; }
[ -f "$LIB" ] || { echo "ERROR: no such library: $LIB" >&2; exit 1; }
command -v blastn >/dev/null || { echo "ERROR: blastn not on PATH" >&2; exit 1; }

OUT="${OUT:-overassembly_seqtest/${FAM//\//_}}"
mkdir -p "$OUT"
Q="$OUT/query.fa"

# ---------------------------------------------------------------- 0. extract
# Match the id exactly, or the id up to the first '#' (RepeatMasker style),
# or the first whitespace-delimited token.
awk -v fam="$FAM" '
  /^>/ { id=substr($0,2); split(id,w," "); t=w[1]; split(t,h,"#");
         keep = (t==fam || h[1]==fam || id==fam) ? 1 : 0;
         if (keep) print ">" t; next }
  keep { print }
' "$LIB" > "$Q"

if [ ! -s "$Q" ]; then
  echo "ERROR: '$FAM' not found in $LIB" >&2
  echo "       ids look like:" >&2
  grep '^>' "$LIB" | head -3 >&2
  exit 1
fi

L=$(awk '!/^>/ {n+=length($0)} END {print n}' "$Q")
NSEQ=$(grep -c '^>' "$Q")
[ "$NSEQ" -eq 1 ] || { echo "ERROR: '$FAM' matched $NSEQ sequences" >&2; exit 1; }

echo "================================================================"
echo "family    : $FAM"
echo "library   : $LIB"
echo "length L  : $L bp"
echo "outdir    : $OUT"
echo "================================================================"

FMT='6 qstart qend sstart send sstrand length pident mismatch gaps evalue bitscore'

# ------------------------------------------------------- 1. self-alignment
blastn -query "$Q" -subject "$Q" -word_size 7 -dust no -strand both \
       -evalue 10 -outfmt "$FMT" > "$OUT/self.tsv"

echo
echo "--- 1. SELF-ALIGNMENT (off-diagonal HSPs >= ${MINLEN} bp) ---"
echo "    offset = qstart - sstart ; offset*k/L should be ~1 for a k-unit array"
awk -F'\t' -v L="$L" -v m="$MINLEN" 'BEGIN{
    printf "%-18s %-18s %-6s %6s %7s %9s %8s %8s %8s %8s\n",
           "query","subject","strand","len","pident","offset","off/L","off*2/L","off*3/L","off*4/L"
  }
  {
    if ($6 >= L*0.98 && $1==$3) next          # the main diagonal
    if ($6 < m) next
    off = $1 - $3
    printf "%-18s %-18s %-6s %6d %7.1f %9d %8.3f %8.2f %8.2f %8.2f\n",
           $1"-"$2, $3"-"$4, $5, $6, $7, off, off/L, off*2/L, off*3/L, off*4/L
  }' "$OUT/self.tsv" | head -40

# ------------------------------------------------------------- 2. TIR scan
echo
echo "--- 2. TERMINAL INVERTED REPEAT scan (5' end vs revcomp of 3' end) ---"
python3 - "$Q" <<'PYEOF'
import sys
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">")).upper()
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = seq.translate(comp)[::-1]
n = min(120, len(seq) // 2)
mm = 0; last_ok = 0; runs = []
for i in range(n):
    if seq[i] != rc[i]:
        mm += 1
    else:
        last_ok = i + 1
    runs.append((i + 1, mm))
# longest prefix with <=1, <=2, <=3 mismatches
for tol in (0, 1, 2, 3):
    best = max((p for p, m in runs if m <= tol), default=0)
    print(f"    longest terminal inverted repeat with <= {tol} mismatch(es): {best} bp")
print(f"    5' first 30 : {seq[:30]}")
print(f"    3' rc  30   : {rc[:30]}")
print("    (a TIR element gives a short perfect-ish match here AND a short")
print("     opposite-strand terminal HSP in section 1 -- neither is evidence")
print("     of tandem over-assembly.)")
PYEOF

# ------------------------------------------------------ 3. cross-alignment
if [ "${#CROSS[@]}" -gt 0 ]; then
  echo
  echo "--- 3. CROSS-ALIGNMENT against other tools' libraries ---"
  for C in "${CROSS[@]}"; do
    [ -f "$C" ] || { echo "    (skipping missing $C)"; continue; }
    B="$OUT/cross_$(basename "$C" | tr -c 'A-Za-z0-9._-' '_').tsv"
    blastn -query "$Q" -subject "$C" -word_size 11 -dust no -strand both \
           -evalue 1e-5 -outfmt "6 sseqid slen qstart qend sstart send sstrand length pident bitscore" \
           > "$B" || true
    echo "    vs $(basename "$C")  ->  $B"
    # subjects that hit the query MORE THAN ONCE at >=100 bp are the signal:
    # one short entry matching the long consensus k times.
    awk -F'\t' -v L="$L" '$8>=100 {n[$1]++; len[$1]=$2;
        blk[$1]=blk[$1] sprintf(" q%d-%d(%s)",$3,$4,$7)}
      END{
        printf "      %-42s %7s %6s  %s\n","subject","slen","nHSP","query blocks"
        k=0
        for (s in n) { k++
          printf "      %-42s %7d %6d %s\n", substr(s,1,42), len[s], n[s], blk[s] }
        if (k==0) print "      (no subject matched at >=100 bp)"
      }' "$B" | sort -k3,3nr | head -20
    echo "      ^ over-assembly predicts a subject of length ~L/k hitting the"
    echo "        query k times, at consecutive non-overlapping query blocks."
  done
fi

# ------------------------------------------------------------- 4. verdict
echo
echo "--- 4. VERDICT per k (same-strand off-diagonal at offset ~ L/k) ---"
echo "    a k-unit array should give a same-strand HSP of ~ (k-1)*L/k bp"
echo "    at offset ~ L/k, i.e. tandem_cover ~ 1.00"
for K in ${KS//,/ }; do
  awk -F'\t' -v L="$L" -v k="$K" -v tol="$TOL" -v m="$MINLEN" '
    BEGIN { P = L/k; want = (k-1)*P; best = 0; bo = 0 }
    {
      if ($5 != "plus") next
      if ($6 < m) next
      off = $1 - $3
      if (off <= 0) next
      if (off > L*0.98) next
      for (j = 1; j < k; j++) {
        tgt = j*P
        if (off >= tgt*(1-tol) && off <= tgt*(1+tol)) {
          if ($6 > best) { best = $6; bo = off }
        }
      }
    }
    END {
      cov = (want > 0) ? best/want : 0
      verdict = (cov >= 0.50) ? "SUPPORTS k=" k \
              : (cov >= 0.20) ? "weak/partial support for k=" k \
              : "NO support for k=" k
      printf "    k=%-2d unit=%7.1f  best same-strand HSP at offset ~L/k: %5d bp (offset %d)  tandem_cover=%.2f  %s\n",
             k, P, best, bo, cov, verdict
    }' "$OUT/self.tsv"
done

echo
echo "raw HSP tables in $OUT/"
