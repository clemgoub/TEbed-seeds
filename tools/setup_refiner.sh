#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# setup_refiner.sh — install Dfam's Refiner (RepeatModeler 2.x) standalone
#                    on macOS ARM (Apple silicon), no root, no compilation.
#
# WHAT THIS BUYS YOU
#   Refiner builds a TE consensus + a Stockholm seed alignment from a FASTA of
#   genomic copies.  It is the reference implementation used by RepeatModeler
#   and it is what step 4.7 of RESUME.md compares MAFFT against.
#
# THE "STANDALONE" CLAIM — VERIFIED, WITH TWO CAVEATS (measured 2026-08-24)
#   Refiner is a Perl script and needs NO compilation, but it is NOT
#   self-contained.  It needs three things beyond RepeatModeler itself:
#
#   1. RepeatMasker's Perl modules.  Refiner line 99 does
#        use lib $RepModelConfig::configuration->{'REPEATMASKER_DIR'}->{'value'};
#      then `use SearchResult / SearchResultCollection / WUBlastSearchEngine /
#      NCBIBlastSearchEngine / SeqDBI / SimpleBatcher / FastaDB;`  (confirmed by
#      reading Refiner lines 98-111).  All six live in RepeatMasker, are pure
#      Perl, and need no build.  MultAln.pm, NeedlemanWunschGotohAlgorithm.pm
#      and SequenceSimilarityMatrix.pm ship with RepeatModeler itself.
#
#      CAVEAT A — the REPEATMASKER_DIR *environment variable does not work* for
#      this.  `use lib` runs at COMPILE time, while RepModelConfig's
#      resolveConfiguration() (which applies environment_override) runs at RUN
#      time, so `use lib` sees the hardcoded default
#      /usr/local/RepeatMasker-4.2.4-Dfam-4.0_RB and the compile dies with
#      "Can't locate SearchResultCollection.pm in @INC".  Verified.
#      The fix used here is PERL5LIB, which Perl honours at compile time.
#
#   2. rmblastn — NCBI BLAST+ patched by Dfam/ISB.  STOCK blastn CANNOT
#      SUBSTITUTE.  Refiner's generated command line (captured with
#      `Refiner -debug 9`) is:
#        rmblastn -num_alignments 9999999 -db X -query X -gapopen 20 \
#          -gapextend 5 -complexity_adjust -word_size 7 -xdrop_ungap 300 \
#          -xdrop_gap_final 150 -xdrop_gap 75 -min_raw_gapped_score 150 \
#          -dust no -outfmt="6 score perc_sub perc_query_gap perc_db_gap \
#          qseqid qstart qend qlen sstrand sseqid sstart send slen kdiv \
#          cpg_kdiv transi transv cpg_sites qseq sseq" -num_threads 4 \
#          -matrix comparison.matrix
#      Against Homebrew blastn 2.16.0+ that dies immediately with
#        Error: Unknown argument: "complexity_adjust"
#      and -matrix and -mask_level (used in the one-vs-all pass) are likewise
#      rejected: stock blastn has no arbitrary nucleotide substitution-matrix
#      support at all.  Worse, the custom outfmt fields kdiv / cpg_kdiv /
#      transi / transv / cpg_sites are silently DROPPED by stock blastn rather
#      than erroring, so even a hypothetical flag-stripping shim would feed the
#      Perl parser truncated rows.  Substituting blastn was tested end to end:
#      Refiner exits with "ERROR from search engine" and
#      "No consensus could be derived from the input sequences".
#
#      CAVEAT B — so a real rmblastn is mandatory.  Fortunately Dfam ships a
#      prebuilt arm64 macOS binary (rmblast-2.17.1+-arm64-macosx.tar.gz,
#      ~177 MB), which this script downloads.  No NCBI toolkit build needed.
#
#   3. Nothing else.  MAFFT, TRF, RECON, RepeatScout, CD-HIT, LTR_retriever,
#      genometools, NINJA, UCSC tools are needed by the full RepeatModeler
#      pipeline, NOT by Refiner.  Scoring matrices come from
#      RepeatModeler/Matrices, which is in the RepeatModeler clone.
#
# WHAT IT DOWNLOADS
#   - https://github.com/Dfam-consortium/RepeatModeler   (shallow git clone)
#   - https://github.com/Dfam-consortium/RepeatMasker    (shallow git clone)
#   - https://www.repeatmasker.org/rmblast/rmblast-2.17.1+-arm64-macosx.tar.gz
#
# WHAT IT ASSUMES
#   - macOS on Apple silicon (arm64).  For Intel macOS or Linux, change
#     RMBLAST_TARBALL below to the matching build on repeatmasker.org.
#   - system perl (tested: /usr/bin/perl 5.34.1) — only core modules are used.
#   - git, curl, tar on PATH.  ~1.2 GB free in $REFINER_PREFIX.
#
# USAGE
#   bash tools/setup_refiner.sh                 # install into ./vendor/refiner
#   REFINER_PREFIX=/some/dir bash tools/setup_refiner.sh
#   Afterwards:
#     source <prefix>/refiner_env.sh   &&  Refiner copies.fa
#   or just:
#     <prefix>/bin/Refiner copies.fa
#
# OUTPUTS OF A REFINER RUN (given input copies.fa)
#   copies.fa.refiner_cons  — FASTA consensus
#   copies.fa.refiner.stk   — Stockholm seed alignment ('.' gaps, #=GC RF,
#                             ids like c1:1-247_+).  NOTE: it carries only
#                             #=GF ID/CC/BM/SQ, so `stk lint` reports exactly
#                             four errors — missing DE, AU, TP, OC — which the
#                             TEbed-seeds pipeline must add itself.
# ---------------------------------------------------------------------------
set -euo pipefail

REFINER_PREFIX="${REFINER_PREFIX:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/vendor/refiner}"
RMBLAST_VER="2.17.1"
RMBLAST_TARBALL="rmblast-${RMBLAST_VER}+-arm64-macosx.tar.gz"
RMBLAST_URL="https://www.repeatmasker.org/rmblast/${RMBLAST_TARBALL}"

log() { printf '[setup_refiner] %s\n' "$*" >&2; }

if [[ "$(uname -s)" != "Darwin" || "$(uname -m)" != "arm64" ]]; then
  log "WARNING: this script hardcodes the arm64 macOS rmblast build."
  log "         uname is $(uname -s)/$(uname -m).  Edit RMBLAST_TARBALL."
fi

mkdir -p "$REFINER_PREFIX"
cd "$REFINER_PREFIX"

# --- 1. RepeatModeler (supplies Refiner, MultAln.pm, Matrices/) -------------
if [[ -d RepeatModeler/.git ]]; then
  log "RepeatModeler already present, skipping clone"
else
  log "cloning RepeatModeler"
  git clone --depth 1 https://github.com/Dfam-consortium/RepeatModeler.git RepeatModeler
fi

# --- 2. RepeatMasker (supplies the .pm modules Refiner 'use lib's) ----------
if [[ -d RepeatMasker/.git ]]; then
  log "RepeatMasker already present, skipping clone"
else
  log "cloning RepeatMasker (modules only; no configure, no Libraries build)"
  git clone --depth 1 https://github.com/Dfam-consortium/RepeatMasker.git RepeatMasker
fi
for m in SearchResult.pm SearchResultCollection.pm WUBlastSearchEngine.pm \
         NCBIBlastSearchEngine.pm SeqDBI.pm SimpleBatcher.pm FastaDB.pm \
         SearchEngineI.pm Matrix.pm; do
  [[ -f "RepeatMasker/$m" ]] || { log "FATAL: RepeatMasker/$m missing"; exit 1; }
done

# --- 3. rmblast (prebuilt; NOT compiled) -----------------------------------
if [[ -x "rmblast-${RMBLAST_VER}/bin/rmblastn" ]]; then
  log "rmblast already present, skipping download"
else
  log "downloading $RMBLAST_TARBALL (~177 MB)"
  curl -fL --retry 3 -o "$RMBLAST_TARBALL" "$RMBLAST_URL"
  log "extracting"
  tar xzf "$RMBLAST_TARBALL"
  rm -f "$RMBLAST_TARBALL"
fi
# curl-downloaded files are not quarantined, but be explicit for the case where
# the tarball was fetched by a browser instead.
xattr -dr com.apple.quarantine "rmblast-${RMBLAST_VER}" 2>/dev/null || true

RMBLAST_BIN="$REFINER_PREFIX/rmblast-${RMBLAST_VER}/bin"
"$RMBLAST_BIN/rmblastn" -version >/dev/null || { log "FATAL: rmblastn will not run"; exit 1; }
log "rmblastn: $("$RMBLAST_BIN/rmblastn" -version | head -1)"

# --- 4. environment file ---------------------------------------------------
cat > "$REFINER_PREFIX/refiner_env.sh" <<EOF
# source this before calling Refiner
export REFINER_PREFIX="$REFINER_PREFIX"
# PERL5LIB, not REPEATMASKER_DIR: Refiner's 'use lib' runs before
# RepModelConfig::resolveConfiguration() applies the env override.
export PERL5LIB="\$REFINER_PREFIX/RepeatMasker\${PERL5LIB:+:\$PERL5LIB}"
# RMBLAST_DIR *is* honoured — it is read after resolveConfiguration().
export RMBLAST_DIR="\$REFINER_PREFIX/rmblast-${RMBLAST_VER}/bin"
# Harmless/unused by Refiner but keeps RepModelConfig quiet if reused elsewhere.
export REPEATMASKER_DIR="\$REFINER_PREFIX/RepeatMasker"
export PATH="\$REFINER_PREFIX/bin:\$PATH"
EOF

mkdir -p "$REFINER_PREFIX/bin"
cat > "$REFINER_PREFIX/bin/Refiner" <<EOF
#!/usr/bin/env bash
# Wrapper: sets the environment Refiner needs, then execs it.
set -euo pipefail
source "$REFINER_PREFIX/refiner_env.sh"
exec /usr/bin/perl "$REFINER_PREFIX/RepeatModeler/Refiner" "\$@"
EOF
chmod +x "$REFINER_PREFIX/bin/Refiner"

# --- 5. smoke test ---------------------------------------------------------
log "smoke test"
SMOKE="$(mktemp -d)"
python3 - "$SMOKE/t.fa" <<'PY'
import random, sys
random.seed(7)
base = "".join(random.choice("ACGT") for _ in range(300))
with open(sys.argv[1], "w") as fh:
    for i in range(12):
        s = list(base)
        for _ in range(20):
            j = random.randrange(len(s)); s[j] = random.choice("ACGT")
        # plain header, NO description field — Refiner reserves it
        fh.write(">c%d\n%s\n" % (i + 1, "".join(s)))
PY
( cd "$SMOKE" && "$REFINER_PREFIX/bin/Refiner" t.fa >refiner.log 2>&1 ) || true
if [[ -s "$SMOKE/t.fa.refiner_cons" ]]; then
  log "OK — consensus of length $(tail -n +2 "$SMOKE/t.fa.refiner_cons" | tr -d '\n' | wc -c | tr -d ' ') bp built"
  log "OK — Stockholm: $(ls "$SMOKE"/t.fa.refiner.stk 2>/dev/null || echo MISSING)"
else
  log "FAIL — no consensus produced.  Log:"
  cat "$SMOKE/refiner.log" >&2
  exit 1
fi
rm -rf "$SMOKE"

cat >&2 <<EOF

[setup_refiner] done.  To use:

    source $REFINER_PREFIX/refiner_env.sh
    Refiner /path/to/copies.fa

  Input rules (enforced by Refiner, not by this script):
    * plain FASTA, headers must be a bare identifier with NO description field
    * one sequence per genomic copy, already reverse-complemented for '-' copies

EOF
