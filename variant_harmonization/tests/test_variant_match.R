# =============================================================================
# test_variant_match.R
#
# Test suite for variant_match.R and its two wrappers. Base R only, no test
# framework required. Run from anywhere:
#
#   Rscript variant_harmonization/tests/test_variant_match.R
#
# Exits non-zero if any assertion fails, so it can be wired straight into CI.
#
# Where a test encodes a decision rather than an obvious truth, the reason is
# stated inline -- several of these look wrong until you work through the
# genetics (A/G against C/T really is one variant on two strands, not two
# variants), and a future reader should not "fix" them.
# =============================================================================

# --- locate the scripts under test, relative to this file --------------------
.this_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) > 0L) dirname(normalizePath(f[1])) else getwd()
}
SRC <- normalizePath(file.path(.this_dir(), ".."))
source(file.path(SRC, "variant_match.R"))
source(file.path(SRC, "harmonize_pgs.R"))
source(file.path(SRC, "annotate_pgs_with_gwas.R"))

# --- harness -----------------------------------------------------------------
.pass <- 0L; .fail <- 0L; .failed <- character(0); .section <- ""
section <- function(x) { .section <<- x; cat("\n", x, "\n", sep = "") }
ok <- function(label, cond) {
  if (isTRUE(cond)) {
    .pass <<- .pass + 1L
  } else {
    .fail <<- .fail + 1L
    .failed <<- c(.failed, sprintf("[%s] %s", .section, label))
    cat("  FAIL: ", label, "\n", sep = "")
  }
}
eq <- function(label, a, b) ok(label, isTRUE(all.equal(a, b)))
errs <- function(expr) tryCatch({ expr; NA_character_ },
                                error = function(e) conditionMessage(e))
warns <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr,
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  w
}
quiet <- function(expr) suppressWarnings(expr)

# --- builders ----------------------------------------------------------------
qdf <- function(...) data.frame(..., stringsAsFactors = FALSE)
mkv <- function(chr, pos, ea, oa, freq)
  qdf(chromosome = chr, position = pos, effect_allele = ea,
      other_allele = oa, effect_allele_frequency = freq)
EMPTY <- mkv(character(0), numeric(0), character(0), character(0), numeric(0))
# a single query row matched against a single target row, returning just the qc
one <- function(q, t, ...) quiet(match_variants(q, t, ...))$qc


# =============================================================================
section("string and coordinate helpers")
# =============================================================================
eq("revcomp is reverse, not just complement", .vm_revcomp("AG"), "CT")
eq("revcomp of palindrome", .vm_revcomp("AT"), "AT")
eq("revcomp vectorised", .vm_revcomp(c("A", "CG", NA)), c("T", "CG", NA))
ok("revcomp keeps NA as NA_character_", is.na(.vm_revcomp(NA_character_)))
# paste(NA, collapse="") would yield the literal "NA", which round-trips through
# a TSV indistinguishably from missing -- hence the explicit NA branch.
ok("revcomp never produces the string 'NA'", !identical(.vm_revcomp(NA_character_), "NA"))

eq("chr prefix stripped", .vm_norm_chr("chr1"), "1")
eq("chr prefix case-insensitive", .vm_norm_chr(c("CHR2", "Chr3")), c("2", "3"))
eq("PLINK numeric codes", .vm_norm_chr(c("23", "24", "26")), c("X", "Y", "MT"))
eq("M normalises to MT", .vm_norm_chr("M"), "MT")
# 25 is XY (pseudo-autosomal) in PLINK but MT elsewhere, so mapping it either
# way risks silently merging distinct contigs. It must pass through untouched.
eq("25 deliberately left alone", .vm_norm_chr("25"), "25")
eq("whitespace trimmed", .vm_norm_chr(" 1 "), "1")
# NA in, NA out. Note the type is only guaranteed to be character when at least
# one label is non-NA: ifelse() returns a value shaped like its TEST, so an
# all-NA input yields logical NA. Harmless here -- target rows with a missing
# chromosome are dropped before normalisation, and query rows with one are
# already invalid_input -- so this asserts the NA, not the type.
ok("NA preserved", is.na(.vm_norm_chr(NA)))
eq("NA alongside a real label", .vm_norm_chr(c(NA, "chr1")), c(NA, "1"))
eq("unrecognised contig preserved", .vm_norm_chr("GL000191.1"), "GL000191.1")
ok("unrecognised contig sorts last", is.infinite(.vm_chr_key("GL000191.1")))

# Sorting on the NORMALISED label, so "1"/"chr1" interleave in position order
# rather than forming separate blocks by spelling.
o <- .vm_order_rows(c("chr1", "X", "1", "23", "MT", "GL000191.1", "2"),
                    c(5, 1, 1, 2, 1, 1, 1))
eq("mixed labels interleave by position", o[1:2], c(3L, 1L))
eq("unrecognised contig ordered last", o[7], 6L)
ok("rows with no position sort last within chromosome",
   .vm_order_rows(c("1", "1"), c(NA, 5))[1] == 2L)

ok("bad allele: dash", .vm_bad_allele("-"))
ok("bad allele: empty", .vm_bad_allele(""))
ok("bad allele: N", .vm_bad_allele("N"))
ok("bad allele: literal NA string", .vm_bad_allele("NA"))
ok("bad allele: lowercase (callers upcase first)", .vm_bad_allele("a"))
ok("good allele: single base", !.vm_bad_allele("A"))
ok("good allele: multi base", !.vm_bad_allele("ACGT"))


# =============================================================================
section(".vm_trim parsimonious representation")
# =============================================================================
z <- .vm_trim(c(1000, 1000, 1000, 700, 700),
              c("CT", "CTT", "CTTT", "GAAT", "GAATT"),
              c("C",  "CT",  "CTT",  "GT",   "GTT"))
eq("differently-padded representations collapse", z$a1, c("CT", "CT", "CT", "GAA", "GAA"))
eq("collapsed second alleles", z$a2, c("C", "C", "C", "G", "G"))
eq("positions unmoved when only right-trimming", z$pos, c(1000, 1000, 1000, 700, 700))
eq("idempotent", .vm_trim(z$pos, z$a1, z$a2), z)
# Left-trimming advances the position; right-trimming does not.
z2 <- .vm_trim(2999, "ACT", "AC")
eq("left trim advances position", z2$pos, 3000)
eq("left trim result", c(z2$a1, z2$a2), c("CT", "C"))
# Insertions and deletions must stay distinct: 1000 C/CT is not 1000 CT/C.
z3 <- .vm_trim(1000, "C", "CT")
eq("insertion not collapsed onto deletion", c(z3$a1, z3$a2), c("C", "CT"))
# Symmetric in the two alleles, since we cannot know which is the reference.
za <- .vm_trim(1000, "CTT", "CT"); zb <- .vm_trim(1000, "CT", "CTT")
eq("symmetric in allele order", c(za$a1, za$a2), c(zb$a2, zb$a1))
# Trimming removes the same count from both alleles, so length-changing status
# is invariant -- the property the two trimmed rungs depend on.
set.seed(101); b <- c("A", "C", "G", "T")
r1 <- vapply(1:400, function(i) {
  a1 <- paste(sample(b, sample(1:6, 1), TRUE), collapse = "")
  a2 <- paste(sample(b, sample(1:6, 1), TRUE), collapse = "")
  p  <- paste(sample(b, sample(0:3, 1), TRUE), collapse = "")
  s  <- paste(sample(b, sample(0:3, 1), TRUE), collapse = "")
  t  <- .vm_trim(1000, paste0(p, a1, s), paste0(p, a2, s))
  (nchar(t$a1) != nchar(t$a2)) == (nchar(a1) != nchar(a2))
}, logical(1))
ok("length-changing status invariant under trimming", all(r1))
eq("empty input tolerated", length(.vm_trim(numeric(0), character(0), character(0))$a1), 0L)


# =============================================================================
section("matching ladder: exact rung, SNPs")
# =============================================================================
r <- one(mkv("1", 100, "A", "G", 0.3), mkv("1", 100, "A", "G", 0.3))
ok("same-strand SNP matches", r$status == "match")
eq("method is exact", r$match_method, "exact")
ok("not palindromic", !r$palindromic)
ok("no strand flip", !r$strand_flipped)
ok("no allele flip", !r$allele_flipped)
ok("usable", r$usable)
eq("frequency residual zero", r$freq_residual, 0)
eq("one candidate seen", r$n_target_matches, 1L)

r <- one(mkv("1", 100, "A", "G", 0.3), mkv("1", 100, "G", "A", 0.7))
ok("label swap detected", r$status == "allele_swap" && r$allele_flipped)
ok("label swap is not a strand flip", !r$strand_flipped)
eq("swap residual zero once aligned", r$freq_residual, 0)

r <- one(mkv("1", 100, "A", "G", 0.3), mkv("1", 100, "T", "C", 0.3))
ok("strand flip detected", r$status == "match" && r$strand_flipped)
ok("strand flip alone is not an allele flip", !r$allele_flipped)

# A/G against C/T is ONE variant written on two strands with the labels also
# swapped -- rc(A)=T=other, rc(G)=C=effect -- not two different variants.
r <- one(mkv("1", 100, "A", "G", 0.3), mkv("1", 100, "C", "T", 0.7))
ok("strand flip plus swap", r$status == "allele_swap" &&
                            r$strand_flipped && r$allele_flipped)
eq("strand flip plus swap aligns", r$freq_residual, 0)

# Both alleles must correspond as a PAIR; a different variant at the same site
# must not pass on the effect allele alone.
r <- one(mkv("1", 100, "A", "G", 0.3), mkv("1", 100, "C", "G", 0.3))
ok("different variant at same site is a mismatch", r$status == "allele_mismatch")
eq("mismatch reports zero compatible candidates", r$n_target_matches, 0L)

r <- one(mkv("1", 100, "A", "G", 0.3), mkv("1", 500, "A", "G", 0.3))
ok("nothing at the position is unmatched", r$status == "unmatched")
r <- one(mkv("1", 100, "A", "G", 0.3), mkv("2", 100, "A", "G", 0.3))
ok("other chromosome is unmatched", r$status == "unmatched")
# A chromosome absent from the target must not error on the group lookup.
ok("query chromosome missing from target",
   one(mkv("9", 100, "A", "G", 0.3), mkv("1", 100, "A", "G", 0.3))$status == "unmatched")

# MNVs carry no padding base, so revcomp is exactly the right transformation.
r <- one(mkv("1", 500, "AG", "AC", 0.3), mkv("1", 500, "CT", "GT", 0.3))
ok("MNV strand flip via untrimmed strings",
   r$status == "match" && r$strand_flipped && !r$allele_flipped)
# ...and trimming an equal-length pair must NOT be attempted: query AG/AC trims
# to 501 G/C, which would match this row and invert the orientation silently.
ok("equal-length alleles never reach the trimmed rungs",
   !(one(mkv("1", 500, "AG", "AC", 0.3), mkv("1", 501, "G", "C", 0.3))$status
     %in% c("match", "allele_swap")))

r <- one(mkv("1", 100, "a", "g", 0.3), mkv("1", 100, "A", "G", 0.3))
ok("lowercase alleles upcased", r$status == "match")
ok("chr prefix mismatch tolerated",
   one(mkv("chr1", 100, "A", "G", 0.3), mkv("1", 100, "A", "G", 0.3))$status == "match")
ok("X against 23 tolerated",
   one(mkv("X", 100, "A", "G", 0.3), mkv("23", 100, "A", "G", 0.3))$status == "match")
ok("chr1 telomeric position",
   one(mkv("1", 248956422, "A", "G", 0.3), mkv("1", 248956422, "A", "G", 0.3))$status == "match")

# Multi-allelic site split across rows: the right row must win, and the site is
# not a duplicate (different allele pairs).
t <- rbind(mkv("1", 100, "A", "C", 0.1), mkv("1", 100, "A", "G", 0.3),
           mkv("1", 100, "A", "T", 0.2))
r <- one(mkv("1", 100, "A", "G", 0.3), t)
eq("multi-allelic picks the matching row", r$target_idx, 2L)
eq("multi-allelic counts one compatible candidate", r$n_target_matches, 1L)


# =============================================================================
section("matching ladder: indels and the trimmed rungs")
# =============================================================================
r <- one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 3000, "CT", "C", 0.1))
ok("identical indel representation", r$status == "match" && r$match_method == "exact")
ok("an indel is never palindromic", !r$palindromic)

# 1000 C/CT and 1000 CT/C are the same event with the labels swapped.
r <- one(mkv("1", 3000, "C", "CT", 0.1), mkv("1", 3000, "CT", "C", 0.9))
ok("insertion against deletion is a label swap",
   r$status == "allele_swap" && r$match_method == "exact")

r <- one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 3000, "CTT", "CT", 0.1))
ok("padding difference rescued", r$status == "match")
eq("padding difference uses trimmed rung", r$match_method, "trimmed")
eq("trimmed rung reports no shift", r$pos_shift, 0)

r <- one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 3000, "CT", "CTT", 0.9))
ok("padding difference plus swap", r$status == "allele_swap" &&
                                   r$match_method == "trimmed")

# A tool that naively reverse-complemented a padded indel produces a row that
# is not a match, but is worth distinguishing from a plain absence.
r <- one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 3000, "AG", "G", 0.1))
ok("naively revcomp'd indel gets its own status",
   r$status == "allele_mismatch_indel_revcomp")
ok("revcomp'd indel is not usable", !r$usable)

# Padding differences on both sides, resolved through trimming.
ok("both-side padding difference",
   one(mkv("1", 3000, "CTT", "CT", 0.1), mkv("1", 2999, "ACT", "AC", 0.1))$status == "match")


# =============================================================================
section("matching ladder: trimmed_window rung")
# =============================================================================
# The genuine window case is a repeat tract: identical trimmed alleles, but the
# two files place the event at different offsets inside the tract.
r <- one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 3001, "CT", "C", 0.1))
ok("repeat-tract offset rescued", r$status == "match")
eq("window rung reported", r$match_method, "trimmed_window")
eq("window shift is +1", r$pos_shift, 1)

r <- one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 2999, "CT", "C", 0.1))
eq("window shift is signed", r$pos_shift, -1)

r <- one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 3001, "C", "CT", 0.9))
ok("window plus label swap", r$status == "allele_swap" &&
                             r$match_method == "trimmed_window")

ok("shift beyond tolerance refused",
   !(one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 3002, "CT", "C", 0.1),
         indel_pos_tol = 1)$status %in% c("match", "allele_swap")))
ok("shift within widened tolerance accepted",
   one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 3002, "CT", "C", 0.1),
       indel_pos_tol = 2)$status == "match")
ok("tolerance of zero disables the window rung",
   !(one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 3001, "CT", "C", 0.1),
         indel_pos_tol = 0)$status %in% c("match", "allele_swap")))
eq("wide tolerance reports the full offset",
   one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 2996, "CT", "C", 0.1),
       indel_pos_tol = 5)$pos_shift, -4)

# Inside a repeat tract several distinct indels share a length and sequence, so
# "nearest wins" would quietly pick a different variant. Refuse instead.
t <- rbind(mkv("1", 2999, "CT", "C", 0.1), mkv("1", 3001, "CT", "C", 0.1))
r <- one(mkv("1", 3000, "CT", "C", 0.1), t, indel_pos_tol = 1)
ok("ambiguous window refused", r$status == "indel_window_ambiguous")
eq("ambiguous window counts candidates", r$n_target_matches, 2L)
ok("ambiguous window not usable", !r$usable)

# Stronger rungs must win even when a weaker candidate is also available, and
# the ladder must not be blocked by an unrelated variant at the position.
t <- rbind(mkv("1", 3001, "CT", "C", 0.1), mkv("1", 3000, "CT", "C", 0.4))
r <- one(mkv("1", 3000, "CT", "C", 0.4), t, indel_pos_tol = 1)
ok("exact beats window", r$match_method == "exact" && r$target_idx == 2L)
t <- rbind(mkv("1", 3001, "CT", "C", 0.1), mkv("1", 3000, "CTT", "CT", 0.4))
r <- one(mkv("1", 3000, "CT", "C", 0.4), t, indel_pos_tol = 1)
ok("trimmed beats window", r$match_method == "trimmed" && r$target_idx == 2L)
# An unrelated SNP sitting at the query position must not prevent escalation.
t <- rbind(mkv("1", 3000, "A", "G", 0.2), mkv("1", 3000, "CTT", "CT", 0.1))
r <- one(mkv("1", 3000, "CT", "C", 0.1), t)
ok("unrelated row at the position does not block escalation",
   r$match_method == "trimmed" && r$target_idx == 2L)

# A SNP off by one is NOT rescued: the window rung is indel-only by design.
ok("SNP off-by-one stays unmatched",
   one(mkv("1", 100, "A", "G", 0.3), mkv("1", 99, "A", "G", 0.3),
       indel_pos_tol = 5)$status == "unmatched")


# =============================================================================
section("palindromic variants and frequency inference")
# =============================================================================
# rc(effect) == other is exactly when the written letters stop identifying which
# physical allele is meant.
ok("A/T is palindromic", one(mkv("1", 100, "A", "T", 0.2),
                             mkv("1", 100, "A", "T", 0.2))$palindromic)
ok("AG/CT is the MNV analogue", one(mkv("1", 100, "AG", "CT", 0.2),
                                    mkv("1", 100, "AG", "CT", 0.2))$palindromic)
# AT/GC is NOT ambiguous: each allele is its own reverse complement, so the
# letters mean the same physical allele on either strand.
r <- one(mkv("1", 100, "AT", "GC", 0.2), mkv("1", 100, "AT", "GC", 0.2))
ok("AT/GC is not palindromic", !r$palindromic && r$status == "match")
# AA/GG is not confusable either: rc(AA)=TT, and TT is not GG.
ok("AA/GG is not palindromic", !one(mkv("1", 100, "AA", "GG", 0.2),
                                    mkv("1", 100, "AA", "GG", 0.2))$palindromic)

r <- one(mkv("1", 100, "A", "T", 0.2), mkv("1", 100, "A", "T", 0.2))
ok("concordant frequency implies no hidden flip", !r$strand_flipped)
ok("concordant palindrome resolves", r$status == "match" && !r$allele_flipped)
ok("concordant palindrome not ambiguous", !r$ambiguous)

r <- one(mkv("1", 100, "A", "T", 0.2), mkv("1", 100, "A", "T", 0.8))
ok("discordant frequency implies a hidden strand flip", r$strand_flipped)
ok("hidden flip inverts the orientation", r$allele_flipped &&
                                          r$status == "allele_swap")
eq("hidden flip aligns the frequency", r$freq_residual, 0)

# A swapped label and a hidden strand flip cancel, hence XOR.
r <- one(mkv("1", 100, "A", "T", 0.2), mkv("1", 100, "T", "A", 0.2))
ok("swapped label with hidden flip cancels", !r$allele_flipped && r$strand_flipped)

r <- one(mkv("1", 100, "A", "T", 0.48), mkv("1", 100, "A", "T", 0.48))
ok("frequency near 0.5 is ambiguous", r$ambiguous)
# ambiguous rows are still `usable`: whether to keep them is the analyst's call.
ok("ambiguous palindrome remains usable", r$usable)
ok("threshold respected", !one(mkv("1", 100, "A", "T", 0.30),
                               mkv("1", 100, "A", "T", 0.30),
                               maf_ambig_thresh = 0.08)$ambiguous)

# An exact tie arises whenever either frequency is exactly 0.5; without the
# tolerance the winner would be decided by floating-point noise.
r <- one(mkv("1", 100, "A", "T", 0.5), mkv("1", 100, "A", "T", 0.5))
ok("exact tie resolves to no hidden flip", !r$strand_flipped)
ok("exact tie is always ambiguous", r$ambiguous)
r <- one(mkv("1", 100, "A", "T", 0.5), mkv("1", 100, "A", "T", 0.3))
ok("one frequency at 0.5 still ties", !r$strand_flipped && r$ambiguous)

# Without a frequency the orientation is genuinely unknowable. Leaving
# allele_flipped NA is what stops a wrapper scoring on a coin flip.
r <- one(mkv("1", 100, "A", "T", NA_real_), mkv("1", 100, "A", "T", 0.2))
ok("unresolvable palindrome flagged", r$status == "palindromic_unresolved")
ok("unresolvable palindrome not usable", !r$usable)
ok("unresolvable palindrome leaves orientation NA", is.na(r$allele_flipped))
eq("unresolvable palindrome still names the row", r$target_idx, 1L)
ok("unresolvable palindrome with missing target frequency",
   one(mkv("1", 100, "A", "T", 0.2),
       mkv("1", 100, "A", "T", NA_real_))$status == "palindromic_unresolved")

# For a non-palindromic variant a missing frequency costs nothing: the letters
# already determine the orientation.
r <- one(mkv("1", 100, "A", "G", NA_real_), mkv("1", 100, "A", "G", 0.4))
ok("missing frequency harmless off-palindrome", r$usable)
ok("residual NA when a frequency is missing", is.na(r$freq_residual))
ok("freq_mismatch NA when a frequency is missing", is.na(r$freq_mismatch))


# =============================================================================
section("frequency concordance")
# =============================================================================
r <- one(mkv("1", 100, "A", "G", 0.10), mkv("1", 100, "A", "G", 0.40))
eq("residual computed off-palindrome", r$freq_residual, 0.30)
ok("residual over threshold flags mismatch", r$freq_mismatch)
# freq_mismatch is advisory: the row is still oriented, so it stays usable.
ok("frequency outlier remains usable", r$usable)
ok("residual under threshold does not flag",
   !one(mkv("1", 100, "A", "G", 0.30), mkv("1", 100, "A", "G", 0.35))$freq_mismatch)
eq("aligned target frequency exposed",
   one(mkv("1", 100, "A", "G", 0.3), mkv("1", 100, "G", "A", 0.7))$target_freq_aligned, 0.3)
# For palindromic rows the residual must reproduce the winning orientation's
# distance, min(d_same, d_flip), since the orientation was chosen to minimise it.
r <- one(mkv("1", 100, "A", "T", 0.20), mkv("1", 100, "A", "T", 0.75))
eq("palindromic residual is the winning distance", r$freq_residual, 0.05)

# The filter the header recommends must actually select what it claims.
q <- rbind(mkv("1", 100, "A", "G", 0.3), mkv("1", 200, "A", "T", 0.48),
           mkv("1", 300, "A", "G", 0.1), mkv("1", 400, "A", "T", NA_real_))
t <- rbind(mkv("1", 100, "A", "G", 0.3), mkv("1", 200, "A", "T", 0.48),
           mkv("1", 300, "A", "G", 0.6), mkv("1", 400, "A", "T", 0.2))
r <- one(q, t)
keep <- r$usable & !r$ambiguous & (!r$freq_mismatch | is.na(r$freq_mismatch))
eq("recommended filter keeps only the clean row", which(keep), 1L)


# =============================================================================
section("invalid rows and input validation")
# =============================================================================
q <- mkv(c("1", "1", "1", "1", "1"), c(100, NA, 100, 100, 100),
         c("A", "A", "-", "A", "NA"), c("G", "G", "G", "A", "G"), 0.3)
r <- one(q, mkv("1", 100, "A", "G", 0.3))
ok("valid row unaffected", r$status[1] == "match")
ok("missing position invalid", r$status[2] == "invalid_input")
ok("non-ACGT allele invalid", r$status[3] == "invalid_input")
# effect == other would make revcomp(effect) == other true for a
# self-complementary allele, routing a one-allele row into the palindrome branch.
ok("effect == other invalid", r$status[4] == "invalid_input")
# The string "NA" survives is.na() and round-trips through a TSV
# indistinguishably from missing, so it must be rejected as sequence.
ok("literal 'NA' allele invalid", r$status[5] == "invalid_input")
ok("invalid rows report no candidate count", all(is.na(r$n_target_matches[2:5])))
ok("invalid rows are not usable", !any(r$usable[2:5]))
ok("invalid rows have NA orientation", all(is.na(r$allele_flipped[2:5])))
ok("missing chromosome invalid",
   one(mkv(NA, 100, "A", "G", 0.3), mkv("1", 100, "A", "G", 0.3))$status == "invalid_input")
# Whitespace-padded alleles must not be silently accepted as sequence.
ok("padded allele rejected rather than matched",
   one(mkv("1", 100, " A", "G", 0.3), mkv("1", 100, "A", "G", 0.3))$status == "invalid_input")

# Frequencies drive the palindromic strand inference, so a percent-coded column
# would silently invert strand calls rather than fail.
e <- errs(match_variants(mkv("1", 100, "A", "G", 30), mkv("1", 100, "A", "G", 0.3)))
ok("percent-coded frequency rejected", grepl("proportion", e))
ok("range named in the message", grepl("30", e))
ok("missing-data sentinel rejected",
   grepl("proportion", errs(match_variants(mkv("1", 100, "A", "G", -9),
                                           mkv("1", 100, "A", "G", 0.3)))))
ok("target frequency also checked",
   grepl("proportion", errs(match_variants(mkv("1", 100, "A", "G", 0.3),
                                           mkv("1", 100, "A", "G", 99)))))
ok("all-NA frequency column accepted",
   !is.na(one(mkv("1", 100, "A", "G", NA_real_),
              mkv("1", 100, "A", "G", NA_real_))$status))

for (col in c("chromosome", "position", "effect_allele", "other_allele",
              "effect_allele_frequency")) {
  q <- mkv("1", 100, "A", "G", 0.3); q[[col]] <- NULL
  ok(sprintf("missing required column rejected: %s", col),
     !is.na(errs(match_variants(q, mkv("1", 100, "A", "G", 0.3)))))
}
ok("effect_weight required by harmonize_pgs",
   !is.na(errs(harmonize_pgs(mkv("1", 100, "A", "G", 0.3), mkv("1", 100, "A", "G", 0.3)))))
ok("beta required by the annotator",
   !is.na(errs(annotate_pgs_with_gwas(mkv("1", 100, "A", "G", 0.3),
                                      mkv("1", 100, "A", "G", 0.3)))))

# A NA tolerance would otherwise surface from deep inside the row loop, and a
# nonsense threshold would silently mark every palindrome ambiguous.
bad <- list(indel_pos_tol = NA, indel_pos_tol = -1, indel_pos_tol = Inf,
            maf_ambig_thresh = 8, maf_ambig_thresh = -0.1,
            freq_resid_thresh = -1, freq_resid_thresh = 2,
            freq_tie_tol = "x", unmatched_warn_frac = c(0.1, 0.2))
for (i in seq_along(bad)) {
  arg <- bad[i]
  e <- errs(do.call(match_variants,
                    c(list(mkv("1", 100, "A", "G", 0.3), mkv("1", 100, "A", "G", 0.3)), arg)))
  ok(sprintf("parameter rejected: %s = %s", names(arg),
             paste(format(arg[[1]]), collapse = ",")), !is.na(e))
}
# hi = Inf means "unbounded above", but is.finite() still rejects an actual Inf,
# so the message must not advertise a range that includes it.
ok("unbounded parameter message does not claim to accept Inf",
   grepl(">= 0", errs(match_variants(mkv("1", 100, "A", "G", 0.3),
                                     mkv("1", 100, "A", "G", 0.3), indel_pos_tol = Inf))))


# =============================================================================
section("warnings")
# =============================================================================
# Two rows describing the SAME variant are both allele-compatible, and the tie
# breaks on row order -- making pass-through columns depend on input sort order.
w <- warns(match_variants(mkv("1", 100, "A", "G", 0.3),
                          rbind(mkv("1", 100, "A", "G", 0.3), mkv("1", 100, "G", "A", 0.7))))
ok("duplicate site and allele pair warned", any(grepl("duplicate", w)))
# The check runs on the trimmed representation, so padding-only duplicates count.
ok("padding-only duplicate warned",
   any(grepl("duplicate", warns(match_variants(mkv("1", 1000, "CT", "C", 0.1),
       rbind(mkv("1", 1000, "CT", "C", 0.1), mkv("1", 1000, "CTT", "CT", 0.1)))))))
# A multi-allelic site split across rows has different allele pairs and is not
# a duplicate.
ok("multi-allelic split not warned",
   !any(grepl("duplicate", warns(match_variants(mkv("1", 100, "A", "G", 0.3),
       rbind(mkv("1", 100, "A", "G", 0.3), mkv("1", 100, "A", "C", 0.1)))))))

w <- warns(match_variants(mkv("1", 100, "A", "G", 0.3),
                          rbind(mkv("1", 100, "N", "G", 0.3), mkv("1", 100, "A", "G", 0.3))))
ok("dropped target rows warned", any(grepl("dropped", w)))
ok("degenerate target row warned",
   any(grepl("dropped", warns(match_variants(mkv("1", 100, "A", "G", 0.3),
       rbind(mkv("1", 100, "A", "A", 0.3), mkv("1", 100, "A", "G", 0.3)))))))

# A high unmatched rate almost always means a build or labelling disagreement,
# not genuinely absent variants.
w <- warns(match_variants(mkv("1", 100, "A", "G", 0.3), mkv("1", 900, "A", "G", 0.3)))
ok("high unmatched rate warned", any(grepl("unmatched", w)))
ok("warning names both label sets", any(grepl("chromosome labels", w)))
ok("unmatched warning suppressible",
   length(warns(match_variants(mkv("1", 100, "A", "G", 0.3), mkv("1", 900, "A", "G", 0.3),
                               unmatched_warn_frac = 1))) == 0L)
ok("coercion to numeric warns rather than passing silently",
   any(grepl("could not be coerced",
       warns(match_variants(qdf(chromosome = "1", position = "xyz", effect_allele = "A",
                                other_allele = "G", effect_allele_frequency = 0.3),
                            mkv("1", 100, "A", "G", 0.3))))))


# =============================================================================
section("degenerate inputs")
# =============================================================================
ok("empty query", nrow(quiet(match_variants(EMPTY, mkv("1", 100, "A", "G", 0.3)))$qc) == 0L)
ok("empty target", quiet(match_variants(mkv("1", 100, "A", "G", 0.3), EMPTY))$qc$status == "unmatched")
ok("both empty", nrow(quiet(match_variants(EMPTY, EMPTY))$qc) == 0L)
# Every target row dropped leaves an empty position index; the group lookup must
# return NULL rather than erroring.
ok("every target row dropped",
   quiet(match_variants(mkv("1", 100, "A", "G", 0.3),
                        mkv("1", 100, "-", "G", 0.3)))$qc$status == "unmatched")
ok("every query row invalid",
   quiet(match_variants(mkv(NA, NA, NA, NA, NA_real_),
                        mkv("1", 100, "A", "G", 0.3)))$qc$status == "invalid_input")
EMPTY_PGS <- cbind(EMPTY, effect_weight = numeric(0))
ok("empty harmonize", nrow(quiet(harmonize_pgs(EMPTY_PGS, mkv("1", 100, "A", "G", 0.3)))) == 0L)
ok("empty annotate", nrow(quiet(annotate_pgs_with_gwas(
     EMPTY_PGS, cbind(mkv("1", 100, "A", "G", 0.3), beta = 0.1)))) == 0L)
# target_idx indexes the RETURNED target, not the caller's original.
m <- quiet(match_variants(mkv("1", 100, "A", "G", 0.3),
                          rbind(mkv("1", 100, "N", "G", 0.3), mkv("1", 100, "A", "G", 0.3))))
eq("target_idx indexes the filtered target", nrow(m$target), 1L)
eq("target_idx resolves correctly", m$target$effect_allele[m$qc$target_idx], "A")


# =============================================================================
section("column types and coercion")
# =============================================================================
# as.numeric() on a FACTOR returns level codes, not labels, so every coercion
# must route through as.character() first.
q <- data.frame(chromosome = factor("1"), position = factor("1000"),
                effect_allele = factor("A"), other_allele = factor("G"),
                effect_allele_frequency = factor("0.3"))
r <- one(q, mkv("1", 1000, "A", "G", 0.3))
ok("factor columns survive coercion", r$status == "match")
eq("factor position not read as a level code", r$freq_residual, 0)
ok("integer positions", one(mkv("1", 100L, "A", "G", 0.3),
                            mkv("1", 100L, "A", "G", 0.3))$status == "match")
ok("character positions coerced",
   one(qdf(chromosome = "1", position = "100", effect_allele = "A",
           other_allele = "G", effect_allele_frequency = "0.3"),
       mkv("1", 100, "A", "G", 0.3))$status == "match")
ok("unparseable position becomes invalid rather than erroring",
   quiet(match_variants(qdf(chromosome = "1", position = "xyz", effect_allele = "A",
                            other_allele = "G", effect_allele_frequency = 0.3),
                        mkv("1", 100, "A", "G", 0.3)))$qc$status == "invalid_input")
if (requireNamespace("tibble", quietly = TRUE)) {
  tb <- tibble::as_tibble(qdf(chromosome = "1", position = 100, effect_allele = "A",
                              other_allele = "G", effect_allele_frequency = 0.3,
                              effect_weight = 0.5))
  tr <- tibble::as_tibble(qdf(chromosome = "1", position = 100, effect_allele = "A",
                              other_allele = "G", effect_allele_frequency = 0.3,
                              rsid = "rs1"))
  h <- harmonize_pgs(tb, tr)
  ok("tibble input to harmonize_pgs", nrow(h) == 1L && h$usable && h$ref_rsid == "rs1")
  a <- annotate_pgs_with_gwas(tb, tibble::as_tibble(
         qdf(chromosome = "1", position = 100, effect_allele = "A", other_allele = "G",
             effect_allele_frequency = 0.3, beta = 0.4)))
  ok("tibble input to the annotator", isTRUE(all.equal(a$gwas_beta_aligned, 0.4)))
} else {
  cat("  (tibble not installed; skipping 2 tests)\n")
}


# =============================================================================
section("harmonize_pgs")
# =============================================================================
pgs <- qdf(chromosome = c("1", "1", "2", "2"), position = c(100, 200, 300, 400),
           effect_allele = c("A", "A", "A", "A"), other_allele = c("G", "G", "T", "G"),
           effect_allele_frequency = c(0.3, 0.3, 0.2, 0.3),
           effect_weight = c(0.5, -0.5, 0.1, 0.2))
ref <- qdf(chromosome = c("1", "1", "2", "2"), position = c(100, 200, 300, 400),
           effect_allele = c("A", "C", "A", "G"), other_allele = c("G", "T", "T", "A"),
           effect_allele_frequency = c(0.3, 0.7, 0.2, 0.7),
           rsid = paste0("rs", 1:4), info = c(0.9, 0.8, 0.7, 0.6))
h <- harmonize_pgs(pgs, ref)
eq("all rows returned", nrow(h), 4L)
ok("all rows usable", all(h$usable))
ok("reference columns passed through", all(h$ref_rsid == paste0("rs", 1:4)))
ok("output sorted by chromosome then position", identical(h$position, c(100, 200, 300, 400)))
eq("row names reset", rownames(h), as.character(1:4))
# plink2 --score identifies the effect allele by matching letters, so naming one
# of the two true alleles resolves direction WITHOUT re-signing the weight.
ok("weights never re-signed", all(h$harmonized_weight == h$effect_weight))
i <- which(h$position == 200)
ok("strand-flipped row detected", h$strand_flipped[i] && h$allele_flipped[i])
eq("harmonized allele taken from the reference row", h$harmonized_effect_allele[i], "T")
eq("weight preserved through a strand flip", h$harmonized_weight[i], -0.5)
i <- which(h$position == 400)
eq("swapped row effect allele", h$harmonized_effect_allele[i], "A")
eq("swapped row other allele", h$harmonized_other_allele[i], "G")
ok("harmonized alleles are the reference row's own pair",
   all(pmin(h$harmonized_effect_allele, h$harmonized_other_allele) ==
       pmin(h$ref_effect_allele, h$ref_other_allele)))
# A harmonized frequency column would only restate a frequency the caller has.
ok("no harmonized_freq column", !"harmonized_freq" %in% names(h))
ok("raw reference frequency still available", "ref_effect_allele_frequency" %in% names(h))

# Rows with no resolved orientation must emit NA, not a coin-flip allele.
pgs1 <- qdf(chromosome = "1", position = 100, effect_allele = "A", other_allele = "T",
            effect_allele_frequency = NA_real_, effect_weight = 0.5)
h1 <- harmonize_pgs(pgs1, mkv("1", 100, "A", "T", 0.2))
ok("unresolved row emits NA allele", is.na(h1$harmonized_effect_allele))
ok("unresolved row emits NA weight", is.na(h1$harmonized_weight))
# Built by subset assignment rather than ifelse(), so the type survives even
# when nothing resolved.
ok("allele column stays character when nothing resolves",
   is.character(h1$harmonized_effect_allele))
ok("weight column stays numeric when nothing resolves", is.numeric(h1$harmonized_weight))


# =============================================================================
section("annotate_pgs_with_gwas")
# =============================================================================
gwas <- qdf(chromosome = c("1", "1"), position = c(100, 200),
            effect_allele = c("A", "C"), other_allele = c("G", "T"),
            effect_allele_frequency = c(0.3, 0.7), beta = c(0.8, 0.4),
            se = c(0.1, 0.2))
a <- annotate_pgs_with_gwas(pgs[1:2, ], gwas)
eq("beta passed through when aligned", a$gwas_beta_aligned[1], 0.8)
eq("beta re-signed when the labels swap", a$gwas_beta_aligned[2], -0.4)
eq("frequency reversed when the labels swap",
   a$gwas_effect_allele_frequency_aligned[2], 0.3)
ok("GWAS columns passed through", all(a$gwas_se == c(0.1, 0.2)))
eq("raw beta retained for auditing", a$gwas_beta, c(0.8, 0.4))
ok("beta not duplicated into the pass-through block",
   sum(names(a) == "gwas_beta") == 1L)
a1 <- annotate_pgs_with_gwas(pgs1, qdf(chromosome = "1", position = 100,
        effect_allele = "A", other_allele = "T", effect_allele_frequency = 0.2, beta = 0.9))
ok("unresolved row emits NA beta", is.na(a1$gwas_beta_aligned))
ok("unresolved row emits NA frequency", is.na(a1$gwas_effect_allele_frequency_aligned))
ok("beta column stays numeric when nothing resolves", is.numeric(a1$gwas_beta_aligned))


# =============================================================================
section("output column collisions")
# =============================================================================
base_pgs <- qdf(chromosome = "1", position = 100, effect_allele = "A",
                other_allele = "G", effect_allele_frequency = 0.3, effect_weight = 0.5)
base_ref <- qdf(chromosome = "1", position = 100, effect_allele = "A",
                other_allele = "G", effect_allele_frequency = 0.3, rsid = "rs1")
ok("clean inputs warn about nothing",
   length(warns(harmonize_pgs(base_pgs, base_ref))) == 0L)
ok("clean inputs warn about nothing in the annotator",
   length(warns(annotate_pgs_with_gwas(base_pgs, cbind(base_ref, beta = 0.4)))) == 0L)
# A reference column named "position_matched" is written to
# ref_position_matched, then overwritten there by the raw reference position --
# so the caller's values do not survive anywhere and must be reported.
w <- warns(harmonize_pgs(base_pgs, cbind(base_ref, position_matched = "KEEP")))
ok("displaced pass-through column warned", any(grepl("dropped", w)))
ok("displaced column named", any(grepl("ref_position_matched", w)))
h <- quiet(harmonize_pgs(base_pgs, cbind(base_ref, position_matched = "KEEP")))
ok("displaced values genuinely absent from the output",
   !any(vapply(h, function(c) any(as.character(c) == "KEEP", na.rm = TRUE), logical(1))))
# A query column sharing an output name is replaced.
p <- base_pgs; p$usable <- "S"; p$harmonized_weight <- "S"
w <- warns(harmonize_pgs(p, base_ref))
ok("clobbered query columns warned", any(grepl("replaced", w)))
ok("both clobbered columns named", any(grepl("usable", w) & grepl("harmonized_weight", w)))
ok("output wins over the query column", is.logical(quiet(harmonize_pgs(p, base_ref))$usable))
# ...including one colliding with a pass-through name rather than an output name.
p <- base_pgs; p$ref_rsid <- "S"
ok("query column colliding with a pass-through name warned",
   any(grepl("ref_rsid", warns(harmonize_pgs(p, base_ref)))))
eq("pass-through wins over the query column",
   quiet(harmonize_pgs(p, base_ref))$ref_rsid, "rs1")
# The realistic trigger: re-running a wrapper on output it produced earlier.
ok("re-harmonizing already-harmonized output is reported",
   length(unlist(regmatches(
     w <- warns(harmonize_pgs(harmonize_pgs(base_pgs, base_ref), base_ref)),
     gregexpr("`[^`]+`", w)))) >= 15L)
ok("annotator reports displaced pass-through columns too",
   any(grepl("gwas_position_matched",
       warns(annotate_pgs_with_gwas(base_pgs,
             cbind(base_ref, beta = 0.4, position_matched = "K"))))))
# Column order: pass-through block, then QC block, then raw target, then derived.
h <- harmonize_pgs(base_pgs, cbind(base_ref, info = 0.9))
ok("pass-through block precedes the QC block",
   which(names(h) == "ref_rsid") < which(names(h) == "match_status"))
ok("raw reference columns follow the QC block",
   which(names(h) == "ref_effect_allele") > which(names(h) == "n_ref_matches"))
eq("derived scoring column last", names(h)[length(names(h))], "harmonized_weight")


# =============================================================================
section("cross-cutting invariants")
# =============================================================================
# Assembled to hit every status at once.
q <- rbind(mkv("1", 100, "A", "G", 0.3),      # match
           mkv("1", 200, "A", "G", 0.3),      # allele_mismatch
           mkv("1", 300, "A", "G", 0.3),      # allele_mismatch_indel_revcomp
           mkv("1", 400, "A", "T", NA_real_), # palindromic_unresolved
           mkv("1", 500, "CT", "C", 0.1),     # indel_window_ambiguous
           mkv("1", 600, "-",  "G", 0.3),     # invalid_input
           mkv("1", 700, "CT", "C", 0.1),     # match via trimmed_window
           mkv("1", 800, "A", "G", 0.3))      # unmatched
t <- rbind(mkv("1", 100, "A", "G", 0.3), mkv("1", 200, "C", "G", 0.3),
           mkv("1", 300, "AG", "G", 0.1), mkv("1", 400, "A", "T", 0.2),
           mkv("1", 499, "CT", "C", 0.1), mkv("1", 501, "CT", "C", 0.1),
           mkv("1", 701, "CT", "C", 0.1), mkv("2", 800, "A", "G", 0.3))
r <- quiet(match_variants(q, t, indel_pos_tol = 1))$qc
ok("all documented statuses reachable",
   all(c("match", "allele_mismatch", "palindromic_unresolved",
         "indel_window_ambiguous", "invalid_input", "unmatched") %in% r$status))
# `usable` means exactly "orientation resolved".
ok("usable iff the orientation is known", identical(r$usable, !is.na(r$allele_flipped)))
# The annotator takes its aligned frequency from the matcher, which NAs on
# is.na(allele_flipped) rather than on !usable. The two rules must coincide.
ok("the matcher's NA rule matches the wrappers' usable rule",
   identical(is.na(r$allele_flipped), !r$usable))
ok("pos_shift is populated exactly where the row resolved",
   identical(!is.na(r$pos_shift), r$usable))
ok("match_method populated wherever a target row was named",
   all(!is.na(r$match_method[!is.na(r$target_idx)])))
ok("target_idx NA wherever nothing matched",
   all(is.na(r$target_idx[r$status %in% c("unmatched", "invalid_input")])))
ok("qc has one row per query row", nrow(r) == nrow(q))
# freq_residual is computed from target_freq_aligned, by construction.
ok("residual consistent with the aligned frequency",
   all(abs(r$freq_residual - abs(q$effect_allele_frequency - r$target_freq_aligned))
       < 1e-12, na.rm = TRUE))
ok("query returned in input order, coerced",
   identical(quiet(match_variants(q, t))$query$position, q$position))


# =============================================================================
section("randomised end-to-end invariants")
# =============================================================================
set.seed(11)
b <- c("A", "C", "G", "T"); rc1 <- c(A = "T", C = "G", G = "C", T = "A")
lets <- function(a, flip) ifelse(flip, rc1[a], a)

# --- SNPs: independent strand and label perturbation on each side ------------
N <- 4000
A1 <- sample(b, N, TRUE); A2 <- sample(b, N, TRUE)
drop <- A1 == A2 | rc1[A1] == A2        # degenerate, and palindromes (tested below)
A1 <- A1[!drop]; A2 <- A2[!drop]; N <- length(A1)
f   <- round(runif(N, 0.05, 0.45), 4)   # away from 0.5, so orientation is decidable
bt  <- round(rnorm(N), 4)
pos <- sample(1e7, N); chr <- as.character(sample(1:22, N, TRUE))
t_swap <- sample(c(TRUE, FALSE), N, TRUE); t_strand <- sample(c(TRUE, FALSE), N, TRUE)
q_swap <- sample(c(TRUE, FALSE), N, TRUE); q_strand <- sample(c(TRUE, FALSE), N, TRUE)
tgt <- qdf(chromosome = chr, position = pos,
           effect_allele = lets(ifelse(t_swap, A2, A1), t_strand),
           other_allele  = lets(ifelse(t_swap, A1, A2), t_strand),
           effect_allele_frequency = ifelse(t_swap, 1 - f, f),
           beta = ifelse(t_swap, -bt, bt))
pgsr <- qdf(chromosome = chr, position = pos,
            effect_allele = lets(ifelse(q_swap, A2, A1), q_strand),
            other_allele  = lets(ifelse(q_swap, A1, A2), q_strand),
            effect_allele_frequency = ifelse(q_swap, 1 - f, f),
            effect_weight = round(rnorm(N), 4))
key <- paste(chr, pos)
h <- harmonize_pgs(pgsr, tgt); h <- h[match(key, paste(h$chromosome, h$position)), ]
ok("randomised SNPs all resolve", all(h$usable))
ok("randomised SNPs all resolve on the exact rung", all(h$match_method == "exact"))
ok("strand flip recovered exactly", all(h$strand_flipped == (q_strand != t_strand)))
ok("allele flip recovered exactly", all(h$allele_flipped == (q_swap != t_swap)))
ok("harmonized effect allele is the query's physical allele in target letters",
   all(h$harmonized_effect_allele == lets(ifelse(q_swap, A2, A1), t_strand)))
ok("harmonized other allele likewise",
   all(h$harmonized_other_allele == lets(ifelse(q_swap, A1, A2), t_strand)))
ok("randomised weights never re-signed", all(h$harmonized_weight == h$effect_weight))
ok("randomised residuals vanish", max(h$freq_residual) < 1e-9)
ok("no spurious frequency mismatches", !any(h$freq_mismatch))
ok("no spurious ambiguity", !any(h$ambiguous))
ok("no spurious shifts", all(h$pos_shift == 0))
a <- annotate_pgs_with_gwas(pgsr, tgt); a <- a[match(key, paste(a$chromosome, a$position)), ]
ok("randomised betas signed for the query effect allele",
   max(abs(a$gwas_beta_aligned - ifelse(q_swap, -bt, bt))) < 1e-12)
ok("randomised frequencies aligned to the query effect allele",
   max(abs(a$gwas_effect_allele_frequency_aligned - a$effect_allele_frequency)) < 1e-12)

# --- palindromes: orientation carried only by the frequency ------------------
set.seed(12)
M   <- 2000
P1  <- sample(b, M, TRUE); P2 <- rc1[P1]
fp  <- round(runif(M, 0.05, 0.35), 4)   # far from 0.5, so the inference is decidable
btp <- round(rnorm(M), 4)
pp  <- sample(1e7, M); cp <- as.character(sample(1:22, M, TRUE))
ps  <- sample(c(TRUE, FALSE), M, TRUE)  # target label order
qs  <- sample(c(TRUE, FALSE), M, TRUE)  # query physical allele choice
hid <- sample(c(TRUE, FALSE), M, TRUE)  # hidden strand difference
tgp <- qdf(chromosome = cp, position = pp,
           effect_allele = ifelse(ps, P2, P1), other_allele = ifelse(ps, P1, P2),
           effect_allele_frequency = ifelse(ps, 1 - fp, fp),
           beta = ifelse(ps, -btp, btp))
pgp <- qdf(chromosome = cp, position = pp,
           effect_allele = lets(ifelse(qs, P2, P1), hid),
           other_allele  = lets(ifelse(qs, P1, P2), hid),
           effect_allele_frequency = ifelse(qs, 1 - fp, fp),
           effect_weight = round(rnorm(M), 4))
ap <- annotate_pgs_with_gwas(pgp, tgp)
ap <- ap[match(paste(cp, pp), paste(ap$chromosome, ap$position)), ]
ok("randomised palindromes all flagged", all(ap$palindromic))
ok("randomised palindromes all resolve", all(ap$usable))
ok("randomised palindromes unambiguous away from 0.5", !any(ap$ambiguous))
ok("randomised palindrome betas signed correctly",
   max(abs(ap$gwas_beta_aligned - ifelse(qs, -btp, btp))) < 1e-12)
ok("randomised palindrome residuals vanish", max(ap$freq_residual) < 1e-9)

# --- indels: random padding on each side, random label order ----------------
set.seed(13)
K    <- 3000
anch <- sample(b, K, TRUE)
ins  <- vapply(sample(1:5, K, TRUE), function(k) paste(sample(b, k, TRUE), collapse = ""), "")
padq <- vapply(sample(0:3, K, TRUE), function(k) paste(sample(b, k, TRUE), collapse = ""), "")
padt <- vapply(sample(0:3, K, TRUE), function(k) paste(sample(b, k, TRUE), collapse = ""), "")
fi   <- round(runif(K, 0.05, 0.45), 4); bti <- round(rnorm(K), 4)
pi_  <- sample(1e7, K); ci <- as.character(sample(1:22, K, TRUE))
tsw  <- sample(c(TRUE, FALSE), K, TRUE); qsw <- sample(c(TRUE, FALSE), K, TRUE)
tgi <- qdf(chromosome = ci, position = pi_,
           effect_allele = ifelse(tsw, paste0(anch, ins, padt), paste0(anch, padt)),
           other_allele  = ifelse(tsw, paste0(anch, padt), paste0(anch, ins, padt)),
           effect_allele_frequency = ifelse(tsw, 1 - fi, fi),
           beta = ifelse(tsw, -bti, bti))
pgi <- qdf(chromosome = ci, position = pi_,
           effect_allele = ifelse(qsw, paste0(anch, ins, padq), paste0(anch, padq)),
           other_allele  = ifelse(qsw, paste0(anch, padq), paste0(anch, ins, padq)),
           effect_allele_frequency = ifelse(qsw, 1 - fi, fi),
           effect_weight = round(rnorm(K), 4))
ki <- paste(ci, pi_)
ai <- annotate_pgs_with_gwas(pgi, tgi); ai <- ai[match(ki, paste(ai$chromosome, ai$position)), ]
ok("randomised indels all resolve", all(ai$usable))
ok("randomised indels exercise the trimmed rung", any(ai$match_method == "trimmed"))
ok("an indel is never palindromic", all(!ai$palindromic))
# rc() preserves length, so reverse-complement matching is refused for indels
# and strand is fixed a priori.
ok("indels are never strand flipped", all(!ai$strand_flipped))
ok("indel label order recovered exactly", all(ai$allele_flipped == (qsw != tsw)))
ok("indel betas signed correctly",
   max(abs(ai$gwas_beta_aligned - ifelse(qsw, -bti, bti))) < 1e-12)
ok("indel residuals vanish", max(ai$freq_residual) < 1e-9)
hi <- harmonize_pgs(pgi, tgi); hi <- hi[match(ki, paste(hi$chromosome, hi$position)), ]
ok("indel harmonized pair is the target row's own pair",
   all(pmin(hi$harmonized_effect_allele, hi$harmonized_other_allele) ==
       pmin(hi$ref_effect_allele, hi$ref_other_allele)))


# =============================================================================
section("regressions")
# =============================================================================
# The duplicate prefilter packs chromosome and position into one double. A fixed
# 1e9 multiplier makes the two fields alias once positions reach 1e9, so rows on
# DIFFERENT chromosomes hashed equal and were reported as duplicates of one
# another. The definitive key is built from the chromosome and position instead.
w <- warns(match_variants(mkv("A", 1e9 + 1, "A", "G", 0.3),
                          rbind(mkv("A", 1e9 + 1, "A", "G", 0.3),
                                mkv("B", 1,       "A", "G", 0.3))))
ok("no false duplicate when positions reach 1e9", !any(grepl("duplicate", w)))
# ...while a real duplicate at such a position is still caught.
ok("real duplicate still caught at 1e9",
   any(grepl("duplicate", warns(match_variants(mkv("A", 1e9 + 1, "A", "G", 0.3),
       rbind(mkv("A", 1e9 + 1, "A", "G", 0.3), mkv("A", 1e9 + 1, "G", "A", 0.7)))))))

# pos_shift is a signed base-pair count, not a flag. Filtering it as a logical
# would index by row position instead of masking.
r <- one(mkv("1", 3000, "CT", "C", 0.1), mkv("1", 3001, "CT", "C", 0.1))
ok("pos_shift is numeric", is.numeric(r$pos_shift))
eq("pos_shift reports the offset, not TRUE", r$pos_shift, 1)
ok("pos_shift is not a logical in disguise", !is.logical(r$pos_shift))

# .vm_trim advances only the rows still trimming. Recomputing the comparison
# across the full vector each iteration cost (rows x depth of the deepest row),
# so ONE long allele taxed every other row: 100k rows with a single 20 kb allele
# took ~66 s. The bound is deliberately loose; the point is the shape, not the
# constant.
long <- paste(sample(b, 20000, TRUE), collapse = "")
a1 <- c(paste0(long, "T"), rep("CT", 1e5 - 1L))
a2 <- c(long,              rep("C",  1e5 - 1L))
el <- system.time(z <- .vm_trim(rep(1000, 1e5), a1, a2))[["elapsed"]]
cat(sprintf("  (one 20 kb allele among 100k rows: %.3f s)\n", el))
ok("one long allele does not tax every other row", el < 5)
eq("ordinary rows still trimmed correctly beside a long one", z$a1[2], "CT")
# The long row is a 1 bp insertion buried under a 20 kb shared prefix, so it
# must trim down to the minimal pair: one allele a single base, the other that
# base plus one. (Which base is whatever the prefix ended with, not "T" -- the
# left trim stops when the shorter allele has one base left.)
eq("the long row trims to a minimal 1 bp indel",
   nchar(z$a1[1]) - nchar(z$a2[1]), 1L)
eq("the long row is fully trimmed", min(nchar(z$a1[1]), nchar(z$a2[1])), 1L)
eq("the long row's alleles still agree on their shared base",
   substring(z$a1[1], 1L, 1L), z$a2[1])

# harmonize_pgs never re-signs a weight, so a harmonized frequency column could
# only restate one the caller already has.
ok("harmonize_pgs emits no harmonized_freq",
   !"harmonized_freq" %in% names(harmonize_pgs(base_pgs, base_ref)))
# The annotator's aligned frequency is the matcher's, not a second copy.
m <- match_variants(base_pgs, cbind(base_ref, beta = 0.4))
a <- annotate_pgs_with_gwas(base_pgs, cbind(base_ref, beta = 0.4))
eq("annotator reuses the matcher's aligned frequency",
   a$gwas_effect_allele_frequency_aligned, m$qc$target_freq_aligned)


# =============================================================================
cat(sprintf("\n%s\n=== %d passed, %d failed ===\n", strrep("-", 60), .pass, .fail))
if (.fail > 0L) {
  cat(paste0("  ", .failed, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("OK\n")
