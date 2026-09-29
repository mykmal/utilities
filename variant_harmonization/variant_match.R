# =============================================================================
# variant_match.R
#
# match_variants() resolves each row of a `query` table (a PGS scoring file)
# against a `target` table (a reference panel or a GWAS summary statistics
# file) and reports, per query row, which target row it corresponds to and how
# the two rows' allele labels relate. It deliberately knows nothing about
# weights, betas, or scoring files: downstream wrappers (harmonize_pgs.R,
# annotate_pgs_with_gwas.R) add the domain-specific columns on top.
#
# -----------------------------------------------------------------------------
# Quality control flags
#
# `allele_flipped` answers, is the target row's effect allele the same physical
# allele that the query calls its other allele? That is, does the target's
# effect estimate/frequency need reversing to speak about the query's effect
# allele? Every downstream use needs exactly this: the harmonizer uses it to
# pick which of the target row's two alleles to emit, and the annotator uses it
# to sign the beta and to take 1 - frequency.
#
# `strand_flipped` answers, do the two files write this site on opposite
# strands? Meanwhile, `palindromic` answers, is this site palindromic?
#
# `ambiguous` denotes palindromic variants whose strand could not be confidently
# resolved using frequency information, and `freq_mismatch` denotes variants
# whose harmonized effect allele frequencies are very different.
#
# `pos_shift` is the SIGNED number of base pairs that an indel had to be shifted
# in order to find a match: the matched target position minus the query
# position, both in trimmed coordinates. It is therefore 0 on the "exact" and
# "trimmed" rungs, and non-zero only on "trimmed_window", where its sign says
# which side of the query the matched target row sits on. It is NA exactly where
# `usable` is FALSE. NOTE: this column is NUMERIC, not logical -- select shifted
# rows with `pos_shift != 0`, never with `df[df$pos_shift, ]`, which would index
# by row position instead of masking.
#
# `target_freq_aligned` is the target's effect allele frequency re-expressed for
# the QUERY's effect allele (that is, 1 - f wherever `allele_flipped`), and is
# the quantity `freq_residual` is computed from. It is exposed so that wrappers
# needing an aligned frequency reuse this one rather than recomputing the same
# flip against a second copy of the rule.
#
# `usable` is the overall flag meaning "orientation resolved". In practice,
# downstream applications will often want to filter the output on
# `usable & !ambiguous & (!freq_mismatch | is.na(freq_mismatch))`.
#
# -----------------------------------------------------------------------------
# The matching ladder
#
# Candidates are sought in order of decreasing evidence strength, and the rung
# that produced the match is reported in `match_method`:
#
#   "exact"          Allele strings agree verbatim (directly or, for
#                    equal-length alleles, as reverse complements) at the exact
#                    query position. Strongest.
#   "trimmed"        Length-changing variants only. The two rows' parsimonious
#                    representations agree, at the same trimmed position. This
#                    reconciles padding differences. Example: query 3000 CT/C
#                    is the same 1 bp deletion as target 3000 CTT/CT.
#   "trimmed_window" Same as "trimmed", but the trimmed positions differ by up
#                    to indel_pos_tol. Must be unique or the row is refused.
#                    Reports the offset in pos_shift.
#
# Why are the two trimmed rungs restricted to length-changing variants?
#
# Trimming discards the shared flanking bases, and for equal-length alleles
# those bases are what identify the strand. Worked example: query AG/AC at 500
# trims to 501 G/C. A target row reading 501 G/C (same strand) and one reading
# 500 CT/GT (opposite strand, trims to 501 C/G) both match that trimmed key,
# the first as same_letter and the second as a label swap. Yet the truth is
# allele_flipped = FALSE in both cases, because the second is a strand flip.
# Trimming would therefore invert the orientation on the second, silently.
# Untrimmed, that same target row reads CT/GT at 500, matches by reverse
# complement, and is correctly called a strand flip.
#
# No such hazard exists for length-changing variants: rc() preserves length, so
# rc_ea == oa is impossible and an indel is never palindromic. Moreover,
# reverse-complement matching is already refused for them (see class 3 below).
# Strand and palindromic status are therefore fixed a priori for every row the
# trimmed rungs can reach, so trimming there loses nothing. Equal-length rows
# are matched on untrimmed strings only, and their strand and palindromic
# decisions are made from those untrimmed strings.
#
# Cases that are not handled (would need more than string arithmetic):
#
# Two representations placing the same indel at different positions inside a
# repeat tract (a left- vs right-aligned deletion in a homopolymer) trim to
# different positions and so only match via "trimmed_window", within
# indel_pos_tol. Beyond that they are reported as allele_mismatch. Resolving
# them properly needs either left-alignment against a reference FASTA, or a
# pairwise equivalence test that reconstructs the local reference from the two
# rows' own allele strings. Thus, `allele_mismatch` on an indel should be read
# as "possibly a representation difference", not "absent".
# =============================================================================

.VM_REQ_COLS <- c("chromosome", "position", "effect_allele",
                  "other_allele", "effect_allele_frequency")

# Defensive numeric coercion. NOTE: as.numeric() applied directly to a FACTOR
# returns the integer level codes, not the labels -- as.numeric(factor("1000"))
# is 1, not 1000 -- so every numeric coercion routes through as.character()
# first. Coercions that introduce new NAs warn rather than passing silently.
.vm_as_num <- function(x, what) {
  if (is.numeric(x)) return(x)
  y <- suppressWarnings(as.numeric(as.character(x)))
  n_new <- sum(is.na(y) & !is.na(x))
  if (n_new > 0L)
    warning(sprintf("%s: %d value(s) could not be coerced to numeric and became NA.",
                    what, n_new), call. = FALSE)
  y
}

# Frequencies drive the palindromic strand inference, so a percent-coded column
# (30 instead of 0.30) would silently invert strand calls rather than fail.
# Reject out-of-range values outright. The error message names the observed
# range to aid in manual QC.
.vm_chk_freq <- function(x, what) {
  bad <- !is.na(x) & (x < 0 | x > 1)
  if (any(bad))
    stop(sprintf(paste0("%s must be a proportion in [0, 1]; found %d value(s) outside ",
                        "that range (observed %.4g to %.4g). Percentages must be divided ",
                        "by 100 first, and missing data sentinels (-9, -99, 999) must be ",
                        "recoded to NA first."),
                 what, sum(bad), min(x[bad]), max(x[bad])), call. = FALSE)
  x
}

# Scalar parameter validation. Without this, a NA tolerance surfaces as
# "argument is not interpretable as logical" from deep inside the row loop, and
# a nonsense threshold (say maf_ambig_thresh = 8 for "8%") silently marks every
# palindromic variant ambiguous.
#
# An unbounded-above parameter is passed hi = Inf, but is.finite() still rejects
# an actual Inf argument, so the range is phrased as ">= lo" in that case rather
# than as "[lo, Inf]" -- which would name a value the check does not accept.
.vm_chk_param <- function(x, what, lo, hi) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lo || x > hi)
    stop(sprintf("%s must be a single finite number %s; got %s.", what,
                 if (is.finite(hi)) sprintf("in [%g, %g]", lo, hi)
                 else sprintf(">= %g", lo),
                 paste(utils::capture.output(dput(x)), collapse = "")),
         call. = FALSE)
  x
}

# Reverse-complement (not just base-complement): revcomp("AG") is "CT", not
# "TC". NA must come back as NA_character_, because paste(NA, collapse = "")
# yields the literal string "NA", which would then compare equal to an allele
# coded "NA" and could leak into output columns.
.vm_revcomp <- function(x) {
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    paste(rev(strsplit(chartr("ACGT", "TGCA", s), "")[[1]]), collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

# Chromosome label normalization. Without stripping the "chr" prefix, a PGS
# Catalog file using "chr1" against a panel using "1" matches nothing and
# reports every row as `unmatched` with no other signal. 23/24/26 are PLINK's
# numeric codes for X/Y/MT. 25 is deliberately left alone: it is XY
# (pseudo-autosomal) in PLINK but MT in some other conventions, so mapping it
# either way risks silently merging distinct contigs.
# Normalization runs on the handful of DISTINCT labels and is mapped back by
# match(), so the string passes cost nothing on a panel with tens of millions
# of rows but only ~25 distinct chromosomes.
.vm_norm_chr <- function(x) {
  u  <- unique(x)
  nu <- toupper(trimws(as.character(u)))
  nu <- sub("^CHR", "", nu)
  nu <- ifelse(nu == "23", "X",  nu)
  nu <- ifelse(nu == "24", "Y",  nu)
  nu <- ifelse(nu == "26", "MT", nu)
  nu <- ifelse(nu == "M",  "MT", nu)
  nu[match(x, u)]
}

# Sort key for ALREADY-NORMALISED labels. Unrecognized contigs collapse to Inf
# and are separated by the label itself downstream.
.vm_chr_key <- function(nx) {
  nx  <- ifelse(nx == "X", "23", ifelse(nx == "Y", "24", ifelse(nx == "MT", "26", nx)))
  key <- suppressWarnings(as.numeric(nx))
  ifelse(is.na(key), Inf, key)
}

# Row order for output: chromosome, then position. The secondary key is the
# NORMALISED label, not the raw one, so "1"/"chr1" and "X"/"23"/"chrX" rows are
# genuinely interleaved in position order instead of being grouped into
# separate blocks by spelling. method = "radix" makes the character ordering
# locale-independent (C collation), so the same inputs give the same row order
# on every machine. Rows with no resolvable position sort last within their
# chromosome (order()'s default na.last = TRUE).
.vm_order_rows <- function(chr, pos) {
  nx <- .vm_norm_chr(chr)
  order(.vm_chr_key(nx), nx, pos, method = "radix")
}

# Parsimonious ("trimmed") representation, computed from the two allele strings
# ALONE -- no reference sequence needed. Right-trim shared trailing bases, then
# left-trim shared leading bases (advancing the position), each stopping while
# both alleles still have at least one base left. This is the trimming half of
# `bcftools norm` / `vt normalize`; the left-ALIGNMENT half, which shifts an
# indel to the leftmost equivalent placement inside a repeat tract, is the part
# that genuinely requires the reference FASTA and is NOT done here.
#
# Trimming reconciles differently-padded representations of one event:
#   1000 CT/C, 1000 CTT/CT, 1000 CTTT/CTT   all -> 1000 CT/C
#   700 GAAT/GT, 700 GAATT/GTT              both -> 700 GAA/G
# Insertions and deletions stay distinct (1000 C/CT does not collapse onto
# 1000 CT/C), the operation is idempotent, and it is symmetric in the two
# alleles -- so it does not matter which one is the reference allele, which we
# have no way of knowing here.
#
# Order matters: right before left. Trimming never changes the LENGTH DIFFERENCE
# between the alleles (it removes the same count from both), so whether a
# variant is length-changing is invariant under trimming.
#
# Vectorised over rows. `act` carries the indices of the rows STILL trimming and
# shrinks on every iteration, so the total cost is the sum of the per-row trim
# depths rather than (number of rows) x (depth of the deepest row). That
# distinction is what stops one long allele from taxing every other row:
# recomputing the comparison across the full vector each iteration made a single
# 5 kb allele among 100k length-changing rows cost 16 s instead of 0.03 s, and a
# 50 kb one cost 166 s. Panels that carry long indels or SVs (gnomAD, TOPMed,
# HRC) hit this routinely. Each loop terminates after (min allele length - 1)
# iterations, which is 0 for the common case of an indel padded with a single
# anchor base.
.vm_trim <- function(pos, a1, a2) {
  if (length(a1) == 0L) return(list(pos = pos, a1 = a1, a2 = a2))
  n1 <- nchar(a1); n2 <- nchar(a2)
  k   <- integer(length(a1))
  act <- seq_along(a1)
  repeat {
    e1 <- n1[act] - k[act]; e2 <- n2[act] - k[act]
    act <- act[e1 > 1L & e2 > 1L &
               substring(a1[act], e1, e1) == substring(a2[act], e2, e2)]
    if (length(act) == 0L) break
    k[act] <- k[act] + 1L
  }
  a1 <- substring(a1, 1L, n1 - k); a2 <- substring(a2, 1L, n2 - k)
  n1 <- n1 - k; n2 <- n2 - k
  j   <- integer(length(a1))
  act <- seq_along(a1)
  repeat {
    s   <- j[act] + 1L
    act <- act[(n1[act] - j[act]) > 1L & (n2[act] - j[act]) > 1L &
               substring(a1[act], s, s) == substring(a2[act], s, s)]
    if (length(act) == 0L) break
    j[act] <- j[act] + 1L
  }
  list(pos = pos + j, a1 = substring(a1, j + 1L), a2 = substring(a2, j + 1L))
}

# Alleles must be non-empty ACGT strings. Without this, "-", "I"/"D", "N",
# other IUPAC codes, "" and the literal string "NA" all pass as sequence, get
# revcomp'd as if they were bases, and are emitted into the harmonised allele
# columns where no downstream tool can use them. A string "NA" is the worst of
# these: it survives is.na() and round-trips through a TSV indistinguishably
# from missing.
.vm_bad_allele <- function(a) is.na(a) | !grepl("^[ACGT]+$", a)

# Output-column collision reporting for the wrappers. Both build their result by
# starting from the query table and assigning a fixed set of QC, raw-target and
# derived columns over it, with the target's pass-through columns prefixed. Two
# names can therefore be claimed twice, and in both cases the caller's data is
# the copy that loses:
#
#   - a QUERY column sharing a name with an output column is replaced;
#   - a TARGET pass-through column whose PREFIXED name collides with an output
#     column is dropped (it is written first, then overwritten in place).
#
# The precedence itself is deliberate -- the function's own output has to win, or
# the QC columns could not be trusted -- but it should not happen quietly. The
# realistic trigger is re-running a wrapper on output it produced earlier, which
# collides on every output column at once.
.vm_warn_collisions <- function(query_names, passthrough_names, own_names,
                                query_label, target_label) {
  q <- function(x) paste(sprintf("`%s`", x), collapse = ", ")
  displaced <- intersect(passthrough_names, own_names)
  clobbered <- intersect(query_names, c(passthrough_names, own_names))
  if (length(displaced) > 0L)
    warning(sprintf(paste0("%s: pass-through column(s) %s collide with this function's ",
                           "own output after prefixing and were dropped. Rename them in ",
                           "%s to keep their values."),
                    target_label, q(displaced), target_label), call. = FALSE)
  if (length(clobbered) > 0L)
    warning(sprintf(paste0("%s: column(s) %s share a name with this function's output ",
                           "and were replaced. Rename them in %s to keep their values."),
                    query_label, q(clobbered), query_label), call. = FALSE)
  invisible(NULL)
}


# =============================================================================
# match_variants()
#
# Returns a list:
#   $qc     data.frame, one row per query row (in query order), columns:
#             status, match_method, palindromic, strand_flipped, allele_flipped,
#             ambiguous, target_freq_aligned, freq_mismatch, freq_residual,
#             pos_shift, usable, n_target_matches, target_idx
#   $query  the query table with its shared columns coerced (character
#           chromosome, numeric position, uppercase alleles, numeric frequency)
#   $target the target table after dropping unusable rows, coerced the same
#           way. `target_idx` indexes THIS table, not the caller's original.
#
# status is one of:
#   invalid_input, unmatched, allele_mismatch, allele_mismatch_indel_revcomp,
#   indel_window_ambiguous, palindromic_unresolved, match, allele_swap
# =============================================================================
match_variants <- function(query, target,
                           indel_pos_tol       = 1,
                           maf_ambig_thresh    = 0.08,
                           freq_resid_thresh   = 0.15,
                           freq_tie_tol        = 1e-8,
                           unmatched_warn_frac = 0.5,
                           query_label         = "query",
                           target_label        = "target") {
  stopifnot(all(.VM_REQ_COLS %in% names(query)), all(.VM_REQ_COLS %in% names(target)))
  .vm_chk_param(indel_pos_tol,       "indel_pos_tol",       0, Inf)
  .vm_chk_param(maf_ambig_thresh,    "maf_ambig_thresh",    0, 0.5)
  .vm_chk_param(freq_resid_thresh,   "freq_resid_thresh",   0, 1)
  .vm_chk_param(freq_tie_tol,        "freq_tie_tol",        0, 1)
  .vm_chk_param(unmatched_warn_frac, "unmatched_warn_frac", 0, 1)

  ql <- query_label; tl <- target_label

  query$chromosome  <- as.character(query$chromosome)
  target$chromosome <- as.character(target$chromosome)
  query$position    <- .vm_as_num(query$position,  sprintf("%s$position", ql))
  target$position   <- .vm_as_num(target$position, sprintf("%s$position", tl))
  query$effect_allele    <- toupper(as.character(query$effect_allele))
  query$other_allele     <- toupper(as.character(query$other_allele))
  target$effect_allele   <- toupper(as.character(target$effect_allele))
  target$other_allele    <- toupper(as.character(target$other_allele))
  query$effect_allele_frequency <- .vm_chk_freq(
    .vm_as_num(query$effect_allele_frequency, sprintf("%s$effect_allele_frequency", ql)),
    sprintf("%s$effect_allele_frequency", ql))
  target$effect_allele_frequency <- .vm_chk_freq(
    .vm_as_num(target$effect_allele_frequency, sprintf("%s$effect_allele_frequency", tl)),
    sprintf("%s$effect_allele_frequency", tl))

  # Target rows that can never resolve a match are dropped upfront: missing
  # chromosome/position, and non-ACGT alleles (which would otherwise be
  # revcomp'd as sequence and could match a query row's equally malformed
  # allele, producing a "match" that names an allele no genotype file
  # contains). Degenerate rows with effect == other are dropped for the same
  # reason. Dropping also keeps NA out of the allele comparisons below.
  t_drop <- is.na(target$chromosome) | is.na(target$position) |
            .vm_bad_allele(target$effect_allele) | .vm_bad_allele(target$other_allele) |
            target$effect_allele == target$other_allele
  if (any(t_drop))
    warning(sprintf(paste0("%s: dropped %d of %d row(s) with missing coordinates, ",
                           "non-ACGT alleles, or effect_allele == other_allele."),
                    tl, sum(t_drop), length(t_drop)), call. = FALSE)
  target <- target[!t_drop, , drop = FALSE]

  # Index the target by chromosome AND position. Sorting positions within each
  # chromosome lets candidate lookup below use binary search (findInterval)
  # instead of a linear scan of the whole chromosome, which would cost
  # O(n_query x n_target_per_chr) -- the dominant cost at panel scale. Only row
  # indices and one numeric column are split, not the whole table, so no other
  # target column is re-copied per chromosome group. Matched rows are looked up
  # by absolute index into `target` afterward.
  t_chr_norm     <- .vm_norm_chr(target$chromosome)
  o              <- order(t_chr_norm, target$position, method = "radix")
  chr_norm_sort  <- t_chr_norm[o]
  tgt_idx_by_chr <- split(o, chr_norm_sort)                   # absolute, position-sorted
  tgt_pos_by_chr <- split(target$position[o], chr_norm_sort)  # ascending within group

  # Parsimonious representations of the target's LENGTH-CHANGING rows, plus a
  # second position index over just those rows keyed on the TRIMMED position
  # (trimming moves positions, so the untrimmed index cannot serve this). Only
  # length-changing rows are trimmed: the trimmed rungs never consult any other
  # row, and trimming equal-length alleles would discard the flanking bases that
  # carry their strand information (see the header).
  t_len_chg <- nchar(target$effect_allele) != nchar(target$other_allele)
  t_tr_pos  <- target$position
  t_tr_ea   <- target$effect_allele
  t_tr_oa   <- target$other_allele
  li <- which(t_len_chg)
  if (length(li) > 0L) {
    z <- .vm_trim(target$position[li], target$effect_allele[li],
                  target$other_allele[li])
    t_tr_pos[li] <- z$pos; t_tr_ea[li] <- z$a1; t_tr_oa[li] <- z$a2
  }
  o_tr <- li[order(t_chr_norm[li], t_tr_pos[li], method = "radix")]
  tr_idx_by_chr <- split(o_tr, t_chr_norm[o_tr])
  tr_pos_by_chr <- split(t_tr_pos[o_tr], t_chr_norm[o_tr])

  # Duplicate-target check. Two target rows describing the SAME variant (same
  # site, same unordered allele pair -- so A/G and G/A count as duplicates,
  # because they are two representations of one variant) are both
  # allele-compatible with a query row, and the tie between them breaks on
  # target ROW ORDER. That makes the pass-through columns (INFO score, rsid,
  # se, ...) depend on how the file happened to be sorted. Multi-allelic sites
  # split across rows (A/G plus A/C) have different allele pairs and are NOT
  # flagged.
  #
  # The check runs on the TRIMMED representation, so two rows that differ only
  # in padding (1000 CT/C and 1000 CTT/CT) are recognized as the duplicates they
  # are -- the trimmed rungs below would otherwise treat them as two competing
  # candidates for the same query row.
  #
  # The cheap numeric key (chromosome index x 1e9 + position) is a COLLISION-
  # TOLERANT prefilter whose only job is to narrow the expensive string paste to
  # co-located rows instead of running it on every panel row. Rows that really
  # are co-located always hash equal, so no duplicate can be missed; a hash
  # collision -- which happens once positions reach 1e9, where the chromosome
  # and position fields start to alias, as they do on non-human assemblies --
  # merely admits a few extra rows into the subset. The definitive key below is
  # therefore built from the chromosome and position THEMSELVES rather than from
  # the hash, so those extras cannot be reported as duplicates of one another.
  if (nrow(target) > 1L) {
    chr_i    <- match(t_chr_norm, unique(t_chr_norm))
    colocate <- chr_i * 1e9 + t_tr_pos
    at_site  <- duplicated(colocate) | duplicated(colocate, fromLast = TRUE)
    if (any(at_site)) {
      a1  <- pmin(t_tr_ea[at_site], t_tr_oa[at_site])
      a2  <- pmax(t_tr_ea[at_site], t_tr_oa[at_site])
      key <- paste(t_chr_norm[at_site], t_tr_pos[at_site], a1, a2, sep = ":")
      n_dup <- sum(duplicated(key))
      if (n_dup > 0L)
        warning(sprintf(paste0("%s: %d row(s) duplicate an earlier row's site and allele ",
                               "pair (%d variant(s) affected). Ties between duplicates ",
                               "break on row order, so pass-through columns may depend on ",
                               "input sort order -- de-duplicate %s first if that matters."),
                        tl, n_dup, length(unique(key[duplicated(key)])), tl), call. = FALSE)
    }
  }

  n              <- nrow(query)
  status         <- character(n)
  match_method   <- rep(NA_character_, n)
  palindromic    <- rep(NA, n)
  strand_flipped <- rep(NA, n)
  allele_flipped <- rep(NA, n)
  ambiguous      <- rep(NA, n)
  pos_shift      <- rep(NA_real_, n)      # signed bp, 0 unless the window rung fired
  target_idx     <- rep(NA_integer_, n)
  n_matches      <- rep(NA_integer_, n)   # allele-compatible candidates seen

  # Rows that cannot be safely inspected are flagged before the loop: missing
  # coordinates, non-ACGT alleles, and the degenerate effect == other case.
  # That last one matters more than it looks: for a self-complementary allele
  # such as "AT", effect == other makes revcomp(effect) == other true, so the
  # row would be routed into the palindromic branch and given a
  # frequency-inferred orientation for a variant that has only one allele.
  q_ea <- query$effect_allele; q_oa <- query$other_allele
  invalid <- is.na(query$chromosome) | is.na(query$position) |
             .vm_bad_allele(q_ea) | .vm_bad_allele(q_oa) |
             (!is.na(q_ea) & !is.na(q_oa) & q_ea == q_oa)
  q_chr_norm <- .vm_norm_chr(query$chromosome)

  # revcomp is a per-string strsplit/vapply, so it runs once over the DISTINCT
  # allele strings and is mapped back by match() -- a scoring file with a
  # million rows typically has only a handful of distinct allele strings.
  u_all  <- unique(c(q_ea, q_oa))
  u_rc   <- .vm_revcomp(u_all)
  rc_ea_v <- u_rc[match(q_ea, u_all)]
  rc_oa_v <- u_rc[match(q_oa, u_all)]

  # Query-side parsimonious representations, again only for length-changing
  # rows. Invalid rows are excluded because nchar(NA) is 2 and substring() would
  # happily operate on the string "NA".
  q_len_chg <- !invalid & nchar(q_ea) != nchar(q_oa)
  q_tr_pos  <- query$position
  q_tr_ea   <- q_ea
  q_tr_oa   <- q_oa
  qi <- which(q_len_chg)
  if (length(qi) > 0L) {
    z <- .vm_trim(query$position[qi], q_ea[qi], q_oa[qi])
    q_tr_pos[qi] <- z$pos; q_tr_ea[qi] <- z$a1; q_tr_oa[qi] <- z$a2
  }
  # Reverse complements of the trimmed alleles, needed only for the class-3
  # diagnostic (a naively reverse-complemented indel) on the trimmed rungs.
  u_tr    <- unique(c(q_tr_ea, q_tr_oa))
  u_tr_rc <- .vm_revcomp(u_tr)
  rc_tr_ea_v <- u_tr_rc[match(q_tr_ea, u_tr)]
  rc_tr_oa_v <- u_tr_rc[match(q_tr_oa, u_tr)]

  # Resolve every candidate position range with ONE vectorised findInterval
  # call per chromosome, not one call per query row. findInterval validates its
  # lookup vector with is.unsorted()/anyNA() on each invocation, which is
  # O(n_target_per_chr) and would otherwise dominate the entire runtime --
  # worse than the linear scan it replaces. Ranges are half-open: the matching
  # rows for row i are idx_group[(lo + 1):hi], empty when hi == lo.
  exact_lo <- integer(n); exact_hi <- integer(n)
  tr_lo    <- integer(n); tr_hi    <- integer(n)   # trimmed position, exact
  trw_lo   <- integer(n); trw_hi   <- integer(n)   # trimmed position, +/- tol
  for (ch in intersect(unique(q_chr_norm), names(tgt_pos_by_chr))) {
    rows <- which(q_chr_norm == ch)
    pg   <- tgt_pos_by_chr[[ch]]
    q    <- query$position[rows]
    exact_lo[rows] <- findInterval(q - 0.5, pg)
    exact_hi[rows] <- findInterval(q + 0.5, pg)
  }
  # Same trick against the trimmed index. Bounds are computed for every row in
  # the group, but only consulted for length-changing rows, for which q_tr_pos
  # is the trimmed position (it equals the raw position elsewhere).
  for (ch in intersect(unique(q_chr_norm), names(tr_pos_by_chr))) {
    rows <- which(q_chr_norm == ch)
    pg   <- tr_pos_by_chr[[ch]]
    q    <- q_tr_pos[rows]
    tr_lo[rows] <- findInterval(q - 0.5, pg)
    tr_hi[rows] <- findInterval(q + 0.5, pg)
    if (indel_pos_tol > 0) {
      trw_lo[rows] <- findInterval(q - indel_pos_tol - 0.5, pg)
      trw_hi[rows] <- findInterval(q + indel_pos_tol + 0.5, pg)
    }
  }

  t_ea <- target$effect_allele; t_oa <- target$other_allele
  t_pos <- target$position

  # Classify candidates. Both alleles must correspond as a PAIR -- checking the
  # effect allele alone would let a different variant at the same site (C/T vs
  # G/T at a tri-allelic position) pass as a false match.
  #   1 = direct (same strand)
  #   2 = reverse complement, equal-length alleles only -- an acceptable match
  #   3 = reverse complement, length-changing alleles -- NOT accepted as a
  #       match, recorded only to distinguish a plain allele_mismatch from an
  #       input mangled by a tool that naively reverse-complemented padded
  #       indel strings
  classify <- function(cand, ea, oa, rc_ea, rc_oa, len_changing, tea, toa) {
    rea <- tea[cand]; roa <- toa[cand]
    direct <- (ea == rea & oa == roa) | (ea == roa & oa == rea)
    rcomp  <- (rc_ea == rea & rc_oa == roa) | (rc_ea == roa & rc_oa == rea)
    cls <- rep(NA_integer_, length(cand))
    cls[rcomp]  <- if (len_changing) 3L else 2L
    cls[direct] <- 1L
    cls
  }

  for (i in seq_len(n)) {
    if (invalid[i]) { status[i] <- "invalid_input"; next }

    ea <- q_ea[i]; oa <- q_oa[i]
    rc_ea <- rc_ea_v[i]; rc_oa <- rc_oa_v[i]
    pos <- query$position[i]

    # Length-CHANGING (not merely multi-base) is the predicate that matters.
    # Alleles of unequal length carry a VCF padding/anchor base, and the anchor
    # is drawn from the opposite side of the event on the opposite strand, so
    # the two strands' representations are not reverse complements of each
    # other and don't even share a position. Equal-length alleles (SNPs and
    # MNVs such as AC/GT) carry no anchor, so revcomp is exactly the right
    # transformation for those and is still applied.
    len_changing <- nchar(ea) != nchar(oa)

    # Strand ambiguity. The test is whether revcomp maps the effect allele's
    # REPRESENTATION onto the OTHER allele's representation, because that is
    # precisely when the written letters stop identifying which physical
    # allele is meant:
    #   A/T   -> revcomp("A") == "T"  : ambiguous (the classic palindrome)
    #   AG/CT -> revcomp("AG") == "CT": ambiguous (the MNV analogue)
    #   AT/GC -> revcomp("AT") == "AT": NOT ambiguous. Each allele is its own
    #            reverse complement, so "AT" denotes the same physical allele
    #            on either strand and the letters remain unambiguous.
    #   AA/GG -> revcomp("AA") == "TT": NOT ambiguous, TT is not confusable
    #            with GG.
    # rc_ea == oa implies equal length, so a palindromic pair is never an indel.
    palin <- rc_ea == oa

    idx_group <- tgt_idx_by_chr[[q_chr_norm[i]]]
    if (is.null(idx_group)) { status[i] <- "unmatched"; n_matches[i] <- 0L; next }

    # RUNG 1 -- verbatim allele strings at the exact position.
    lo <- exact_lo[i]; hi <- exact_hi[i]
    cand <- if (hi > lo) idx_group[(lo + 1L):hi] else integer(0)
    cls  <- classify(cand, ea, oa, rc_ea, rc_oa, len_changing, t_ea, t_oa)
    ok   <- which(cls == 1L | cls == 2L)
    method       <- "exact"
    cmp_pos      <- t_pos                 # position vector used for tie-breaking
    q_cmp_pos    <- pos
    saw_any      <- length(cand) > 0L
    saw_rc_indel <- any(cls == 3L, na.rm = TRUE)

    # The trimmed rungs escalate only for length-changing variants -- see the
    # header for why trimming equal-length alleles is unsafe. Each rung fires
    # only when no allele-COMPATIBLE candidate was found so far, not merely when
    # no candidate exists: requiring an empty position made the fallback nearly
    # dead code on a dense panel, where some unrelated variant sits at almost
    # every position and blocked the retry.
    tr_group <- if (len_changing) tr_idx_by_chr[[q_chr_norm[i]]] else NULL

    # RUNG 2 -- parsimonious representations, same trimmed position. Rescues
    # pure padding differences (3000 CT/C against 3000 CTT/CT).
    if (length(ok) == 0L && !is.null(tr_group)) {
      lo <- tr_lo[i]; hi <- tr_hi[i]
      cand_t <- if (hi > lo) tr_group[(lo + 1L):hi] else integer(0)
      if (length(cand_t) > 0L) {
        cls_t <- classify(cand_t, q_tr_ea[i], q_tr_oa[i], rc_tr_ea_v[i], rc_tr_oa_v[i],
                          len_changing, t_tr_ea, t_tr_oa)
        ok_t  <- which(cls_t == 1L | cls_t == 2L)
        saw_any      <- TRUE
        saw_rc_indel <- saw_rc_indel || any(cls_t == 3L, na.rm = TRUE)
        if (length(ok_t) > 0L) {
          cand <- cand_t; cls <- cls_t; ok <- ok_t; method <- "trimmed"
          cmp_pos <- t_tr_pos; q_cmp_pos <- q_tr_pos[i]
        }
      }
    }

    # RUNG 3 -- parsimonious representations within indel_pos_tol of the trimmed
    # position. What it rescues is anchor-base vs. first-changed-base coordinate
    # conventions, 0- vs. 1-based off-by-ones, and two equally non-canonical
    # placements inside a tandem repeat. What it does NOT rescue is a left- vs.
    # right-aligned placement whose offset exceeds the tolerance. Widening
    # indel_pos_tol only extends reach into repeat tracts, which is also where
    # unrelated same-length indels cluster -- hence the uniqueness requirement
    # below.
    via_window <- FALSE
    if (length(ok) == 0L && !is.null(tr_group) && indel_pos_tol > 0) {
      lo <- trw_lo[i]; hi <- trw_hi[i]
      cand_w <- if (hi > lo) tr_group[(lo + 1L):hi] else integer(0)
      if (length(cand_w) > 0L) {
        cls_w <- classify(cand_w, q_tr_ea[i], q_tr_oa[i], rc_tr_ea_v[i], rc_tr_oa_v[i],
                          len_changing, t_tr_ea, t_tr_oa)
        ok_w  <- which(cls_w == 1L | cls_w == 2L)
        saw_any      <- TRUE
        saw_rc_indel <- saw_rc_indel || any(cls_w == 3L, na.rm = TRUE)
        # Any compatible row here is necessarily at a SHIFTED trimmed position:
        # rows at the exact trimmed position were already classified in rung 2
        # under the identical predicate and found incompatible.
        if (length(ok_w) > 0L) {
          cand <- cand_w; cls <- cls_w; ok <- ok_w
          method <- "trimmed_window"; via_window <- TRUE
          cmp_pos <- t_tr_pos; q_cmp_pos <- q_tr_pos[i]
        }
      }
    }

    n_matches[i] <- length(ok)
    if (length(ok) == 0L) {
      # 0 candidates inspected means nothing is recorded near that position;
      # candidates that all failed the allele test mean something is there but
      # it is a different variant (or a differently-represented one).
      status[i] <- if (!saw_any) "unmatched"
                   else if (saw_rc_indel) "allele_mismatch_indel_revcomp"
                   else "allele_mismatch"
      next
    }
    # A window match must be UNIQUE. Inside a repeat tract several distinct
    # indels can share a length and sequence, so "nearest wins" would quietly
    # pick a different variant; refuse rather than guess.
    if (via_window && length(ok) > 1L) { status[i] <- "indel_window_ambiguous"; next }
    # Nearest position wins; at equal position a same-strand (direct) match
    # beats a reverse-complement one, so a genuine exact-letter match is never
    # overridden by a coincidental complement match at the same site. Remaining
    # exact ties are true duplicates and break on row order -- warned about
    # above rather than silently arbitrated.
    k <- cand[ok[order(abs(cmp_pos[cand[ok]] - q_cmp_pos), cls[ok])[1L]]]

    if (method != "exact") {
      # Reached only for length-changing variants, where rc_ea == oa is
      # impossible and reverse-complement matching is refused. Strand and
      # palindromic status are therefore fixed a priori, and the only thing the
      # trimmed strings decide is the label order -- which is exactly the
      # information trimming preserves. The untrimmed strings remain the sole
      # basis for the strand and palindromic decisions of every row that can
      # actually be palindromic or strand-flipped (see the header).
      palindromic[i]    <- FALSE
      strand_flipped[i] <- FALSE
      allele_flipped[i] <- q_tr_ea[i] == t_tr_oa[k] && q_tr_oa[i] == t_tr_ea[k]
      ambiguous[i]      <- FALSE
      status[i]         <- if (allele_flipped[i]) "allele_swap" else "match"
      match_method[i]   <- method
      pos_shift[i]      <- cmp_pos[k] - q_cmp_pos
      target_idx[i]     <- k
      next
    }

    rea <- t_ea[k]; roa <- t_oa[k]
    same_letter <- ea == rea && oa == roa
    swap_letter <- ea == roa && oa == rea
    same_comp   <- rc_ea == rea && rc_oa == roa
    swap_comp   <- rc_ea == roa && rc_oa == rea
    direct_ok   <- same_letter || swap_letter
    palindromic[i] <- palin

    if (palin) {
      # For a palindromic pair, raw letters tell us WITHOUT ambiguity whether
      # the target reports the same effect/other order as the query
      # (same_letter) or the opposite order (swap_letter) -- that part is a
      # plain string comparison. What letters cannot tell us is whether there
      # is ALSO a hidden strand mismatch on top, since for a self-complementary
      # pair "same order" and "opposite order, but secretly opposite strand"
      # produce identical letters either way.
      #
      # Frequency resolves that hidden part. `f_t_same` re-expresses the
      # target's frequency as "whichever target allele is reported using the
      # SAME LETTER as the query's effect allele" -- same_comp is NOT a
      # stand-in for this (for a palindromic pair it coincides with
      # swap_letter, not same_letter), so the baseline must be same_letter.
      f_q      <- query$effect_allele_frequency[i]
      f_t_same <- if (same_letter) target$effect_allele_frequency[k]
                  else 1 - target$effect_allele_frequency[k]
      if (!is.na(f_q) && !is.na(f_t_same)) {
        d_same <- abs(f_q - f_t_same); d_flip <- abs(f_q - (1 - f_t_same))
        # Compared with a tolerance: an exact tie arises whenever either
        # frequency is exactly 0.5, and without the tolerance the winner is
        # decided by floating-point rounding noise (0.5 - 0.7 evaluates
        # marginally smaller in magnitude than 0.5 - 0.3). Ties resolve
        # deterministically to "no hidden flip", and are always ambiguous.
        hidden <- d_flip < d_same - freq_tie_tol
        strand_flipped[i] <- hidden
        # A swapped label and a hidden strand flip cancel, hence XOR.
        allele_flipped[i] <- xor(swap_letter, hidden)
        # Ambiguity is capped by how close EITHER frequency sits to 0.5;
        # mirrors TwoSampleMR's default MAF > 0.42 ambiguity cutoff
        # (maf_ambig_thresh = 0.08 <-> MAF > 0.5 - 0.08 = 0.42).
        ambiguous[i] <- min(abs(f_q - 0.5), abs(f_t_same - 0.5)) < maf_ambig_thresh
      } else {
        # Genuinely unresolvable. Leaving allele_flipped NA (rather than
        # defaulting to FALSE) is what forces the wrappers to emit NA instead
        # of scoring on an allele that is wrong about half the time.
        status[i] <- "palindromic_unresolved"
        match_method[i] <- method
        target_idx[i] <- k   # still reported, so the row can be inspected
        next
      }
    } else if (direct_ok) {
      strand_flipped[i] <- FALSE
      allele_flipped[i] <- swap_letter
      ambiguous[i]      <- FALSE
    } else {
      # Reached only via cls == 2, so same_comp || swap_comp holds here.
      strand_flipped[i] <- TRUE
      allele_flipped[i] <- swap_comp
      ambiguous[i]      <- FALSE
    }

    status[i]       <- if (allele_flipped[i]) "allele_swap" else "match"
    match_method[i] <- method
    pos_shift[i]    <- cmp_pos[k] - q_cmp_pos
    target_idx[i]   <- k
  }

  # ---------------------------------------------------------------------------
  # Frequency concordance, for EVERY resolved row rather than only palindromic
  # ones. Align the target's frequency onto the query's effect allele and take
  # the residual. A large residual means a wrong match, an effect allele
  # mislabelled upstream, or frequencies from different ancestries -- and none
  # of those are specific to palindromes.
  #
  # For palindromic rows this reproduces exactly the winning orientation's
  # distance, min(d_same, d_flip): the orientation was chosen to minimise it,
  # so aligning by the chosen allele_flipped recovers the same number.
  #
  # `aligned` is returned as target_freq_aligned rather than kept local, because
  # a wrapper that needs the target's frequency stated for the query's effect
  # allele needs precisely this vector. Recomputing it wrapper-side duplicates
  # both the flip and the NA rule, and lets the wrapper's copy drift out of step
  # with the residual reported next to it.
  # ---------------------------------------------------------------------------
  t_freq_raw <- target$effect_allele_frequency[target_idx]
  aligned    <- t_freq_raw
  fl         <- !is.na(allele_flipped) & allele_flipped
  aligned[fl] <- 1 - t_freq_raw[fl]
  aligned[is.na(allele_flipped)] <- NA_real_
  freq_residual <- abs(query$effect_allele_frequency - aligned)
  freq_mismatch <- freq_residual > freq_resid_thresh   # NA if either freq missing

  # `usable` means only "the orientation is resolved", i.e. allele_flipped is
  # known, so the weight/beta can be aligned at all. It deliberately does NOT
  # fold in the QC flags, because whether to keep an ambiguous palindrome or a
  # frequency outlier is the analyst's call, not this function's. The
  # conservative filter is `usable & !ambiguous & !freq_mismatch` (with
  # !freq_mismatch guarded for NA where a frequency was missing).
  usable <- status %in% c("match", "allele_swap") & !is.na(allele_flipped)

  # A very high unmatched rate almost always means the two inputs disagree
  # about chromosome labelling or genome build, not that the variants are
  # genuinely absent, so surface it rather than returning a quietly empty
  # result. allele_mismatch is reported alongside because a build mismatch
  # produces both: positions that hit nothing, and positions that hit an
  # unrelated variant.
  n_considered <- sum(status != "invalid_input")
  n_unmatched  <- sum(status == "unmatched")
  n_mismatch   <- sum(status %in% c("allele_mismatch", "allele_mismatch_indel_revcomp"))
  if (n_considered > 0L && n_unmatched / n_considered > unmatched_warn_frac) {
    lab <- function(x) paste(utils::head(sort(unique(x[!is.na(x)])), 5L), collapse = ", ")
    warning(sprintf(paste0("%d of %d resolvable %s rows (%.1f%%) are unmatched (a further ",
                           "%d have allele mismatches). Check that %s and %s use the same ",
                           "genome build and chromosome labels (%s: %s; %s: %s)."),
                    n_unmatched, n_considered, ql, 100 * n_unmatched / n_considered,
                    n_mismatch, ql, tl,
                    ql, lab(query$chromosome), tl, lab(target$chromosome)),
            call. = FALSE)
  }

  list(
    qc = data.frame(
      status              = status,
      match_method        = match_method,
      palindromic         = palindromic,
      strand_flipped      = strand_flipped,
      allele_flipped      = allele_flipped,
      ambiguous           = ambiguous,
      target_freq_aligned = aligned,
      freq_mismatch       = freq_mismatch,
      freq_residual       = freq_residual,
      pos_shift           = pos_shift,
      usable              = usable,
      n_target_matches    = n_matches,
      target_idx          = target_idx,
      stringsAsFactors    = FALSE
    ),
    query  = query,
    target = target
  )
}
