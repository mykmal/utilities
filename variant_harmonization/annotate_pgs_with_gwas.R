# =============================================================================
# annotate_pgs_with_gwas.R
#
# Attach GWAS summary statistics to a polygenic score (PGS) weight file,
# aligning GWAS effect sizes and frequencies to the PGS effect alleles. For most
# use cases, filter the output on
# `usable & !ambiguous & (!freq_mismatch | is.na(freq_mismatch))`.
#
# Requires variant_match.R to be sourced first. All variant matching lives in
# match_variants(); this wrapper only signs the beta and reverses the frequency.
#
# Output QC columns from match_variants():
#   match_status, match_method, usable, palindromic, strand_flipped,
#   allele_flipped, ambiguous, freq_mismatch, freq_residual, pos_shift,
#   n_gwas_matches
#
# pos_shift is NUMERIC (signed base pairs, 0 unless the window rung fired), so
# select shifted rows with `pos_shift != 0` rather than treating it as a flag.
#
# match_method records which rung of the matching ladder resolved the row:
#   "exact" (verbatim alleles match at the exact position),
#   "trimmed" (parsimonious indel representations agree at the exact position),
#   "trimmed_window" (parsimonious indel representations agree within a window)
#
# Raw GWAS values at matched rows, for auditing the decision:
#   gwas_position_matched, gwas_effect_allele, gwas_other_allele, gwas_beta,
#   gwas_effect_allele_frequency
#
# GWAS values aligned to the PGS effect allele:
#   gwas_beta_aligned, gwas_effect_allele_frequency_aligned
# =============================================================================
annotate_pgs_with_gwas <- function(pgs, gwas,
                                   indel_pos_tol       = 1,
                                   maf_ambig_thresh    = 0.08,
                                   freq_resid_thresh   = 0.15,
                                   freq_tie_tol        = 1e-8,
                                   unmatched_warn_frac = 0.5) {
  stopifnot("beta" %in% names(gwas))
  gwas$beta <- .vm_as_num(gwas$beta, "gwas$beta")

  m <- match_variants(pgs, gwas,
                      indel_pos_tol       = indel_pos_tol,
                      maf_ambig_thresh    = maf_ambig_thresh,
                      freq_resid_thresh   = freq_resid_thresh,
                      freq_tie_tol        = freq_tie_tol,
                      unmatched_warn_frac = unmatched_warn_frac,
                      query_label = "pgs", target_label = "gwas")
  qc   <- m$qc
  pgs  <- m$query    # coerced: character chromosome, numeric position, uppercase alleles
  gwas <- m$target   # filtered + coerced; qc$target_idx indexes THIS table
  gi   <- qc$target_idx

  out <- pgs

  # Pass through any other GWAS columns (se, pvalue, n, rsid, ...) verbatim;
  # these aren't orientation-dependent so no flip is applicable.
  extra_cols <- setdiff(names(gwas), c(.VM_REQ_COLS, "beta"))
  add <- list()
  for (col in extra_cols) add[[paste0("gwas_", col)]] <- gwas[[col]][gi]

  # Align the beta with the PGS effect allele. Betas on rows with no resolved
  # orientation get set to NA instead of being passed through. The column is
  # built by subset assignment rather than ifelse() so that it stays numeric
  # even when nothing resolved.
  raw_beta     <- gwas$beta[gi]
  do_flip      <- !is.na(qc$allele_flipped) & qc$allele_flipped
  beta_aligned <- raw_beta
  beta_aligned[do_flip]    <- -raw_beta[do_flip]
  beta_aligned[!qc$usable] <- NA_real_

  # This function's own output. Assembled as a list and merged below rather than
  # assigned one column at a time, so that the names it claims can be compared
  # against the names already present instead of having to be restated in a
  # hard-coded list that would drift as columns are added.
  own <- list(
    match_status   = qc$status,
    match_method   = qc$match_method,
    usable         = qc$usable,
    palindromic    = qc$palindromic,
    strand_flipped = qc$strand_flipped,
    allele_flipped = qc$allele_flipped,
    ambiguous      = qc$ambiguous,
    freq_mismatch  = qc$freq_mismatch,
    freq_residual  = qc$freq_residual,
    pos_shift      = qc$pos_shift,
    n_gwas_matches = qc$n_target_matches,
    # The raw, unaligned GWAS values for matched variants. These are included in
    # the matching engine's required columns, so they are not covered by the
    # pass-through above, but can be useful for manual verification.
    gwas_position_matched        = gwas$position[gi],
    gwas_effect_allele           = gwas$effect_allele[gi],
    gwas_other_allele            = gwas$other_allele[gi],
    gwas_beta                    = gwas$beta[gi],
    gwas_effect_allele_frequency = gwas$effect_allele_frequency[gi],
    gwas_beta_aligned            = beta_aligned,
    # The aligned frequency is exactly the flip match_variants() already applied
    # to compute freq_residual, so it is taken from there rather than recomputed.
    # That keeps one implementation of both the flip and the "NA where the
    # orientation is unresolved" rule, and guarantees this column and the
    # residual reported beside it can never disagree.
    gwas_effect_allele_frequency_aligned = qc$target_freq_aligned
  )

  # Merge, with this function's own output taking precedence over a colliding
  # pass-through column and over a colliding PGS column, and report both cases.
  # own[] assigns in place where the name already exists and appends otherwise,
  # so the column order is the pass-through block followed by the QC block.
  .vm_warn_collisions(names(out), names(add), names(own), "pgs", "gwas")
  add[names(own)] <- own
  for (nm in names(add)) out[[nm]] <- add[[nm]]

  # Sort by chromosome and position, and then reset row names
  out <- out[.vm_order_rows(out$chromosome, out$position), , drop = FALSE]
  rownames(out) <- NULL
  out
}
