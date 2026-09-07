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
  # these aren't orientation-dependent so no flip is applicable. They are
  # assigned before the columns below so that a GWAS column whose name collides
  # after prefixing (e.g., a GWAS column named "position_matched") cannot
  # overwrite this function's own QC output.
  extra_cols <- setdiff(names(gwas), c(.VM_REQ_COLS, "beta"))
  for (col in extra_cols) out[[paste0("gwas_", col)]] <- gwas[[col]][gi]

  out$match_status    <- qc$status
  out$match_method    <- qc$match_method
  out$usable          <- qc$usable
  out$palindromic     <- qc$palindromic
  out$strand_flipped  <- qc$strand_flipped
  out$allele_flipped  <- qc$allele_flipped
  out$ambiguous       <- qc$ambiguous
  out$freq_mismatch   <- qc$freq_mismatch
  out$freq_residual   <- qc$freq_residual
  out$pos_shift       <- qc$pos_shift
  out$n_gwas_matches  <- qc$n_target_matches

  # Add the raw, unaligned GWAS values for matched variants. These are
  # included in the matching engine's required columns, so they are not covered
  # by the pass-through above, but can be useful for manual verification.
  out$gwas_position_matched        <- gwas$position[gi]
  out$gwas_effect_allele           <- gwas$effect_allele[gi]
  out$gwas_other_allele            <- gwas$other_allele[gi]
  out$gwas_beta                    <- gwas$beta[gi]
  out$gwas_effect_allele_frequency <- gwas$effect_allele_frequency[gi]

  # Align with the PGS effect allele. Betas and allele frequencies on rows with
  # no resolved orientation get set to NA instead of being passed through.
  # Columns are built by subset assignment rather than ifelse() so that they
  # stay numeric even when nothing resolved.
  raw_beta <- gwas$beta[gi]
  raw_eaf  <- gwas$effect_allele_frequency[gi]
  do_flip  <- !is.na(qc$allele_flipped) & qc$allele_flipped
  beta_aligned <- raw_beta
  eaf_aligned  <- raw_eaf
  beta_aligned[do_flip] <- -raw_beta[do_flip]
  eaf_aligned[do_flip]  <- 1 - raw_eaf[do_flip]
  beta_aligned[!qc$usable] <- NA_real_
  eaf_aligned[!qc$usable]  <- NA_real_
  out$gwas_beta_aligned                    <- beta_aligned
  out$gwas_effect_allele_frequency_aligned <- eaf_aligned

  # Sort by chromosome and position, and then reset row names
  out <- out[.vm_order_rows(out$chromosome, out$position), , drop = FALSE]
  rownames(out) <- NULL
  out
}
