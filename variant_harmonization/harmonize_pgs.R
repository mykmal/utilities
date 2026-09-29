# =============================================================================
# harmonize_pgs.R
#
# Align a polygenic score (PGS) weight file with a genotype dataset, accounting
# for strand flips but keeping the score's original effect directions. For most
# use cases, filter the output on
# `usable & !ambiguous & (!freq_mismatch | is.na(freq_mismatch))`.
#
# Requires variant_match.R to be sourced first. All variant matching lives in
# match_variants(); this wrapper only turns the results into a scoring file.
#
# Output QC columns from match_variants():
#   match_status, match_method, usable, palindromic, strand_flipped,
#   allele_flipped, ambiguous, freq_mismatch, freq_residual, pos_shift,
#   n_ref_matches
#
# pos_shift is NUMERIC (signed base pairs, 0 unless the window rung fired), so
# select shifted rows with `pos_shift != 0` rather than treating it as a flag.
#
# match_method records which rung of the matching ladder resolved the row:
#   "exact" (verbatim alleles match at the exact position),
#   "trimmed" (parsimonious indel representations agree at the exact position),
#   "trimmed_window" (parsimonious indel representations agree within a window)
#
# Raw genotype dataset values at matched rows, for auditing the decision:
#   ref_position_matched, ref_effect_allele, ref_other_allele,
#   ref_effect_allele_frequency
#
# Scoring columns:
#   harmonized_effect_allele, harmonized_other_allele, harmonized_weight
#
# No harmonized frequency column is emitted. This wrapper never re-signs a
# weight, so such a column could only restate a frequency the caller already
# has: the PGS's own effect_allele_frequency (unchanged, since it already refers
# to the PGS effect allele) or ref_effect_allele_frequency flipped wherever
# allele_flipped. Callers wanting the latter can take it from match_variants()'s
# target_freq_aligned directly.
# =============================================================================
harmonize_pgs <- function(pgs, ref,
                          indel_pos_tol       = 1,
                          maf_ambig_thresh    = 0.08,
                          freq_resid_thresh   = 0.15,
                          freq_tie_tol        = 1e-8,
                          unmatched_warn_frac = 0.5) {
  stopifnot("effect_weight" %in% names(pgs))
  pgs$effect_weight <- .vm_as_num(pgs$effect_weight, "pgs$effect_weight")

  m <- match_variants(pgs, ref,
                      indel_pos_tol       = indel_pos_tol,
                      maf_ambig_thresh    = maf_ambig_thresh,
                      freq_resid_thresh   = freq_resid_thresh,
                      freq_tie_tol        = freq_tie_tol,
                      unmatched_warn_frac = unmatched_warn_frac,
                      query_label = "pgs", target_label = "ref")
  qc  <- m$qc
  pgs <- m$query    # coerced: character chromosome, numeric position, uppercase alleles
  ref <- m$target   # filtered + coerced; qc$target_idx indexes THIS table
  ri  <- qc$target_idx

  out <- pgs

  # Pass through any other reference columns (rsid, INFO score, ...) verbatim;
  # these aren't orientation-dependent so no flip is applicable.
  extra_cols <- setdiff(names(ref), .VM_REQ_COLS)
  add <- list()
  for (col in extra_cols) add[[paste0("ref_", col)]] <- ref[[col]][ri]

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
    n_ref_matches  = qc$n_target_matches,
    # The raw, unaligned reference values for matched variants. These are
    # included in the matching engine's required columns, so they are not
    # covered by the pass-through above, but can be useful for manual
    # verification.
    ref_position_matched        = ref$position[ri],
    ref_effect_allele           = ref$effect_allele[ri],
    ref_other_allele            = ref$other_allele[ri],
    ref_effect_allele_frequency = ref$effect_allele_frequency[ri]
  )

  # Only allele letters are changed, and only where the reference codes the site
  # on the opposite strand from the PGS. plink2 --score identifies the effect
  # allele by matching letters against the genotype file, so as long as
  # harmonized_effect_allele correctly names one of the two true alleles at the
  # site, direction is resolved correctly without re-signing the weight.
  #
  # Note that harmonized alleles are taken from the matched reference row,
  # rather than from the PGS. Since the variant matching engine accounts for
  # indel padding differences, allele representations in the PGS may not match
  # those in the reference despite being equivalent.
  #
  # Alleles on rows with no resolved orientation get set to NA instead of being
  # passed through. Columns are built by subset assignment rather than ifelse()
  # so that they stay as character type even when nothing resolved.
  rea_v <- ref$effect_allele[ri]
  roa_v <- ref$other_allele[ri]
  swap  <- !is.na(qc$allele_flipped) & qc$allele_flipped
  hm_ea <- rea_v; hm_ea[swap] <- roa_v[swap]
  hm_oa <- roa_v; hm_oa[swap] <- rea_v[swap]
  hm_ea[is.na(qc$allele_flipped)] <- NA_character_
  hm_oa[is.na(qc$allele_flipped)] <- NA_character_
  own$harmonized_effect_allele <- hm_ea
  own$harmonized_other_allele  <- hm_oa

  # The weight on unusable rows gets an NA too
  w <- pgs$effect_weight; w[!qc$usable] <- NA_real_
  own$harmonized_weight <- w

  # Merge, with this function's own output taking precedence over a colliding
  # pass-through column and over a colliding PGS column, and report both cases.
  # own[] assigns in place where the name already exists and appends otherwise,
  # so the column order is the pass-through block followed by the QC block.
  .vm_warn_collisions(names(out), names(add), names(own), "pgs", "ref")
  add[names(own)] <- own
  for (nm in names(add)) out[[nm]] <- add[[nm]]

  # Sort by chromosome and position, and then reset row names
  out <- out[.vm_order_rows(out$chromosome, out$position), , drop = FALSE]
  rownames(out) <- NULL
  out
}
