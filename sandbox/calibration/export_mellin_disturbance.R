# Export the Mellin et al. (2019) / COTS_Model disease and unknown disturbance
# coefficients and LTMP event records for the background-mortality builder
# (crash_mechanism_protocol.md, 2026-10-13).
#
# Source: github.com/sammatthews990/COTS_Model, RData/COTSMod_bckp.RData, with the model
# equations in COTSModel_COTSFunctions.R::doCoralDisturbance. Disease and unknown events
# act on log coral cover: log(C) += B.DISEASE (or B.UNKNOWN) in an event year.
# COTSMod.R simulates events at 1% incidence per cell-year each, replaced by LTMP records
# where available.
#
#   Rscript sandbox/calibration/export_mellin_disturbance.R <COTSMod_bckp.RData> <output dir>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
rdata <- args[1]
output <- args[2]
if (dir.exists(output)) stop("Refusing to overwrite ", output)
dir.create(output, recursive = TRUE)

e <- new.env()
load(rdata, envir = e)
coefficients <- data.frame(
  term = c("disease", "unknown", "storms", "bleaching", "cots"),
  log_cover_effect_mean = c(e$disease.mn.sd[1], e$unknown.mn.sd[1], e$storms.mn.sd[1],
                            e$bleaching.mn.sd[1], e$COTS.mn.sd[1]),
  log_cover_effect_sd = c(e$disease.mn.sd[2], e$unknown.mn.sd[2], e$storms.mn.sd[2],
                          e$bleaching.mn.sd[2], e$COTS.mn.sd[2]))
write.csv(coefficients, file.path(output, "mellin_coefficients.csv"), row.names = FALSE)

# Reef-year LTMP records: an event if any pixel of the reef records one; observed if
# any pixel is non-missing.
records <- function(table, kind) {
  years <- as.integer(sub("^X", "", names(table)[-(1:5)]))
  values <- as.matrix(table[, -(1:5)])
  out <- list()
  for (reef in unique(table$REEF_ID)) {
    rows <- which(table$REEF_ID == reef)
    block <- values[rows, , drop = FALSE]
    observed <- colSums(!is.na(block)) > 0
    event <- colSums(block > 0, na.rm = TRUE) > 0
    keep <- observed
    if (any(keep)) out[[length(out) + 1]] <- data.frame(
      reef_id = reef, year = years[keep], kind = kind, event = as.integer(event[keep]))
  }
  do.call(rbind, out)
}
events <- rbind(records(e$data.ltmp.disease, "disease"), records(e$data.ltmp.unknown, "unknown"))
write.csv(events, file.path(output, "mellin_ltmp_events.csv"), row.names = FALSE)

observed <- function(table) {
  v <- as.matrix(table[, -(1:5)])
  c(cell_years = sum(!is.na(v)), incidence = mean(v[!is.na(v)] > 0))
}
summary <- rbind(disease = observed(e$data.ltmp.disease), unknown = observed(e$data.ltmp.unknown))
write.csv(data.frame(kind = rownames(summary), summary, row.names = NULL),
          file.path(output, "mellin_ltmp_incidence.csv"), row.names = FALSE)
writeLines(c(
  paste0("source_rdata=", normalizePath(rdata)),
  paste0("source_rdata_md5=", unname(tools::md5sum(rdata))),
  "simulated_incidence_per_cell_year=0.01 (COTSMod.R: sampling::srswor(round(n * 0.01), n) for disease and unknown)",
  "effect=log(cover) += B in an event year (COTSModel_COTSFunctions.R::doCoralDisturbance)"),
  file.path(output, "provenance.txt"))
cat("exported", nrow(events), "LTMP reef-year records\n")
print(coefficients)
print(summary)
