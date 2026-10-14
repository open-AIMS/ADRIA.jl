# Per-reef GAMs of the wet-season chlorophyll-a anomaly on Nino 3.4, to back-cast the
# chl-a index to years without satellite data (lead's decision, 2026-10-15, ninth round: the
# lead's earlier workflow predicted unsurveyed years from the Nino 3.4 - chl-a anomaly
# relationship with a GAM per reef).
#
# Inputs:
#   - extract_chla_wet_season.py output (MODIS-Aqua OC3, wet seasons 2003-2024; the IMOS
#     product ends 2024-12-25, so 2025-2026 are also predicted);
#   - NOAA CPC monthly Nino 3.4 anomaly (https://www.cpc.ncep.noaa.gov/data/indices/sstoi.indices).
# ENSO predictors per season `year`: mean Nino 3.4 anomaly over Dec(year-1)-Apr(year) ("wet")
# and Sep-Nov(year-1) ("pre"). Per reef and response (log anomaly of the seasonal mean and of
# the seasonal 90th percentile): mgcv::gam(y ~ s(wet, k = 4)) and gam(y ~ s(pre, k = 4)); the
# predictor is chosen once for the whole domain by leave-one-season-out cross-validated R^2
# (pooled over reefs). Predictions for seasons without satellite values are flagged.
#
#   Rscript build_chla_enso_backcast.R [<chla dir>] [<output dir>]

suppressPackageStartupMessages({library(mgcv)})
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "..", ".."))
chla_dir <- if (length(args) >= 1) args[1] else file.path(root, "sandbox/data/derived/chla_wet_season_lizard")
out_dir <- if (length(args) >= 2) args[2] else file.path(root, "sandbox/data/derived/chla_index_lizard_1991_2026")
if (dir.exists(out_dir)) stop("Refusing to overwrite ", out_dir)

chla <- read.csv(file.path(chla_dir, "chla_wet_season_reef_year.csv"), colClasses = c(RME_UNIQUE_ID = "character"))
url <- "https://www.cpc.ncep.noaa.gov/data/indices/sstoi.indices"
raw <- readLines(url)
nino <- read.table(text = raw, header = FALSE, skip = 1)[, c(1, 2, 10)]
names(nino) <- c("yr", "mon", "anom")
season_mean <- function(y, months_prev, months_cur) {
  v <- c(nino$anom[nino$yr == y - 1 & nino$mon %in% months_prev], nino$anom[nino$yr == y & nino$mon %in% months_cur])
  if (length(v) < length(months_prev) + length(months_cur)) NA else mean(v)
}
years <- 1991:2026
enso <- data.frame(year = years,
                   wet = sapply(years, season_mean, months_prev = 12, months_cur = 1:4),
                   pre = sapply(years, season_mean, months_prev = 9:11, months_cur = integer(0)))

d <- merge(chla, enso, by = "year")
obs <- d[is.finite(d$log_anom_mean) & is.finite(d$log_anom_p90), ]
cv_r2 <- function(response, predictor) {
  pred <- rep(NA_real_, nrow(obs))
  for (yr in unique(obs$year)) {
    for (reef in unique(obs$RME_UNIQUE_ID)) {
      train <- obs[obs$RME_UNIQUE_ID == reef & obs$year != yr, ]
      test <- which(obs$RME_UNIQUE_ID == reef & obs$year == yr)
      if (length(test) == 0 || nrow(train) < 8) next
      f <- gam(as.formula(paste(response, "~ s(", predictor, ", k = 4)")), data = train)
      pred[test] <- predict(f, newdata = obs[test, ])
    }
  }
  ok <- is.finite(pred)
  y <- obs[[response]][ok]
  1 - sum((y - pred[ok])^2) / sum((y - mean(y))^2)
}
scores <- expand.grid(response = c("log_anom_mean", "log_anom_p90"), predictor = c("wet", "pre"), stringsAsFactors = FALSE)
scores$cv_r2 <- mapply(cv_r2, scores$response, scores$predictor)
print(scores)
choice <- sapply(c("log_anom_mean", "log_anom_p90"), function(r) {
  s <- scores[scores$response == r, ]; s$predictor[which.max(s$cv_r2)]
})

grid <- expand.grid(RME_UNIQUE_ID = unique(chla$RME_UNIQUE_ID), year = years, stringsAsFactors = FALSE)
grid <- merge(grid, enso, by = "year")
grid <- merge(grid, chla[, c("RME_UNIQUE_ID", "reef_id", "year", "mean", "p90", "log_anom_mean", "log_anom_p90", "n")],
              by = c("RME_UNIQUE_ID", "year"), all.x = TRUE)
for (r in names(choice)) {
  grid[[paste0(r, "_gam")]] <- NA_real_
  for (reef in unique(grid$RME_UNIQUE_ID)) {
    train <- obs[obs$RME_UNIQUE_ID == reef, ]
    if (nrow(train) < 8) next
    f <- gam(as.formula(paste(r, "~ s(", choice[[r]], ", k = 4)")), data = train)
    rows <- grid$RME_UNIQUE_ID == reef
    grid[rows, paste0(r, "_gam")] <- predict(f, newdata = grid[rows, ])
  }
  grid[[paste0(r, "_index")]] <- ifelse(is.finite(grid[[r]]), grid[[r]], grid[[paste0(r, "_gam")]])
}
grid$source <- ifelse(is.finite(grid$log_anom_mean), "MODIS-Aqua OC3", "Nino3.4 GAM")
grid <- grid[order(grid$RME_UNIQUE_ID, grid$year), ]
dir.create(out_dir, recursive = TRUE)
write.csv(grid, file.path(out_dir, "chla_index_reef_year.csv"), row.names = FALSE)
write.csv(enso, file.path(out_dir, "nino34_seasons.csv"), row.names = FALSE)
write.csv(scores, file.path(out_dir, "gam_cv_r2.csv"), row.names = FALSE)
writeLines(c(paste("nino34_source:", url), paste("retrieved:", Sys.time()),
             paste("chla_input:", file.path(chla_dir, "chla_wet_season_reef_year.csv")),
             paste("chla_input_md5:", tools::md5sum(file.path(chla_dir, "chla_wet_season_reef_year.csv"))),
             paste("predictor_choice:", paste(names(choice), choice, sep = "=", collapse = "; ")),
             "model: mgcv::gam(y ~ s(x, k = 4)) per reef; choice by pooled leave-one-season-out CV R^2",
             paste("script_md5:", tools::md5sum(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))),
           file.path(out_dir, "provenance.txt"))
agg <- aggregate(cbind(log_anom_mean_index, log_anom_p90_index) ~ year + source, data = grid, FUN = median)
print(agg[order(agg$year), ], row.names = FALSE)
