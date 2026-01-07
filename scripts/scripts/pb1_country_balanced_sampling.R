## =========================================================
##  COUNTRY-BALANCED FASTA SAMPLING (ROUND-ROBIN)
##  Purpose:
##    Create a country-balanced subset from a large FASTA file
##    to reduce geographic over-representation before downstream analyses.
##
##  Output:
##    - One FASTA file (subset)
##    - In-R summary tables (country counts + flags)
## =========================================================

## ---------- (A) WORKING DIRECTORY ----------
setwd("C:/Users/HP/Desktop/PB1_DATASET/PB1_2024/H1N1_PB1_N.AMERICA_2024")

## ---------- (B) MAIN FUNCTION ----------
sample_round_robin_fill_FASTAonly <- function(
    fasta_in,
    fasta_out,
    n_total   = 300,
    start_cap = 5,
    hard_cap  = 50,
    seed      = 42,
    flag_low  = 1,
    flag_high = 15
) {
  set.seed(seed)

  ## ---- 1) Read FASTA ----
  dna <- Biostrings::readDNAStringSet(fasta_in)
  hdr <- names(dna)

  if (length(dna) < n_total) {
    stop("Total sequences ", length(dna), " < target ", n_total, ". Cannot reach n_total.")
  }

  ## ---- 2) Parse headers (assumes your header format) ----
  ## Expected to contain a strain field like: A/Togo/0629/2023
  parse_header <- function(h) {
    parts  <- stringr::str_split(h, "\\|")[[1]]
    strain <- parts[3]
    sub    <- stringr::str_split(strain, "/")[[1]]

    tibble::tibble(
      country = sub[2],
      year    = as.integer(sub[length(sub)]),
      header  = h
    )
  }

  meta <- dplyr::bind_rows(lapply(hdr, parse_header))

  ## ---- 3) Available sequences per country ----
  avail <- meta %>% dplyr::count(country, name = "available")
  available_counts <- setNames(avail$available, avail$country)

  ## ---- 4) Increase cap to the minimum value that can reach n_total ----
  total_capacity <- function(cap) sum(pmin(available_counts, cap))

  cap <- start_cap
  while (total_capacity(cap) < n_total && cap <= hard_cap) cap <- cap + 1

  if (total_capacity(cap) < n_total) {
    stop("Even hard_cap=", hard_cap, " cannot reach n_total. Increase hard_cap.")
  }

  ## ---- 5) Round-robin sampling across countries ----
  by_country <- split(meta$header, meta$country)
  by_country <- lapply(by_country, sample)  # shuffle within each country

  picked <- character(0)
  used_counts <- setNames(rep(0, length(by_country)), names(by_country))

  while (length(picked) < n_total) {
    progress <- FALSE

    for (cty in names(by_country)) {
      if (length(picked) >= n_total) break
      if (used_counts[cty] >= cap) next
      if (length(by_country[[cty]]) == 0) next

      picked <- c(picked, by_country[[cty]][1])
      by_country[[cty]] <- by_country[[cty]][-1]
      used_counts[cty] <- used_counts[cty] + 1
      progress <- TRUE
    }

    if (!progress) break
  }

  if (length(picked) != n_total) {
    stop(
      "Unexpected: selected ", length(picked), " != ", n_total,
      ". This suggests a header parse or capacity issue."
    )
  }

  ## ---- 6) Write the sampled FASTA (single output file) ----
  dna_out <- dna[picked]
  Biostrings::writeXStringSet(dna_out, fasta_out, format = "fasta")

  ## ---- 7) Build in-R sanity-check tables ----
  picked_meta <- meta %>% dplyr::filter(header %in% picked)

  country_table <- picked_meta %>%
    dplyr::count(country, sort = TRUE)

  flagged <- country_table %>%
    dplyr::mutate(
      flag = dplyr::case_when(
        n <= flag_low  ~ "LOW",
        n >= flag_high ~ "HIGH",
        TRUE           ~ ""
      )
    ) %>%
    dplyr::filter(flag != "")

  summary_stats <- country_table %>%
    dplyr::summarise(
      countries = dplyr::n(),
      min_n  = min(n),
      q1     = quantile(n, 0.25),
      median = median(n),
      mean   = mean(n),
      q3     = quantile(n, 0.75),
      max_n  = max(n)
    )

  message(
    "✅ FASTA saved: ", fasta_out,
    " | selected countries: ", summary_stats$countries,
    " | cap used: ", cap
  )

  list(
    fasta_out      = fasta_out,
    cap_used       = cap,
    country_table  = country_table,
    flagged        = flagged,
    summary        = summary_stats
  )
}

## =======================
##  (C) RUN (EXAMPLE)
## =======================

res <- sample_round_robin_fill_FASTAonly(
  fasta_in   = "N.AMERICA_2024_PB1_final.fasta",
  fasta_out  = "SAMPLED_600_N.AMERICA_2024.fasta",
  n_total    = 600,
  start_cap  = 5,
  hard_cap   = 50,
  seed       = 123,
  flag_high  = 15
)

print(res$summary)
print(res$country_table)
print(res$flagged)
View(res$country_table)
