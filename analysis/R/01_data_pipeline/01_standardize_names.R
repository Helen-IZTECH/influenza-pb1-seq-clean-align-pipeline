# 01_standardize_names.R
# Hierarchical folder + FASTA file standardization for the PB1 dataset.
# Walks a raw dataset folder, detects (year, continent, processing status)
# from each file's path/name, and copies it into a clean tree:
#   PB1_DATASET_STANDARDIZED / PB1_{YEAR} / H1N1_PB1_{CONTINENT}_{YEAR} / files
# Originals are never modified (copy only).

# ---- CONFIG (edit these two paths) ------------------------------------------
base_dir   <- "data/PB1_DATASET"               # raw input tree
output_dir <- "data/PB1_DATASET_STANDARDIZED"  # standardized output tree
dry_run    <- FALSE                            # TRUE = simulate, print only, copy nothing

# ---- Continent name variants -> canonical tag ------------------------------
continents_map <- list(
  "S\\.AMERICA" = "SAMERICA",
  "S_AMERICA"  = "SAMERICA",
  "SAMERICA"   = "SAMERICA",
  "N\\.AMERICA" = "NAMERICA",
  "N_AMERICA"  = "NAMERICA",
  "NAMERICA"   = "NAMERICA",
  "AFRICA"     = "AFRICA",
  "ASIA"       = "ASIA",
  "EUROPE"     = "EUROPE",
  "OCEANIA"    = "OCEANIA"
)

# Detect year, continent and processing status from a file path/name
detect_metadata <- function(file_path) {
  path_upper <- toupper(file_path)
  file_upper <- toupper(basename(file_path))

  # 1. Year (2009-2025)
  year_match <- regmatches(path_upper, regexpr("(?<![0-9])(200[9]|201[0-9]|202[0-5])(?![0-9])", path_upper, perl = TRUE))
  year <- ifelse(length(year_match) > 0, year_match, "UNKNOWN_YEAR")

  # 2. Continent
  continent <- "UNKNOWN_CONTINENT"
  for (pattern in names(continents_map)) {
    if (grepl(pattern, path_upper)) {
      continent <- continents_map[[pattern]]
      break
    }
  }

  # 3. Processing status / stage
  if (grepl("SAMPLED[_-]?([0-9]+)", file_upper)) {
    sample_tag <- regmatches(file_upper, regexpr("SAMPLED[_-]?([0-9]+)", file_upper))
    status <- gsub("[_-]", "", sample_tag)              # e.g. SAMPLED300
  } else if (grepl("RAW", file_upper) || grepl("00_RAW", path_upper)) {
    status <- "RAW"
  } else if (grepl("FINAL|CLEAN", file_upper)) {
    status <- "CLEAN"
  } else if (grepl("UNALIGNED", file_upper)) {
    status <- "CDS_UNALIGNED"
  } else if (grepl("ALIGNED|MAFFT|ALIVIEW|ALGN", file_upper)) {
    status <- "ALIGNED"
  } else {
    status <- "PROCESSED"
  }

  return(list(year = year, continent = continent, status = status))
}

run_standardization <- function() {
  if (!dir.exists(base_dir)) {
    stop(paste("ERROR: source directory not found ->", base_dir))
  }

  all_files <- list.files(
    path = base_dir,
    pattern = "\\.(fasta|fa|fna|fas)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  cat(paste("Total files found:", length(all_files), "\n"))
  cat(paste("Mode:", ifelse(dry_run, "DRY RUN (simulation)", "REAL COPY"), "\n"))
  cat(paste(rep("-", 80), collapse = ""), "\n\n")

  success_count <- 0

  for (src_path in all_files) {
    meta <- detect_metadata(src_path)
    ext  <- tolower(paste0(".", tools::file_ext(src_path)))

    # Target folder: {output}/PB1_{YEAR}/H1N1_PB1_{CONTINENT}_{YEAR}
    target_folder <- file.path(
      output_dir,
      paste0("PB1_", meta$year),
      sprintf("H1N1_PB1_%s_%s", meta$continent, meta$year)
    )

    if (!dry_run && !dir.exists(target_folder)) {
      dir.create(target_folder, recursive = TRUE)
    }

    # Standard file name: STATUS_PB1_H1N1_YEAR_CONTINENT.fasta
    new_filename <- sprintf("%s_PB1_H1N1_%s_%s%s", meta$status, meta$year, meta$continent, ext)
    dest_path <- file.path(target_folder, new_filename)

    # Collision guard (same name in same folder -> append a counter)
    counter <- 1
    while (file.exists(dest_path) && !dry_run) {
      new_filename <- sprintf("%s_PB1_H1N1_%s_%s_%d%s", meta$status, meta$year, meta$continent, counter, ext)
      dest_path <- file.path(target_folder, new_filename)
      counter <- counter + 1
    }

    if (dry_run) {
      cat(sprintf("[TEST] %s\n       ---> %s\n", basename(src_path),
                  file.path(basename(dirname(target_folder)), basename(target_folder), new_filename)))
    } else {
      file.copy(from = src_path, to = dest_path, overwrite = FALSE)
      cat(sprintf("[COPIED] %s / %s\n", basename(target_folder), new_filename))
      success_count <- success_count + 1
    }
  }

  if (!dry_run) {
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat(sprintf("Done: %d files copied into the standardized tree.\n", success_count))
    cat(sprintf("Location: %s\n", output_dir))
  }
}

run_standardization()
