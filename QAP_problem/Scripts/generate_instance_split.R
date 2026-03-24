#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
}))

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript generate_instance_split.R [--seed N] [--size 100] [--include-plus]\n",
    "                                   [--root PATH] [--out-dir PATH]\n",
    "\n",
    "Description:\n",
    "  Stratified random split of QAP instances into train/test by family.\n",
    "  By default, uses the four Thomas-recommended families only.\n",
    sep = ""
  )
}

parse_args <- function(args) {
  opts <- list(
    seed = 1L,
    size = "100",
    include_plus = FALSE,
    root = NULL,
    out_dir = NULL
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--seed") {
      i <- i + 1L
      opts$seed <- as.integer(args[[i]])
    } else if (arg == "--size") {
      i <- i + 1L
      opts$size <- args[[i]]
    } else if (arg == "--include-plus") {
      opts$include_plus <- TRUE
    } else if (arg == "--root") {
      i <- i + 1L
      opts$root <- args[[i]]
    } else if (arg == "--out-dir") {
      i <- i + 1L
      opts$out_dir <- args[[i]]
    } else if (arg %in% c("-h", "--help")) {
      usage()
      quit(status = 0L)
    } else {
      stop("Unknown argument: ", arg)
    }
    i <- i + 1L
  }

  if (is.na(opts$seed)) stop("Invalid --seed")
  opts
}

script_path <- function() {
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", file_arg, value = TRUE)
  if (length(file_flag) == 0L) stop("Cannot resolve script path")
  normalizePath(sub("^--file=", "", file_flag[[1L]]), winslash = "/", mustWork = TRUE)
}

repo_root <- function() {
  normalizePath(file.path(dirname(script_path()), ".."), winslash = "/", mustWork = TRUE)
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- if (is.null(opts$root)) repo_root() else normalizePath(opts$root, winslash = "/", mustWork = TRUE)
instances_root <- file.path(project_root, "instances", opts$size)
if (!dir.exists(instances_root)) stop("Instances directory not found: ", instances_root)

families <- c("GridRandom", "GridStructured", "RandomRandom", "RandomStructured")
if (opts$include_plus) {
  families <- c(families, "GridStructuredPlus", "RandomStructuredPlus")
}

split_seed_dir <- if (is.null(opts$out_dir)) {
  file.path(project_root, "Scripts", "splits", paste0("seed_", opts$seed))
} else {
  normalizePath(opts$out_dir, winslash = "/", mustWork = FALSE)
}
dir.create(split_seed_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(opts$seed)
rows <- list()

for (family in families) {
  family_dir <- file.path(instances_root, family)
  if (!dir.exists(family_dir)) stop("Missing family directory: ", family_dir)
  files <- sort(list.files(family_dir, pattern = "\\.dat$", full.names = TRUE))
  n_files <- length(files)
  if (n_files == 0L) stop("No .dat files in: ", family_dir)

  train_n <- floor(n_files / 2L)
  idx <- sample.int(n_files, size = train_n, replace = FALSE)
  split <- rep("test", n_files)
  split[idx] <- "train"

  rows[[family]] <- data.frame(
    family = family,
    split = split,
    instance = normalizePath(files, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

split_df <- do.call(rbind, rows)
split_df <- split_df[order(split_df$family, split_df$split, split_df$instance), ]

train_df <- split_df[split_df$split == "train", , drop = FALSE]
test_df <- split_df[split_df$split == "test", , drop = FALSE]

writeLines(train_df$instance, file.path(split_seed_dir, "train_instances.txt"))
writeLines(test_df$instance, file.path(split_seed_dir, "test_instances.txt"))
write.csv(split_df, file.path(split_seed_dir, "split_manifest.csv"), row.names = FALSE)

summary_df <- aggregate(instance ~ family + split, data = split_df, FUN = length)
names(summary_df)[names(summary_df) == "instance"] <- "count"
write.csv(summary_df, file.path(split_seed_dir, "split_summary.csv"), row.names = FALSE)

cat("Wrote split files to:", split_seed_dir, "\n")
print(summary_df)
