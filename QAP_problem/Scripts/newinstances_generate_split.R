#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
}))

newinstances_root_default <- function() {
  normalizePath(file.path(qap_project_root(), "NewInstances"), winslash = "/", mustWork = TRUE)
}

newinstances_split_paths <- function(split_seed,
                                     out_root = file.path(newinstances_root_default(), "splits")) {
  out_dir <- file.path(out_root, paste0("seed_", as.integer(split_seed)))
  list(
    dir = out_dir,
    train = file.path(out_dir, "train_instances.txt"),
    test = file.path(out_dir, "test_instances.txt"),
    manifest = file.path(out_dir, "split_manifest.csv"),
    summary = file.path(out_dir, "split_summary.csv")
  )
}

list_newinstances_by_group <- function(newinstances_root = newinstances_root_default()) {
  newinstances_root <- normalizePath(newinstances_root, winslash = "/", mustWork = TRUE)
  files <- list.files(newinstances_root, pattern = "\\.dat$", recursive = TRUE, full.names = TRUE)
  files <- normalizePath(sort(files), winslash = "/", mustWork = TRUE)
  if (length(files) == 0L) stop("No .dat files found under: ", newinstances_root)

  rel <- sub(paste0("^", gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", newinstances_root), "/"), "", files)
  source_group <- sub("/.*$", "", rel)

  data.frame(
    source_group = source_group,
    instance = files,
    stringsAsFactors = FALSE
  )
}

generate_newinstances_split <- function(split_seed = 11L,
                                        train_ratio = 0.5,
                                        newinstances_root = newinstances_root_default(),
                                        out_root = file.path(newinstances_root_default(), "splits")) {
  df <- list_newinstances_by_group(newinstances_root = newinstances_root)
  groups <- unique(df$source_group)

  set.seed(as.integer(split_seed))
  parts <- lapply(groups, function(g) {
    sub <- df[df$source_group == g, , drop = FALSE]
    n_total <- nrow(sub)
    n_train <- max(1L, min(n_total - 1L, floor(n_total * train_ratio)))
    train_idx <- sort(sample.int(n_total, size = n_train, replace = FALSE))
    sub$split <- "test"
    sub$split[train_idx] <- "train"
    sub
  })
  manifest <- do.call(rbind, parts)
  rownames(manifest) <- NULL

  paths <- newinstances_split_paths(split_seed = split_seed, out_root = out_root)
  dir.create(paths$dir, recursive = TRUE, showWarnings = FALSE)

  train_instances <- manifest$instance[manifest$split == "train"]
  test_instances <- manifest$instance[manifest$split == "test"]

  if (length(train_instances) > 1L) {
    set.seed(as.integer(split_seed) + 10001L)
    train_instances <- sample(train_instances, length(train_instances), replace = FALSE)
  }
  if (length(test_instances) > 1L) {
    set.seed(as.integer(split_seed) + 20001L)
    test_instances <- sample(test_instances, length(test_instances), replace = FALSE)
  }

  writeLines(train_instances, paths$train)
  writeLines(test_instances, paths$test)
  utils::write.csv(manifest, paths$manifest, row.names = FALSE)

  summary_df <- aggregate(
    as.integer(rep(1L, nrow(manifest))),
    by = list(source_group = manifest$source_group, split = manifest$split),
    FUN = sum
  )
  names(summary_df)[names(summary_df) == "x"] <- "n_instances"
  summary_df$split_seed <- as.integer(split_seed)
  summary_df$train_ratio <- as.numeric(train_ratio)
  summary_df <- summary_df[, c("split_seed", "source_group", "split", "n_instances", "train_ratio")]

  utils::write.csv(summary_df, paths$summary, row.names = FALSE)

  list(
    train = paths$train,
    test = paths$test,
    manifest = paths$manifest,
    summary = paths$summary,
    manifest_df = manifest,
    summary_df = summary_df
  )
}

if (sys.nframe() == 0L) {
  res <- generate_newinstances_split()
  cat("Saved split files under:\n")
  cat(dirname(res$train), "\n")
  print(res$summary_df)
}
