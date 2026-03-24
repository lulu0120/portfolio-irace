#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
}))

script_dir <- local({
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", file_arg, value = TRUE)
  if (length(file_flag) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_flag[[1L]]), winslash = "/", mustWork = TRUE))
  } else {
    frames <- sys.frames()
    ofiles <- vapply(frames, function(env) {
      if (exists("ofile", envir = env, inherits = FALSE)) as.character(get("ofile", envir = env, inherits = FALSE)) else NA_character_
    }, character(1))
    ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
    if (length(ofiles) > 0L) {
      dirname(normalizePath(tail(ofiles, 1L), winslash = "/", mustWork = TRUE))
    } else {
      normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    }
  }
})
source(file.path(script_dir, "qap_common.R"), local = .GlobalEnv)
script_dir <- qap_scripts_dir()

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript qap_runner.R --instance PATH --a INT --p INT --l INT [--t NUM] [--seed INT]\n",
    "  Rscript qap_runner.R --self-test [--n INT]\n",
    sep = ""
  )
}

parse_args <- function(args) {
  opts <- list(
    instance = NULL,
    a = NULL,
    p = NULL,
    l = NULL,
    t = 5,
    seed = NULL,
    self_test = FALSE,
    n = 5L
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--instance") {
      i <- i + 1L
      opts$instance <- args[[i]]
    } else if (arg == "--a") {
      i <- i + 1L
      opts$a <- as.integer(args[[i]])
    } else if (arg == "--p") {
      i <- i + 1L
      opts$p <- as.integer(args[[i]])
    } else if (arg == "--l") {
      i <- i + 1L
      opts$l <- as.integer(args[[i]])
    } else if (arg == "--t") {
      i <- i + 1L
      opts$t <- as.numeric(args[[i]])
    } else if (arg == "--seed") {
      i <- i + 1L
      opts$seed <- as.integer(args[[i]])
    } else if (arg == "--self-test") {
      opts$self_test <- TRUE
    } else if (arg == "--n") {
      i <- i + 1L
      opts$n <- as.integer(args[[i]])
    } else if (arg %in% c("-h", "--help")) {
      usage()
      quit(status = 0L)
    } else {
      stop("Unknown argument: ", arg)
    }
    i <- i + 1L
  }

  opts
}

run_self_test <- function(n) {
  instances_root <- file.path(qap_project_root(), "instances", "100")
  family_dirs <- c("GridRandom", "GridStructured", "RandomRandom", "RandomStructured")
  files <- unlist(lapply(family_dirs, function(family) {
    list.files(file.path(instances_root, family), pattern = "\\.dat$", full.names = TRUE)
  }), use.names = FALSE)
  files <- sort(files)
  if (length(files) == 0L) stop("No candidate instances found for self-test")

  combos <- data.frame(
    a = c(0L, 2L, 7L, 1L, 5L, 3L, 6L, 4L, 0L, 7L),
    p = c(0L, 3L, 9L, 1L, 5L, 8L, 2L, 7L, 4L, 6L),
    l = c(0L, 1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 0L)
  )

  n <- min(n, nrow(combos), length(files))
  for (i in seq_len(n)) {
    res <- run_qap_solver(
      instance = files[[i]],
      a = combos$a[[i]],
      p = combos$p[[i]],
      l = combos$l[[i]],
      t = 5,
      seed = i
    )
    cat(sprintf(
      "[%d] status=%s cost=%s instance=%s a=%d p=%d l=%d\n",
      i, res$status, format(res$cost, scientific = FALSE), basename(files[[i]]),
      combos$a[[i]], combos$p[[i]], combos$l[[i]]
    ))
  }
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

if (isTRUE(opts$self_test)) {
  run_self_test(opts$n)
  quit(status = 0L)
}

required <- c("instance", "a", "p", "l")
missing_required <- required[vapply(required, function(name) is.null(opts[[name]]), logical(1))]
if (length(missing_required) > 0L) {
  usage()
  stop("Missing required arguments: ", paste(missing_required, collapse = ", "))
}

instance_path <- tryCatch(
  normalizePath(opts$instance, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)

if (is.na(instance_path)) {
  cat(format(PENALTY, scientific = FALSE), "\n", sep = "")
  quit(status = 1L)
}

res <- run_qap_solver(
  instance = instance_path,
  a = opts$a,
  p = opts$p,
  l = opts$l,
  t = opts$t,
  seed = opts$seed
)

cat(format(res$cost, scientific = FALSE), "\n", sep = "")
quit(status = if (is.finite(res$cost) && res$cost < PENALTY) 0L else 1L)
