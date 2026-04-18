PENALTY <- 1e20

qap_find_scripts_dir_from <- function(start) {
  if (!is.character(start) || length(start) == 0L || is.na(start) || !nzchar(start)) {
    return(NA_character_)
  }

  current <- normalizePath(start[[1L]], winslash = "/", mustWork = TRUE)
  if (!dir.exists(current)) {
    current <- dirname(current)
  }

  repeat {
    candidate <- file.path(current, "QAP_problem", "Scripts")
    if (dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
    if (basename(current) == "Scripts" && basename(dirname(current)) == "QAP_problem") {
      return(normalizePath(current, winslash = "/", mustWork = TRUE))
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }

  NA_character_
}

qap_rstudio_editor_path <- function() {
  if (!interactive() || !requireNamespace("rstudioapi", quietly = TRUE)) {
    return(NA_character_)
  }

  ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
  if (is.null(ctx)) {
    return(NA_character_)
  }

  path <- ctx$path
  if (!is.character(path) || length(path) == 0L || !nzchar(path) || !file.exists(path)) {
    return(NA_character_)
  }

  normalizePath(path[[1L]], winslash = "/", mustWork = TRUE)
}

qap_script_path <- function() {
  editor_path <- qap_rstudio_editor_path()
  if (is.character(editor_path) && !is.na(editor_path)) {
    return(editor_path)
  }

  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", file_arg, value = TRUE)
  if (length(file_flag) > 0L) {
    return(normalizePath(sub("^--file=", "", file_flag[[1L]]), winslash = "/", mustWork = TRUE))
  }

  frames <- sys.frames()
  ofiles <- vapply(frames, function(env) {
    if (exists("ofile", envir = env, inherits = FALSE)) {
      as.character(get("ofile", envir = env, inherits = FALSE))
    } else {
      NA_character_
    }
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles) > 0L) {
    return(normalizePath(tail(ofiles, 1L), winslash = "/", mustWork = TRUE))
  }

  NA_character_
}

qap_scripts_dir <- function() {
  if (exists("script_dir", envir = .GlobalEnv, inherits = FALSE)) {
    existing <- get("script_dir", envir = .GlobalEnv, inherits = FALSE)
    if (is.character(existing) && length(existing) > 0L && nzchar(existing[[1L]]) && dir.exists(existing[[1L]])) {
      candidate <- normalizePath(existing[[1L]], winslash = "/", mustWork = TRUE)
      if (basename(candidate) == "Scripts" && basename(dirname(candidate)) == "QAP_problem") {
        return(candidate)
      }
    }
  }

  script_file <- qap_script_path()
  if (is.character(script_file) && !is.na(script_file) && file.exists(script_file)) {
    resolved <- qap_find_scripts_dir_from(script_file)
    if (is.character(resolved) && !is.na(resolved)) {
      return(resolved)
    }
  }

  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  resolved <- qap_find_scripts_dir_from(cwd)
  if (is.character(resolved) && !is.na(resolved)) {
    return(resolved)
  }

  stop("Cannot resolve QAP_problem/Scripts. In RStudio, set working directory to QAP_problem/Scripts or the repo root.")
}

qap_project_root <- function(script_dir = NULL) {
  if (is.null(script_dir)) script_dir <- qap_scripts_dir()
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
}

extract_qap_cost <- function(lines) {
  trimmed <- trimws(lines)
  numeric_lines <- grep("^-?[0-9]+(\\.[0-9]+)?$", trimmed, value = TRUE)
  if (length(numeric_lines) == 0L) {
    return(NA_real_)
  }
  as.numeric(numeric_lines[[length(numeric_lines)]])
}

run_qap_solver <- function(instance, a, p, l, t = 5, seed = NULL, solver = NULL, echo_logs = FALSE) {
  if (is.null(solver)) {
    solver <- file.path(qap_project_root(), "Prog", "lsmcqap")
  }
  if (!file.exists(solver)) {
    return(list(cost = PENALTY, status = 127L, output = "solver_missing"))
  }

  args <- c("-i", instance, "-t", format(t, scientific = FALSE), "-a", a, "-p", p, "-l", l)
  if (!is.null(seed) && !is.na(seed)) {
    args <- c(args, "-s", seed)
  }

  output <- tryCatch(
    system2(solver, args = args, stdout = TRUE, stderr = TRUE),
    warning = function(w) {
      attr(w, "status") <- 1L
      structure(conditionMessage(w), status = 1L)
    },
    error = function(e) {
      structure(conditionMessage(e), status = 1L)
    }
  )

  if (is.character(output)) {
    status <- attr(output, "status")
    if (is.null(status) || is.na(status)) status <- 0L
    cost <- extract_qap_cost(output)
    if (!is.finite(cost)) cost <- PENALTY
    if (echo_logs) {
      cat(paste(output, collapse = "\n"), "\n", sep = "")
    }
    return(list(cost = cost, status = status, output = output))
  }

  list(cost = PENALTY, status = 1L, output = "unexpected_runner_state")
}
