# Machine and interpreter capture for the R benchmark drivers (NV-08).
#
# The archived timing tables carry no hardware field, which is why the manuscript
# may only say "a single Apple silicon laptop" and must not assert a CPU model or a
# memory size for them.  This module closes that gap for *future* runs, mirroring
# benchmarks/machine_info.py on the Python side: the same column names
# (os_platform, os_machine, cpu, interpreter, total_memory_bytes, total_memory_gb)
# so a reader can line a `_py` table up with its `_r` counterpart and see that both
# halves came from the same box.
#
# It changes nothing about what is measured: the values are read once, before any
# timing loop, from R's own version and platform objects and from the operating
# system's memory query.  No clock is started, no repeat is added, and no published
# number moves.  Base R only, like the drivers that source it.
#
# Usage from a driver:
#   source(file.path(<this dir>, "machine_info.R"))
#   MI <- machine_fields()
#   cat(describe_machine(MI), "\n")

machine_fields <- function() {
  # Single numeric or NA, never an error: a driver must still produce its table on a
  # machine whose memory query is refused.
  memory <- tryCatch(physical_memory_bytes(), error = function(e) NA_real_)
  list(
    os_platform = paste(
      R.version.string,
      Sys.info()[["sysname"]],
      Sys.info()[["release"]],
      sep = " | "
    ),
    os_machine = Sys.info()[["machine"]],
    cpu = cpu_label(),
    interpreter = paste0("R ", getRversion()),
    total_memory_bytes = memory,
    total_memory_gb = if (is.na(memory)) NA_real_ else round(memory / 1e9, 1)
  )
}

physical_memory_bytes <- function() {
  sysname <- Sys.info()[["sysname"]]
  if (identical(sysname, "Darwin")) {
    out <- system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = FALSE)
    return(as.numeric(trimws(out[[1]])))
  }
  if (identical(sysname, "Linux")) {
    lines <- readLines("/proc/meminfo", warn = FALSE)
    total <- grep("^MemTotal:", lines, value = TRUE)
    if (length(total) == 0) return(NA_real_)
    kb <- as.numeric(gsub("[^0-9]", "", total[[1]]))
    return(kb * 1024)
  }
  # Windows and anything else: no cheap, dependency-free query.
  NA_real_
}

cpu_label <- function() {
  # R has no processor-name primitive in base; the core count is what a reader can
  # actually compare against the Python side, so name the vendor only when it is free.
  info <- Sys.info()
  cores <- tryCatch({
    sysname <- info[["sysname"]]
    if (identical(sysname, "Darwin")) {
      as.integer(trimws(system2("sysctl", c("-n", "hw.ncpu"), stdout = TRUE,
                                stderr = FALSE)[[1]]))
    } else if (identical(sysname, "Linux")) {
      parallel::detectCores()
    } else {
      parallel::detectCores()
    }
  }, error = function(e) NA_integer_)
  if (is.na(cores)) return("unknown")
  paste0(info[["sysname"]], " ", cores, " cores")
}

describe_machine <- function(fields) {
  sprintf(
    "measured on: %s | %s | %s | %s GB physical RAM | %s",
    fields$os_platform,
    fields$os_machine,
    fields$cpu,
    if (is.na(fields$total_memory_gb)) "unknown" else fields$total_memory_gb,
    fields$interpreter
  )
}
