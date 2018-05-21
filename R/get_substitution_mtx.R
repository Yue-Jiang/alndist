#' Get the substitution matrix.
#'
#' @param mtx.name A character string. The name of substitution matrix.
#' Available choices are "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
#' "PAM30", "PAM40", "PAM70", "PAM120", "PAM250".
#' @param gap.char A character for gap. Default is "-".
#' @param gap.gap A numeric value. Custom score for gap vs gap that overrides the default in the substitution matrix. Do not override if set to NA.
#' @param gap.non.gap A numeric value. Custom score for gap vs gap that overrides the default in the substitution matrix. Do not override if set to NA.
#'
#' @return A substitution matrix.
#' @examples
#' \dontrun{
#' get_substitution_mtx("BLOSUM62", gap.char="*")
#' }

get_substitute_mtx <- function(mtx.name="BLOSUM62", gap.char="-", gap.gap=NA, gap.non.gap=NA) {
  if (is.na(gap.char)) {
    warning("gap.char is set to '-'")
    gap.char <- "-"
  }
  mtx <- switch(
    toupper(mtx.name),
    "BLOSUM45"  = readRDS(system.file("extdata/BLOSUM45.rds", package="alndist")),
    "BLOSUM50"  = readRDS(system.file("extdata/BLOSUM50.rds", package="alndist")),
    "BLOSUM62"  = readRDS(system.file("extdata/BLOSUM62.rds", package="alndist")),
    "BLOSUM80"  = readRDS(system.file("extdata/BLOSUM80.rds", package="alndist")),
    "BLOSUM100" = readRDS(system.file("extdata/BLOSUM100.rds", package="alndist")),
    "PAM30"     = readRDS(system.file("extdata/PAM30.rds", package="alndist")),
    "PAM40"     = readRDS(system.file("extdata/PAM40.rds", package="alndist")),
    "PAM70"     = readRDS(system.file("extdata/PAM70.rds", package="alndist")),
    "PAM120"    = readRDS(system.file("extdata/PAM120.rds", package="alndist")),
    "PAM250"    = readRDS(system.file("extdata/PAM250.rds", package="alndist")),
    stop(sprintf("%s is not a valid substitution matrix name", mtx.name))
  )
  if (!is.na(gap.gap)) {
    if (!is.numeric(gap.gap)) stop("gap.gap can only be numeric")
    mtx[rownames(mtx) == "-", colnames(mtx) == "-"] <- gap.gap
  }
  if (!is.na(gap.non.gap)) {
    if (!is.numeric(gap.non.gap)) stop("gap.non.gap can only be numeric")
    mtx[rownames(mtx) == "-", colnames(mtx) != "-"] <- gap.non.gap
    mtx[rownames(mtx) != "-", colnames(mtx) == "-"] <- gap.non.gap
  }
  if (gap.char != "-") {
    if (gap.char %in% colnames(mtx)) stop("gap.char has already been used as a non-gap character, choose another one")
    colnames(mtx)[colnames(mtx) == "-"] <- gap.char
    rownames(mtx)[rownames(mtx) == "-"] <- gap.char
  }
  return(mtx)
}
