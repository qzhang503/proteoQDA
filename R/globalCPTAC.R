#' Copy Mascot \code{.csv} files
#'
#' \code{copy_globalCPTAC} copies the \code{.csv} files from Mascot searches to
#' a target directory.
#' @export
copy_globalCPTAC <- function(dat_dir) {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  filelist <- c("F003481", "F003485", "F003486", "F003487", "F003488", "F003510")
  data(list = filelist, package = "proteoQDA", envir = environment())

  for (i in seq_along(filelist)) {
    fileConn <- file(file.path(dat_dir, paste0(filelist[i], ".csv")))
    df <- get(filelist[i])

    writeLines(df, fileConn)
    close(fileConn)
  }
}


#' Copy \code{expt_smry.xlsx}
#'
#' \code{copy_expt} copies the \code{expt_smry.xlsx} to a target directory.
#' @export
copy_expt <- function(dat_dir) {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  filepath <- system.file("extdata", "expt_smry.xlsx", package = "proteoQ")
  if (nchar(filepath) == 0) stop("Load `library(proteoQ)` first.")
  filepath <- gsub("/", "\\\\", filepath)
  file.copy(from = filepath, to = file.path(dat_dir, "expt_smry.xlsx"))
}


#' Copy \code{frac_smry.xlsx}
#'
#' \code{copy_frac} copies the \code{expt_smry.xlsx} to a target directory.
#' @export
copy_frac <- function(dat_dir) {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  filepath <- system.file("extdata", "frac_smry.xlsx", package = "proteoQ")
  if (nchar(filepath) == 0) stop("Load `library(proteoQ)` first.")
  filepath <- gsub("/", "\\\\", filepath)
  file.copy(from = filepath, to = file.path(dat_dir, "frac_smry.xlsx"))
}

