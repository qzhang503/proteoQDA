#' Copy Mascot \code{.csv} files
#'
#' \code{copy_dat} copies the Mascot outputs of \code{.csv} files to a target
#' directory.
#' @export
copy_csv <- function(dat_dir, filelist) {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  data(list = filelist, package = "proteoQDA", envir = environment())

  for (i in seq_along(filelist)) {
    fileConn <- file(file.path(dat_dir, paste0(filelist[i], ".csv")))
    df <- get(filelist[i])
    writeLines(df, fileConn)
    close(fileConn)
  }
}


#' Copy Mascot \code{.csv} files
#'
#' @export
cptac_csv_1 <- function(dat_dir) {
	copy_csv(dat_dir, filelist = c("F003481", "F003485", "F003486", "F003487", "F003488", "F003510"))
}


#' Copy Mascot \code{.csv} files
#'
#' @export
cptac_csv_2 <- function(dat_dir) {
  copy_csv(dat_dir, filelist = c("F003529", "F003530", "F003531", "F003532", "F003533", "F003534"))
}



#' Copy \code{expt_smry.xlsx}
#'
#' \code{copy_expt} copies the \code{expt_smry.xlsx} to a target directory.
#' @export
copy_expt <- function(dat_dir, filename = "expt_smry.xlsx") {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  filepath <- system.file("extdata", filename, package = "proteoQDA")
  # if (nchar(filepath) == 0) stop("Load `library(proteoQ)` first.")
  filepath <- gsub("/", "\\\\", filepath)
  file.copy(from = filepath, to = file.path(dat_dir, "expt_smry.xlsx"))
}


#' Copy \code{expt_smry.xlsx}
#'
#' @export
cptac_expt_1 <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_cptac_gl.xlsx")
}


#' Copy \code{expt_smry.xlsx}
#'
#' @export
cptac_expt_2 <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_cptac_cmbn.xlsx")
}



#' Copy \code{frac_smry.xlsx}
#'
#' \code{copy_frac} copies the \code{frac_smry.xlsx} to a target directory.
#' @export
copy_frac <- function(dat_dir, filename = "frac_smry.xlsx") {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  filepath <- system.file("extdata", filename, package = "proteoQDA")
  # if (nchar(filepath) == 0) stop("Load `library(proteoQ)` first.")
  filepath <- gsub("/", "\\\\", filepath)
  file.copy(from = filepath, to = file.path(dat_dir, "frac_smry.xlsx"))
}


#' Copy \code{frac_smry.xlsx}
#'
#' @export
cptac_frac_1 <- function(dat_dir) {
  copy_frac(dat_dir, "frac_smry_cptac_gl.xlsx")
}


#' Copy \code{frac_smry.xlsx}
#'
#' @export
cptac_frac_2 <- function(dat_dir) {
  copy_frac(dat_dir, "frac_smry_cptac_cmbn.xlsx")
}



