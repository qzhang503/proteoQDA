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



#' Copy an \code{expt_smry...} file to \code{dat_dir}
#'
#' \code{copy_expt} copies a system file of \code{expt_smry...} to the target
#' directory specified by \code{dat_dir}.
#' @export
copy_expt <- function(dat_dir, from = "expt_smry.xlsx", to = "expt_smry.xlsx") {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  filepath <- system.file("extdata", from, package = "proteoQDA")
  filepath <- gsub("/", "\\\\", filepath)
  file.copy(from = filepath, to = file.path(dat_dir, to))
}


#' Copy an \code{expt_smry...} file to \code{dat_dir}
#'
#' @export
cptac_expt_1 <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_cptac_gl.xlsx", "expt_smry.xlsx")
}


#' Copy an \code{expt_smry...} file to \code{dat_dir}
#'
#' @export
cptac_expt_2 <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_cptac_cmbn.xlsx", "expt_smry.xlsx")
}


#' Copy an \code{expt_smry...} file to \code{dat_dir}
#'
#' @export
cptac_expt_ref_w2 <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_ref_w2.xlsx", "expt_smry_ref_w2.xlsx")
}


#' Copy \code{expt_smry.xlsx}
#'
#' @export
expt_smry_ref_w2_w16 <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_ref_w2_w16.xlsx", "expt_smry_ref_w2_w16.xlsx")
}



#' Copy \code{frac_smry.xlsx}
#'
#' \code{copy_frac} copies the \code{frac_smry.xlsx} to a target directory.
#' @export
copy_frac <- function(dat_dir, from = "frac_smry.xlsx", to = "frac_smry.xlsx") {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  filepath <- system.file("extdata", from, package = "proteoQDA")
  filepath <- gsub("/", "\\\\", filepath)
  file.copy(from = filepath, to = file.path(dat_dir, to))
}


#' Copy \code{frac_smry.xlsx}
#'
#' @export
cptac_frac_1 <- function(dat_dir) {
  copy_frac(dat_dir, "frac_smry_cptac_gl.xlsx", "frac_smry.xlsx")
}


#' Copy \code{frac_smry.xlsx}
#'
#' @export
cptac_frac_2 <- function(dat_dir) {
  copy_frac(dat_dir, "frac_smry_cptac_cmbn.xlsx", "frac_smry.xlsx")
}



