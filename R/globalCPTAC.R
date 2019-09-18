#' Make rda files
foo <- function () {
  # dat_dir <- "c:\\The\\Mascot\\Example"
  # filelist <- c("F003590", "F003591", "F003593", "F003594", "F003595", "F003597")

  dat_dir <- "c:\\The\\Phosphopeptide\\Example"
  filelist <- c("F003607", "F003608", "F003609", "F003610", "F003611", "F003612")

  purrr::walk(filelist, ~ {
    assign(.x, readLines(file.path(dat_dir, paste0(.x, ".csv"))))
    save(list = .x, file = file.path(dat_dir, paste0(.x, ".rda")))
  })
}


#' Make rda files
foo2 <- function () {
  dat_dir <- "c:\\The\\MQ\\Example"
  filelist <- c("modificationSpecificPeptides_bi_1", "modificationSpecificPeptides_bi_2",
                "modificationSpecificPeptides_jhu_1", "modificationSpecificPeptides_jhu2",
                "modificationSpecificPeptides_pnnl_1", "modificationSpecificPeptides_pnnl_2")
  # filelist <- c("peptides_bi_1", "peptides_bi_2", "peptides_jhu_1", "peptides_jhu_2", "peptides_pnnl_1", "peptides_pnnl_2")

  purrr::walk(filelist, ~ {
    assign(.x, read.csv(file.path(dat_dir, paste0(.x, ".txt")), check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#"))
    save(list = .x, file = file.path(dat_dir, paste0(.x, ".rda")))
  })
}


#' Copy Mascot \code{.csv} files
#'
#' \code{copy_csv} copies the Mascot outputs of \code{.csv} files to a target
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
	# copy_csv(dat_dir, filelist = c("F003481", "F003485", "F003486", "F003487", "F003488", "F003510"))
	copy_csv(dat_dir, filelist = c("F003590", "F003591", "F003593", "F003594", "F003595", "F003597"))
}


#' Copy Mascot \code{.csv} files
#'
#' @export
cptac_csv_2 <- function(dat_dir) {
  copy_csv(dat_dir, filelist = c("F003529", "F003530", "F003531", "F003532", "F003533", "F003534"))
}


#' Copy MaxQuant \code{.txt} files
#'
#' \code{copy_mq_txt} copies the MaxQuant outputs of \code{.txt} files to a
#' target directory.
#' @export
copy_mq_txt <- function(dat_dir, filelist) {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  data(list = filelist, package = "proteoQDA", envir = environment())

  for (i in seq_along(filelist)) {
    df <- get(filelist[i])
    write.table(df, file.path(dat_dir, paste0(filelist[i], ".txt")), sep = "\t",
                col.names = TRUE, row.names = FALSE)
  }
}


#' Copy MaxQuant \code{.txt} files
#'
#' @export
cptac_mqpep_txt <- function(dat_dir) {
  copy_mq_txt(dat_dir,
              filelist = c("peptides_bi_1", "peptides_bi_2", "peptides_jhu_1",
                           "peptides_jhu_2", "peptides_pnnl_1", "peptides_pnnl_2"))
}


#' Copy MaxQuant \code{.txt} files
#'
#' @export
cptac_mqpep_txt2 <- function(dat_dir) {
  copy_mq_txt(dat_dir,
              filelist <- c("modificationSpecificPeptides_bi_1", "modificationSpecificPeptides_bi_2",
                            "modificationSpecificPeptides_jhu_1", "modificationSpecificPeptides_jhu2",
                            "modificationSpecificPeptides_pnnl_1", "modificationSpecificPeptides_pnnl_2"))
}




#' Copy an \code{expt_smry...} file to \code{dat_dir}
#'
#' @export
cptac_mqpep_expt <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_cptac_gl_mq.xlsx", "expt_smry.xlsx")
}


#' Copy \code{frac_smry.xlsx}
#'
#' @export
cptac_mqpep_frac <- function(dat_dir) {
  copy_frac(dat_dir, "frac_smry_cptac_gl.xlsx", "frac_smry.xlsx")
}


#' Copy protein tables
#'
#' \code{copy_prn_tbl} copies protein tables to a target directory.
#' @export
copy_prn_tbl <- function(dat_dir, from, to) {
  dir.create(file.path(dat_dir, "Protein"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Calls"), recursive = TRUE, showWarnings = FALSE)

  data(list = from, package = "proteoQDA", envir = environment())

  for (i in seq_along(from)) {
    df <- get(from[i])
    write.table(df, to[i], sep = "\t", col.names = TRUE, row.names = FALSE)
  }
}


#' Copy \code{Protein.txt} and \code{Protein_impNA.txt}
#'
#' @export
cptac_prn_1 <- function(dat_dir) {
  copy_prn_tbl(dat_dir,
               from = c("Protein_cptac_1", "Protein_impNA_cptac_1", "normPrn_pars_cptac_1"),
               to = c(file.path(dat_dir, "Protein\\Protein.txt"),
                      file.path(dat_dir, "Protein\\Protein_impNA.txt"),
                      file.path(dat_dir, "Calls\\normPrn.txt")))
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
cptac_expt_3 <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_raneff.xlsx", "expt_smry.xlsx")
}


#' Copy an \code{expt_smry...} file to \code{dat_dir}
#'
#' @export
cptac_mq_expt_1 <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_cptac_gl_mq.xlsx", "expt_smry.xlsx")
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


#' Copy \code{frac_smry.xlsx}
#'
#' @export
cptac_frac_3 <- function(dat_dir) {
  copy_frac(dat_dir, "frac_smry_cptac_gl.xlsx", "frac_smry.xlsx")
}



#' Copy \code{fasta} files
#'
#' \code{copy_fasta} copies \code{fasta} files to a target directory.
#' @export
copy_fasta <- function(fasta_dir = "~\\proteoQ\\dbs\\fasta\\refseq",
                       from = "refseq_hs_2013_07.fasta", to = "refseq_hs_2013_07.fasta") {
  dir.create(file.path(fasta_dir), recursive = TRUE, showWarnings = FALSE)

  filepath <- system.file("extdata", from, package = "proteoQDA")
  filepath <- gsub("/", "\\\\", filepath)
  file.copy(from = filepath, to = file.path(fasta_dir, to))
}


#' Copy \code{refseq_hs_2013_07.fasta}
#'
#' @export
copy_refseq_hs <- function(fasta_dir = "~\\proteoQ\\dbs\\fasta\\refseq") {
  copy_fasta(fasta_dir, "refseq_hs_2013_07.fasta", "refseq_hs_2013_07.fasta")
}


#' Copy \code{refseq_mm_2013_07.fasta}
#'
#' @export
copy_refseq_mm <- function(fasta_dir = "~\\proteoQ\\dbs\\fasta\\refseq") {
  copy_fasta(fasta_dir, "refseq_mm_2013_07.fasta", "refseq_mm_2013_07.fasta")
}




