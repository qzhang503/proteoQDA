#' Mascot subset rda files
#'
#' @import purrr magrittr
foo_mascot_psmidx <- function () {
  dat_dir <- "~/proteoQ/examples"

  ## global
  filelist <- c("F003590_hdr_rm.csv", "F003591_hdr_rm.csv", "F003593_hdr_rm.csv",
               "F003594_hdr_rm.csv", "F003595_hdr_rm.csv", "F003597_hdr_rm.csv")
  filelist_hdr <- c("F003590_header.txt")

  ## phospho
  # filelist <- c("F003598_hdr_rm.csv", "F003602_hdr_rm.csv", "F003603_hdr_rm.csv",
  #             "F003604_hdr_rm.csv", "F003605_hdr_rm.csv", "F003606_hdr_rm.csv")
  # filelist_hdr <- c("F003598_header.txt")

  ## combined phospho and global
  # filelist <- c("F003607_hdr_rm.csv", "F003608_hdr_rm.csv", "F003609_hdr_rm.csv",
  #               "F003610_hdr_rm.csv", "F003611_hdr_rm.csv", "F003612_hdr_rm.csv")
  # filelist_hdr <- c("F003607_header.txt")

  # data thinning
  purrr::walk(filelist, ~ {
    assign(.x, read.csv(file.path(dat_dir, .x), check.names = FALSE, header = TRUE, comment.char = "#"))

    df <- get(.x)
    r_start <- which(names(df) == "pep_scan_title") + 1
    int_end <- ncol(df)

    TMT_plex <- 10
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C", "I130N", "I130C", "I131")
    df <- df[, ((TMT_plex - 1) * 2 + r_start) : int_end]
    df <- df[, seq(2, 2*TMT_plex, 2)]
    colnames(df) <- col_int

    df[] <- lapply(df, function(.x) {ifelse(.x == -1, NA, .x)})
    df$is_complete <- complete.cases(df)
    df$psm_idx <- 1:nrow(df)

    n_row <- floor(nrow(df)/10)
    df_comp <- df[df$is_complete, ]
    set.seed(123)
    rows <- sample(nrow(df_comp), n_row)
    comp_idx <- df_comp[rows, "psm_idx"]

    n_row2 <- floor(n_row/10)
    df_incomp <- df[!df$is_complete, ]
    set.seed(1234)
    rows2 <- sample(nrow(df_incomp), n_row2)
    incomp_idx <- df_incomp[rows2, "psm_idx"]

    filename <- file.path(dat_dir, gsub("_hdr_rm", "_idx", .x))
    both_idx <- c(comp_idx, incomp_idx)
    both_idx <- sort(both_idx)
    both_idx <- data.frame(idx = both_idx)

    write.table(both_idx, filename, sep = ',', na = "", col.names = TRUE, row.names = FALSE, quote = FALSE)
  })

}

#' Mascot subset rda files
#'
#' @import purrr magrittr
foo_mascot_subset_tenperent_na <- function () {
  dat_dir <- "~/proteoQ/examples"

  ## global
  filelist <- c("F003590", "F003591", "F003593", "F003594", "F003595", "F003597")

  ## phospho
  # filelist <- c("F003598", "F003602", "F003603", "F003604", "F003605", "F003606")

  ## combined phospho and global
  # filelist <- c("F003607", "F003608", "F003609", "F003610", "F003611", "F003612")

  purrr::walk(filelist, ~ {
    assign(.x, readLines(file.path(dat_dir, paste0(.x, ".csv"))))

    data <- get(.x)
    eoh <- grep("prot_hit_num", data)
    hdr <- data[1:eoh]

    rows <- read.csv(file.path(dat_dir, paste0(.x, "_idx.csv")), check.names = FALSE, header = TRUE,
                     comment.char = "#")
    rows <- rows + length(hdr)
    rows <- unlist(rows)

    set.seed(123)
    len <- length(data)
    data_sub <- data[rows]
    data_sub <- append(hdr, data_sub)
    assign(.x, data_sub)
    save(list = .x, file = file.path(dat_dir, paste0(.x, ".rda")), compress = "xz")
  })
}

#' Mascot subset rda files
#'
#' @import purrr magrittr
foo_mascot_subset <- function () {
  dat_dir <- "~/proteoQ/examples"

  ## global
  filelist <- c("F003590", "F003591", "F003593", "F003594", "F003595", "F003597")

  ## phospho
  # filelist <- c("F003598", "F003602", "F003603", "F003604", "F003605", "F003606")

  ## combined phospho and global
  # filelist <- c("F003607", "F003608", "F003609", "F003610", "F003611", "F003612")

  ## F5 and F15
  # F003795: Refseq human and mouse
  # F003797 Refseq human
  # filelist <- c("F003795", "F003797")

  purrr::walk(filelist, ~ {
    assign(.x, readLines(file.path(dat_dir, paste0(.x, ".csv"))))

    data <- get(.x)
    eoh <- grep("prot_hit_num", data)
    hdr <- data[1:eoh]

    set.seed(123)
    len <- length(data)
    data_sub <- data[sample((eoh+1):len, len/10)]
    data_sub <- append(hdr, data_sub)
    assign(.x, data_sub)
    save(list = .x, file = file.path(dat_dir, paste0(.x, ".rda")), compress = "xz")
  })
}


#' MaxQuant subset rda files
#'
#' @import purrr magrittr
foo_mq_subset <- function () {
  dat_dir <- "~/proteoQ/examples"

  ## global
  filelist <- c("msms_bi_1", "msms_jhu_1", "msms_pnnl_1", "msms_bi_2", "msms_jhu_2", "msms_pnnl_2")

  ## phospho
  # filelist <- c("msms_bi_p1", "msms_jhu_p1", "msms_pnnl_p1", "msms_bi_p2", "msms_jhu_p2", "msms_pnnl_p2")

  # data thinning
  purrr::walk(filelist, ~ {
    assign(.x, read.csv(file.path(dat_dir, paste0(.x, ".txt")), sep = "\t",
                        check.names = FALSE, header = TRUE, comment.char = "#"))

    df <- get(.x)
    df <- df[, grepl("Reporter intensity corrected", names(df))]
    df$is_complete <- complete.cases(df)
    df$psm_idx <- 1:nrow(df)

    df_comp <- df[df$is_complete, ]
    set.seed(123)
    rows <- sample(nrow(df_comp), floor(nrow(df_comp)/10))
    comp_idx <- df_comp[rows, "psm_idx"] %>% unclass() %>% unlist()


    df_incomp <- df[!df$is_complete, ]
    set.seed(1234)
    rows2 <- sample(nrow(df_incomp), floor(nrow(df_incomp)/10))
    incomp_idx <- df_incomp[rows2, "psm_idx"] %>% unclass() %>% unlist()

    df <- get(.x)
    df <- df[c(comp_idx, incomp_idx), ]
    df <- df[sample(nrow(df)), ]

    ## poor compression
    # assign(.x, df)
    # save(list = .x, file = file.path(dat_dir, paste0(.x, ".rda")))
    # load(file.path(dat_dir, paste0(.x, ".rda")))

    write.table(df, file.path(dat_dir, paste0(.x, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE)
  })

  ## rda
  # better compression
  purrr::walk(filelist, ~ {
   assign(.x, read.csv(file.path(dat_dir, paste0(.x, ".txt")), check.names = FALSE,
                       header = TRUE, sep = "\t", comment.char = "#"))
   save(list = .x, file = file.path(dat_dir, paste0(.x, ".rda")))
  })
}


#' SpectrumMill subset rda files
#'
#' @import purrr magrittr
foo_sm_subset <- function () {
  dat_dir <- "~/proteoQ/examples"

  ## global
  filelist <- c("PSMexport_bi_1", "PSMexport_bi_2", "PSMexport_jhu_1", "PSMexport_jhu_2", "PSMexport_pnnl_1", "PSMexport_pnnl_2")

  ## phospho
  # filelist <- c("PSMexport_bi_p1", "PSMexport_bi_p2", "PSMexport_jhu_p1", "PSMexport_jhu_p2", "PSMexport_pnnl_p1", "PSMexport_pnnl_p2")

  # data thinning
  purrr::walk(filelist, ~ {
    assign(.x, readr::read_delim(file.path(dat_dir, paste0(.x, ".ssv")), delim = ";"))

    df <- get(.x)
    df <- df[, grepl("^TMT_[0-9]{3}[NC]{0,1}$", names(df))]
    df$is_complete <- complete.cases(df)
    df$psm_idx <- 1:nrow(df)

    df_comp <- df[df$is_complete, ]
    set.seed(123)
    rows <- sample(nrow(df_comp), floor(nrow(df_comp)/10))
    comp_idx <- df_comp[rows, "psm_idx"] %>% unclass() %>% unlist()


    df_incomp <- df[!df$is_complete, ]
    set.seed(1234)
    rows2 <- sample(nrow(df_incomp), floor(nrow(df_incomp)/10))
    incomp_idx <- df_incomp[rows2, "psm_idx"] %>% unclass() %>% unlist()

    df <- get(.x)
    df <- df[c(comp_idx, incomp_idx), ]
    df <- df[sample(nrow(df)), ]

    assign(.x, df)
    save(list = .x, file = file.path(dat_dir, paste0(.x, ".rda")))
    # load(file.path(dat_dir, paste0(.x, ".rda")))
    # write.table(df, file.path(dat_dir, paste0(.x, ".ssv")), sep = ";", col.names = TRUE, row.names = FALSE)
  })

  ## rda
  # purrr::walk(filelist, ~ {
  #  assign(.x, readr::read_delim(file.path(dat_dir, paste0(.x, ".ssv")), delim = ";"))
  #  save(list = .x, file = file.path(dat_dir, paste0(.x, ".rda")))
  # })
}


#' Mascot subset rda files
#'
#' @import purrr magrittr
foo_mascot_fullset <- function () {
  dat_dir <- "~/proteoQ/examples"

  ## global
  # filelist <- c("F003590", "F003591", "F003593", "F003594", "F003595", "F003597")

  ## phospho
  filelist <- c("F003598", "F003602", "F003603", "F003604", "F003605", "F003606")

  ## combined phospho and global
  # filelist <- c("F003607", "F003608", "F003609", "F003610", "F003611", "F003612")

  purrr::walk(filelist, ~ {
    assign(.x, readLines(file.path(dat_dir, paste0(.x, ".csv"))))
    save(list = .x, file = file.path(dat_dir, paste0(.x, ".rda")))
  })

  # load(file.path(dat_dir, paste0(.x, ".rda")))
}


#' Copy Mascot \code{.csv} files
#'
#' \code{copy_csv} copies the Mascot outputs of \code{.csv} files to a target
#' directory.
#'
#' @import rlang purrr dplyr magrittr
#' @param dat_dir A character string to the working directory. The default is to
#'   match the value under the global environment.
#' @param filelist A list of files
copy_csv <- function(dat_dir, filelist) {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  data(list = filelist, package = "proteoQDA", envir = environment())

  for (i in seq_along(filelist)) {
    fileConn <- file(file.path(dat_dir, paste0(filelist[i], ".csv")))
    df <- get(filelist[i])

    for (j in seq_along(df)) df[j] <- gsub("^\\s+", "", df[j])

    writeLines(df, fileConn)
    close(fileConn)
  }
}


#' Copy Mascot global \code{.csv}
#'
#' @inheritParams copy_csv
#' @export
copy_global_mascot <- function(dat_dir) {
	copy_csv(dat_dir, filelist = c("F003590", "F003591", "F003593", "F003594", "F003595", "F003597"))
}


#' Copy Mascot phospho \code{.csv}
#'
#' @inheritParams copy_csv
#' @export
copy_phospho_mascot <- function(dat_dir) {
  copy_csv(dat_dir, filelist = c("F003598", "F003602", "F003603", "F003604", "F003605", "F003606"))
}


#' Copy MaxQuant \code{.txt} files
#'
#' \code{copy_txt} copies the MaxQuant \code{msms.txt} files to a target
#' directory.
#' @inheritParams copy_csv
copy_txt <- function(dat_dir, filelist) {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  data(list = filelist, package = "proteoQDA", envir = environment())

  for (i in seq_along(filelist)) {
    df <- get(filelist[i])
    write.table(df, file.path(dat_dir, paste0(filelist[i], ".txt")), sep = "\t",
                col.names = TRUE, row.names = FALSE)
  }
}


#' Copy MaxQuant global \code{.txt}
#'
#' @inheritParams copy_csv
#' @export
copy_global_maxquant <- function(dat_dir) {
  copy_txt(dat_dir, filelist = c("msms_bi_1", "msms_jhu_1", "msms_pnnl_1", "msms_bi_2", "msms_jhu_2", "msms_pnnl_2"))
}

#' Copy MaxQuant phospho \code{.txt}
#'
#' @inheritParams copy_csv
#' @export
copy_phospho_maxquant <- function(dat_dir) {
  copy_txt(dat_dir, filelist = c("msms_bi_p1", "msms_jhu_p1", "msms_pnnl_p1", "msms_bi_p2", "msms_jhu_p2", "msms_pnnl_p2"))
}


#' Copy Spectrum Mill \code{.ssv} files
#'
#' \code{copy_ssv} copies the Spectrum Mill outputs of \code{.ssv} files to a
#' target directory.
#' @inheritParams copy_csv
copy_ssv <- function(dat_dir, filelist) {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  data(list = filelist, package = "proteoQDA", envir = environment())

  for (i in seq_along(filelist)) {
    df <- get(filelist[i])
    write.table(df, file.path(dat_dir, paste0(filelist[i], ".ssv")), sep = ";",
                col.names = TRUE, row.names = FALSE)
  }
}


#' Copy Spectrum Mill global \code{.ssv}
#'
#' @inheritParams copy_csv
#' @export
copy_global_sm <- function(dat_dir) {
  copy_ssv(dat_dir, filelist = c("PSMexport_bi_1", "PSMexport_bi_2", "PSMexport_jhu_1",
                                 "PSMexport_jhu_2", "PSMexport_pnnl_1", "PSMexport_pnnl_2"))
}

#' Copy Spectrum Mill phospho \code{.ssv}
#'
#' @inheritParams copy_csv
#' @export
copy_phospho_sm <- function(dat_dir) {
  copy_ssv(dat_dir, filelist = c("PSMexport_bi_p1", "PSMexport_bi_p2", "PSMexport_jhu_p1",
                                 "PSMexport_jhu_p2", "PSMexport_pnnl_p1", "PSMexport_pnnl_p2"))
}


#' Copy an \code{expt_smry[...].xlsx} file to \code{dat_dir}
#'
#' \code{copy_expt} copies a system file of \code{expt_smry[...].xlsx} to the
#' target directory specified by \code{dat_dir}.
#' @inheritParams copy_csv
#' @param from The source \code{.xlsx} file name.
#' @param to The target \code{.xlsx} file name.
copy_expt <- function(dat_dir, from = "expt_smry.xlsx", to = "expt_smry.xlsx") {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  filepath <- system.file("extdata", from, package = "proteoQDA")
  file.copy(from = filepath, to = file.path(dat_dir, to))
}

#' Copy a metadata file \code{frac_smry.xlsx}
#'
#' \code{copy_frac} copies the \code{frac_smry[...].xlsx} to a target directory
#' specified by \code{dat_dir}.
#' @inheritParams copy_expt
copy_frac <- function(dat_dir, from = "frac_smry.xlsx", to = "frac_smry.xlsx") {
  dir.create(file.path(dat_dir), recursive = TRUE, showWarnings = FALSE)

  filepath <- system.file("extdata", from, package = "proteoQDA")
  file.copy(from = filepath, to = file.path(dat_dir, to))
}

#' Copy a metadata file \code{expt_smry.xlsx} to \code{dat_dir}
#'
#' @inheritParams copy_expt
#' @export
copy_global_exptsmry <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_cptac_gl.xlsx", "expt_smry.xlsx")
}

#' Copy a metadata file \code{frac_smry.xlsx} to \code{dat_dir}
#'
#' @inheritParams copy_expt
#' @export
copy_global_fracsmry <- function(dat_dir) {
  copy_frac(dat_dir, "frac_smry_cptac_gl.xlsx", "frac_smry.xlsx")
}

#' Copy a metadata file \code{expt_smry.xlsx} file to \code{dat_dir}
#'
#' @inheritParams copy_expt
#' @export
copy_cmbn_exptsmry <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_cptac_cmbn.xlsx", "expt_smry.xlsx")
}

#' Copy a metadata file \code{frac_smry.xlsx} to \code{dat_dir}
#'
#' @inheritParams copy_expt
#' @export
copy_cmbn_fracsmry <- function(dat_dir) {
  copy_frac(dat_dir, "frac_smry_cptac_cmbn.xlsx", "frac_smry.xlsx")
}

#' Copy a metadata file \code{expt_smry.xlsx} to \code{dat_dir}
#'
#' @inheritParams copy_expt
#' @export
copy_w2ref_exptsmry <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_ref_w2.xlsx", "expt_smry_ref_w2.xlsx")
}


#' Copy a metadata file \code{expt_smry.xlsx} to \code{dat_dir}
#'
#' @inheritParams copy_expt
#' @export
copy_w2w16ref_exptsmry <- function(dat_dir) {
  copy_expt(dat_dir, "expt_smry_ref_w2_w16.xlsx", "expt_smry_ref_w2_w16.xlsx")
}


#' Save \code{fasta} files to \code{.rda}
#'
#' @inheritParams copy_fasta
#' @inheritParams copy_expt
#' @examples \donttest{
#' foo_fasta_rda("~\\proteoQ\\dbs\\fasta\\uniprot", "uniprot_hs_2014_07.fasta")
#' foo_fasta_rda("~\\proteoQ\\dbs\\fasta\\refseq", "refseq_hs_2013_07.fasta")
#' }
#'
#' @import purrr magrittr
foo_fasta_rda <- function(db_path = "~\\proteoQ\\dbs\\fasta\\refseq", from = "refseq_hs_2013_07.fasta") {
  filepath <- file.path(db_path, from)
  nm <- gsub("\\.fasta", "", from)
  suppressWarnings(assign(nm, readLines(filepath)))
  save(list = nm, file = file.path(db_path, paste0(nm, ".rda")), compress = "xz")
}


#' Copy \code{fasta} files
#'
#' \code{copy_fasta} copies \code{fasta} files to a target directory.
#' @param db_path Character string; the local path for database(s).
#' @param file Character string; the file name of a target \code{.rda}.
copy_fasta <- function(db_path = "~\\proteoQ\\dbs\\fasta\\uniprot",
                       file = "uniprot_hs_2014_07.rda") {
  dir.create(file.path(db_path), recursive = TRUE, showWarnings = FALSE)
  filepath <- system.file("extdata", file, package = "proteoQDA")
  load(filepath)
  out_nm <- gsub("\\.rda", "", file)
  write(get(out_nm), file.path(db_path, paste0(out_nm, ".fasta")))
}


#' Copy \code{uniprot_hs_2014_07.fasta} to \code{db_path}
#'
#' @inheritParams copy_fasta
#' @export
copy_uniprot_hs <- function(db_path = "~\\proteoQ\\dbs\\fasta\\uniprot") {
  copy_fasta(db_path, "uniprot_hs_2014_07.rda")
}


#' Copy \code{uniprot_hs_2014_07.fasta} to \code{db_path}
#'
#' @inheritParams copy_fasta
#' @export
copy_uniprot_mm <- function(db_path = "~\\proteoQ\\dbs\\fasta\\uniprot") {
  copy_fasta(db_path, "uniprot_mm_2014_07.rda")
}


#' Copy \code{refseq_hs_2013_07.fasta} to \code{db_path}
#'
#' @inheritParams copy_fasta
#' @export
copy_refseq_hs <- function(db_path = "~\\proteoQ\\dbs\\fasta\\refseq") {
  copy_fasta(db_path, "refseq_hs_2013_07.rda")
}


#' Copy \code{refseq_mm_2013_07.fasta} to \code{db_path}
#'
#' @inheritParams copy_fasta
#' @export
copy_refseq_mm <- function(db_path = "~\\proteoQ\\dbs\\fasta\\refseq") {
  copy_fasta(db_path, "refseq_mm_2013_07.rda")
}




#' Mascot subset rda files
#'
#' @import purrr magrittr
foo_mascot_subset_not_working <- function () {
  dat_dir <- file.path("~", "proteoQ", "examples")

  ## global
  filelist <- c("F003590_hdr_rm.csv", "F003591_hdr_rm.csv", "F003593_hdr_rm.csv",
                "F003594_hdr_rm.csv", "F003595_hdr_rm.csv", "F003597_hdr_rm.csv")
  filelist_hdr <- c("F003590_header.txt")

  ## phospho
  # filelist <- c("F003598_hdr_rm.csv", "F003602_hdr_rm.csv", "F003603_hdr_rm.csv",
  #              "F003604_hdr_rm.csv", "F003605_hdr_rm.csv", "F003606_hdr_rm.csv")
  # filelist_hdr <- c("F003598_header.txt")

  ## combined phospho and global
  # filelist <- c("F003607_hdr_rm.csv", "F003608_hdr_rm.csv", "F003609_hdr_rm.csv",
  #               "F003610_hdr_rm.csv", "F003611_hdr_rm.csv", "F003612_hdr_rm.csv")
  # filelist_hdr <- c("F003607_header.txt")

  # data thinning
  purrr::walk(filelist, ~ {
    assign(.x, read.csv(file.path(dat_dir, .x), check.names = FALSE, header = TRUE, comment.char = "#"))

    df <- get(.x)
    r_start <- which(names(df) == "pep_scan_title") + 1
    int_end <- ncol(df)

    TMT_plex <- 10
    col_int <- c("I126", "I127N", "I127C", "I128N", "I128C", "I129N", "I129C", "I130N", "I130C", "I131")
    df <- df[, ((TMT_plex - 1) * 2 + r_start) : int_end]
    df <- df[, seq(2, 2*TMT_plex, 2)]
    colnames(df) <- col_int

    df[] <- lapply(df, function(.x) {ifelse(.x == -1, NA, .x)})
    df$is_complete <- complete.cases(df)
    df$psm_idx <- 1:nrow(df)

    n_row <- floor(nrow(df)/10)
    df_comp <- df[df$is_complete, ]
    set.seed(123)
    rows <- sample(nrow(df_comp), n_row)
    comp_idx <- df_comp[rows, "psm_idx"]

    n_row2 <- floor(n_row/10)
    df_incomp <- df[!df$is_complete, ]
    set.seed(1234)
    rows2 <- sample(nrow(df_incomp), n_row2)
    incomp_idx <- df_incomp[rows2, "psm_idx"]

    df_sub <- get(.x)
    df_sub <- df_sub[c(comp_idx, incomp_idx), ]
    df_sub <- df_sub[sample(nrow(df_sub)), ]


    hdr_mascot <- readLines(file.path(dat_dir, filelist_hdr))
    hdr_mascot <- gsub("\"", "", hdr_mascot, fixed = TRUE)
    hdr_mascot <- paste0(hdr_mascot, "\n")

    hdr_df <- paste0(names(df_sub), collapse = ",")
    hdr_df <- paste0(hdr_df, "\n")
    header <- append(hdr_mascot, hdr_df)

    filename <- file.path(dat_dir, gsub("_hdr_rm", "", .x))
    cat(header, file = filename)
    write.table(df_sub, filename, append = TRUE, sep = ',', na = "", col.names = FALSE, row.names = FALSE, quote = FALSE)
  })

  # rda
  purrr::walk(gsub("_hdr_rm\\.csv", "", filelist), ~ {
    assign(.x, readLines(file.path(dat_dir, paste0(.x, ".csv"))))
    save(list = .x, file = file.path(dat_dir, paste0(.x, ".rda")))
  })
}


