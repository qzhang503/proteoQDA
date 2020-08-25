foo_subcellular <- function () {
  dat_dir <- "~/proteoQ/examples"

  fn <- "uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review--.tab"
  sp <- "hs"

  run_scripts <- FALSE
  if (run_scripts) {
    fn <- "uniprot-filtered-organism__Mus+musculus+(Mouse)+[10090]_+AND+revie--.tab"
    sp <- "mm"

    fn <- "uniprot-filtered-organism__Rattus+norvegicus+(Rat)+[10116]_+AND+re--.tab"
    sp <- "rn"
  }

  scc <- readr::read_tsv(file.path(dat_dir, fn)) %>%
    dplyr::rename(uniprot_acc = Entry, scc = "Subcellular location [CC]") %>%
    dplyr::select(uniprot_acc, scc) %>%
    dplyr::filter(!is.na(scc))

  scc <- scc %>%
    dplyr::mutate(scc = gsub(";.*\\.+?", ".", scc)) %>%
    dplyr::mutate(scc = gsub("\\s*\\{+?.*\\}+?", "", scc)) %>%
    dplyr::mutate(scc = gsub("\\s*Note\\=.*$", "", scc)) %>%
    dplyr::mutate(scc = gsub("\\s*\\[+?.*\\]+?\\:", "", scc)) %>%
    dplyr::mutate(scc = gsub("SUBCELLULAR LOCATION: ", "", scc)) %>%
    dplyr::mutate(scc = gsub("\\.\\s*$", "", scc))

  scc <- as.character(scc$scc) %>%
    strsplit("\\.") %>%
    plyr::ldply(rbind) %>%
    purrr::map(~ gsub("^.*\\,\\s+(.*)", "\\1", .x)) %>%
    dplyr::bind_cols() %>%
    `colnames<-`(paste0("scc_", names(.))) %>%
    dplyr::bind_cols(scc, .) %>%
    dplyr::select(-scc)

  run_scripts <- FALSE
  if (run_scripts) {
    scc <- scc %>%
      tidyr::gather("id", "scc", -uniprot_acc) %>%
      dplyr::select(-id) %>%
      dplyr::filter(!is.na(scc))
  }

  fn <- paste0("scc_", sp)
  assign(fn, scc)

  save(list = fn, file = file.path(dat_dir, paste0(fn, ".rda")), compress = "xz")
}

