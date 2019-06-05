usethis::create_package("C:\\Results\\R\\proteoQDA") 

# save rda
library(purrr)
dat_dir <- c("C:\\Results\\Mertins_Nat_Protoc_2018") # the WORKING directory
filelist = list.files(path = file.path(dat_dir), pattern = "^F[0-9]{6}\\$")

for (i in seq_along(filelist)) {
	fn_prx <- gsub("\\$", "", filelist[i])
	assign(fn_prx, readLines(file.path(dat_dir, filelist[i])))
	save(list = fn_prx, file = file.path("C:\\Results\\R\\proteoQDA\\data", paste0(fn_prx, ".rda")))
}


# copy .csv files
devtools::document(pkg  = "C:\\Results\\R\\proteoQDA")
# devtools::install("C:\\Results\\R\\proteoQDA")
# library(proteoQDA)
dat_dir <- "C:\\my_directory"

# copy global CPTAC data to `dat_dir`
copy_globalCPTAC(dat_dir)

# copy global `expt_smry.xlsx` to `dat_dir`
copy_expt(dat_dir)

# copy global `frac_smry.xlsx` to `dat_dir`
copy_frac(dat_dir)
