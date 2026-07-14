dest_dir <- "data/HPS_PUFs"
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

weeks <- 1:12

urls <- sprintf(
  "https://www2.census.gov/programs-surveys/demo/datasets/hhp/2020/wk%d/HPS_Week%02d_PUF_CSV.zip",
  weeks,
  weeks
)

for (url in urls) {
  zip_file <- file.path(dest_dir, basename(url))

  download.file(url, zip_file, mode = "wb")

  extract_dir <- file.path(
    dest_dir,
    tools::file_path_sans_ext(basename(zip_file))
  )

  dir.create(extract_dir, showWarnings = FALSE)

  unzip(zip_file, exdir = extract_dir)
}
