library(dplyr)
library(stringr)
library(readr)
library(CalCEN)

if (!dir.exists("raw_data/biogrid")) {
    cat("Creating 'raw_data/biogrid'\n")
    dir.create("raw_data/biogrid")
}


biogrid_version <- "4.4.216"


system(paste0("\\
  cd raw_data/biogrid && \\
  wget https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/",
      "BIOGRID-", biogrid_version, "/",
      "BIOGRID-ORGANISM-", biogrid_version, ".tab2.zip && \\
  unzip BIOGRID-ORGANISM-", biogrid_version, ".tab2.zip && \\
  ls | grep -v -e 'Candida_albicans' -e 'Saccharomyces_cerevisiae' | xargs rm")

ca_biogrid <- read_biogrid_tab2(
    fname = paste0(
        "raw_data/biogrid/BIOGRID-ORGANISM-Candida_albicans_SC5314-",
        biogrid_version, ".tab2.txt"),
    taxon = 237561)

sac_biogrid <- read_biogrid_tab2(
    fname = paste0(
        "raw_data/biogrid/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-",
        biogrid_version, ".tab2.txt"),
    taxon = 559292)

save(ca_biogrid, file="intermediate_data/ca_biogrid.Rdata")
save(sac_biogrid, file="intermediate_data/sac_biogrid.Rdata")
