library(tidyverse)
library(CalCEN)

parameters <- CalCEN::load_parameters()



cmd <- paste0(
    "wget -O", parameters$source_data$EBIComplexes$SacComplexes_path, " ",
    file = parameters$source_data$EBIComplexes$SacComplexes_url)
cat(cmd, "\n", sep = "")
system(cmd)

sac_complexes <- readr::read_tsv(
    file = parameters$source_data$EBIComplexes$SacComplexes_path)

save(sac_complexes, file = "intermediate_data/sac_complexes.Rdata")
