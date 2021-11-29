


#Kim, H., Jung, KW., Maeng, S. et al. Network-assisted genetic
#dissection of pathogenicity and drug resistance in the opportunistic
#human pathogenic fungus Cryptococcus neoformans. Sci Rep 5, 8767
#(2015). https:// doi.org/10.1038/srep08767
#
#http://www.inetbio.org/cryptonet/

library(tidyverse)
library(CalCEN)

parameters <- CalCEN::load_parameters()

#CryptoNet V1
cmd <- "mkdir -p raw_data/CryptoNetV1"
cat(cmd, "\n")
system(cmd)

# Integrated network
# [CryptoNet: a network by integration of all data-type-specific networks (CN-CC,CN-CX,CN-DC,CN-GN,CN-PG,HS-HT,HS-LC,SC-CC,SC-CX,SC-DC,SC-GT,SC-HT,SC-LC,SC-TS)]
# 5649 genes x 156506 links
cmd <- paste0(
    "wget  -O \"", parameters$source_data$CryptoNetV1$CryptoNetV1_path, "\" ",
    "\"", parameters$source_data$CryptoNetV1$CryptoNetV1_url, "\"")
cat(cmd, "\n")
system(cmd)
CryptoNetV1 <- readr::read_tsv(
    file = parameters$source_data$CryptoNetV1$CryptoNetV1_path,
    col_names = c("target1", "target2", "score"))
save(CryptoNetV1, file = "intermediate_data/CryptoNetV1.Rdata")


# GS: Gold standard data of positive gene functional associations
# 1414 genes x 43524 links
cmd <- paste0(
    "wget -O ", parameters$source_data$CryptoNetV1$CryptoNetV1_GS_path, " ",
    parameters$source_data$CryptoNetV1$CryptoNetV1_GS_url)
cat(cmd, "\n")
system(cmd)
CryptoNetV1_GS <- readr::read_tsv(
    file = parameters$source_data$CryptoNetV1$CryptoNetV1_GS_path,
    col_names = c("target1", "target2"))
save(CryptoNetV1_GS, file = "intermediate_data/CryptoNetV1_GS.Rdata")


# sub networks
for (subnetwork_type in parameters$source_data$CryptoNetV1$CryptoNetV1_subnetwork_types) {
    cat("Gathering CryptoNetV1 subentwork: ", subnetwork_type, " ... \n", sep = "")
    url <- parameters$source_data$CryptoNetV1$CryptoNetV1_subnetwork_url %>%
        stringr::str_replace("SUBNETWORK", subnetwork_type)
    path <- parameters$source_data$CryptoNetV1$CryptoNetV1_subnetwork_path %>%
        stringr::str_replace("SUBNETWORK", subnetwork_type)
    cmd <- paste0("wget -O \"", path, "\" \"", url, "\"")
    cat(cmd, "\n")
    system(cmd)
}


CryptoNetV1_subnetworks <- parameters$source_data$CryptoNetV1$CryptoNetV1_subnetwork_types %>%
    purrr::map(
        .f = function(subnetwork_type) {
            readr::read_tsv(
                file = parameters$source_data$CryptoNetV1$CryptoNetV1_subnetwork_path %>%
                    stringr::str_replace("SUBNETWORK", subnetwork_type),
                col_names = c("target1", "target2", "score"),
                col_types = readr::cols(
                    target1 = readr::col_character(),
                    target2 = readr::col_character(),
                    score = readr::col_double()))

        })
names(CryptoNetV1_subnetworks) <- parameters$source_data$CryptoNetV1$CryptoNetV1_subnetwork_types
save(
    CryptoNetV1_subnetworks,
    file = paste0("intermediate_data/CryptoNetV1_subnetworks.Rdata"))


 
