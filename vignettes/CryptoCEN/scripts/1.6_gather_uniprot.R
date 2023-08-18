library(CalCEN)
library(httr)

parameters <- CalCEN::load_parameters()

sgd_to_uniprot <- httr::GET(
    "http://www.uniprot.org/uniprot/",
    httr::user_agent(parameters$web_tokens$user_agent),
    query = list(
        query = "organism:559292",
        format = "tab",
        columns = c(
            "id", "entry name", "genes", "database(SGD)") |>
            paste(collapse = ","))) |>
    httr::content() |>
    readr::read_tsv() |>
    dplyr::transmute(
        uniprot_accession = Entry,
        uniprot_entry = `Entry name`,
        sgd_id = `Cross-reference (SGD)` |> stringr::str_replace(";$", ""))

# There are a few proteins with multiple sgd_id values:
#
# r |> dplyr::filter(sgd_id |> stringr::str_detect(";"))
# # A tibble: 5 x 3
#   uniprot_accession uniprot_entry sgd_id               
#   <chr>             <chr>         <chr>                
# 1 P02309            H4_YEAST      S000000213;S000004975
# 2 P02994            EF1A_YEAST    S000006284;S000000322
# 3 P32324            EF2_YEAST     S000005659;S000002793
# 4 P10081            IF4A_YEAST    S000001767;S000003674
# 5 P61830            H3_YEAST      S000000214;S000004976

sgd_to_uniprot |>
    readr::write_tsv(
        paste0("product/sgd_to_uniprot_", CalCEN::date_code(), ".tsv"))

save(sgd_to_uniprot, file = "intermediate_data/sac_to_uniprot.Rdata")
