library(dplyr)

load("intermediate_data/h99_transcript_annotations.Rdata")
load("intermediate_data/H99_to_sac_best_hits.Rdata")
load("intermediate_data/sac_biogrid.Rdata")


sac_physical_network_long <- tidyr::expand_grid(
    target1 = H99_to_sac_best_hits$cnag_id |> unique(),
    target2 = H99_to_sac_best_hits$cnag_id |> unique()) |>
    dplyr::inner_join(
        H99_to_sac_best_hits |>
        dplyr::transmute(
            target1 = cnag_id,
            sac_target1 = sac_target,
            target1_EValue = pmax(H99_to_sac_EValue, sac_to_H99_EValue)),
        by = "target1") |>
    dplyr::inner_join(
        H99_to_sac_best_hits |>
        dplyr::transmute(
            target2 = cnag_id,
            sac_target2 = sac_target,
            target2_EValue = pmax(H99_to_sac_EValue, sac_to_H99_EValue)),
        by = "target2") |>
    dplyr::inner_join(
        sac_biogrid |>
        dplyr::filter(experimental_system_type == "physical") |>
        dplyr::select(
            sac_target1 = feature_name_1,
            sac_target2 = feature_name_2),
        by = c("sac_target1", "sac_target2")) |>
    dplyr::mutate(max_EValue = pmax(target1_EValue, target2_EValue)) |>
    dplyr::arrange(desc(max_EValue)) |>
    dplyr::distinct(target1, target2, .keep_all = TRUE) |>
    dplyr::arrange(target1, target2)


sac_physical_network <- sac_physical_network_long |>
    dplyr::distinct(target1, target2) |>
    as.data.frame() |>
    EGAD::build_binary_network(genes) |>
    EGAD::extend_network()
sac_physical_network[is.na(sac_physical_network)] <- 0


save(
    sac_physical_network_long,
    file = "intermediate_data/sac_physical_network_long.Rdata")
save(
    sac_physical_network,
    file = "intermediate_data/sac_physical_network.Rdata")


sac_genetic_network_long <- tidyr::expand_grid(
    target1 = H99_to_sac_best_hits$cnag_id |> unique(),
    target2 = H99_to_sac_best_hits$cnag_id |> unique()) |>
    dplyr::inner_join(
        H99_to_sac_best_hits |>
        dplyr::transmute(
            target1 = cnag_id,
            sac_target1 = sac_target,
            target1_EValue = pmax(H99_to_sac_EValue, sac_to_H99_EValue)),
        by = "target1") |>
    dplyr::inner_join(
        H99_to_sac_best_hits |>
        dplyr::transmute(
            target2 = cnag_id,
            sac_target2 = sac_target,
            target2_EValue = pmax(H99_to_sac_EValue, sac_to_H99_EValue)),
        by = "target2") |>
    dplyr::inner_join(
        sac_biogrid |>
        dplyr::filter(experimental_system_type == "genetic") |>
        dplyr::select(
            sac_target1 = feature_name_1,
            sac_target2 = feature_name_2),
        by = c("sac_target1", "sac_target2")) |>
    dplyr::mutate(max_EValue = pmax(target1_EValue, target2_EValue)) |>
    dplyr::arrange(desc(max_EValue)) |>
    dplyr::distinct(target1, target2, .keep_all = TRUE) |>
    dplyr::arrange(target1, target2)

sac_genetic_network <- sac_genetic_network_long |>
    dplyr::distinct(target1, target2) |>
    as.data.frame() |>
    EGAD::build_binary_network(genes) |>
    EGAD::extend_network()
sac_genetic_network[is.na(sac_genetic_network)] <- 0

save(
    sac_genetic_network_long,
    file = "intermediate_data/sac_genetic_network_long.Rdata")
save(
    sac_genetic_network,
    file = "intermediate_data/sac_genetic_network.Rdata")
