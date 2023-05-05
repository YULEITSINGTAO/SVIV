#' Project translocation to Cartesian coordinate system
#' @description
#' Visualization of translocations in Cartesian coordinate system.
#'
#' @param translocation_bed dataframe, paired translocation bed dataframe
#' Columns:
#' -Chr_1 and Pos_1: First break point chromosome ID and location.
#' -Chr_2 and Pos_2: Second breal point chromosome ID and location.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' translocation_bed <- data.frame(Chr_1 = paste0("chr", c(7, 8, 9, 9, 13, 15)),
#' Pos_1 = c(8663231, 70602300, 131457166, 33130549, 21746650, 42750778),
#' Chr_2 = paste0("chr", c(5, 1, 2, 6, 11, 2)),
#' Pos_2 = c(37709720, 91853200, 116376668, 43655549, 108585748, 214996194))
#'
#' projectTranslocation(translocation_bed)
#'
projectTranslocation <- function(translocation_bed){

    translocation_bed$Event <- paste0("event_", c(1:nrow(translocation_bed)))
    translocation_bed$Chr_1 <- factor(translocation_bed$Chr_1, levels = paste0("chr", c(1:22, "X", "Y")))
    translocation_bed$Chr_2 <- factor(translocation_bed$Chr_2, levels = paste0("chr", c(1:22, "X", "Y")))

    p <- ggplot2::ggplot(translocation_bed, aes(Pos_1, Pos_2)) + geom_point() +
        facet_grid(vars(Chr_1), vars(Chr_2)) +
        scale_x_continuous(labels  = scales::label_number(scale = 1e-6, suffix = "Mbp", accuracy = 1)) +
        scale_y_continuous(labels  = scales::label_number(scale = 1e-6, suffix = "Mbp", accuracy = 1)) + xlab("Break point 1") + ylab("Break point 2") # for the y axis label


    return(p)
}

#' sortTranslocationBed
#' @description
#' Sort the translocation bed dataframe and output the translocation dataframe following Chr_1 <= Chr_2 and rows sorted from chr1 to chrY as well.
#'
#' @param translocation_bed dataframe, paired translocation bed dataframe
#' Columns:
#' -Chr_1 and Pos_1: First break point chromosome ID and location.
#' -Chr_2 and Pos_2: Second breal point chromosome ID and location.
#'
#' @return sorted translocation bed dataframe
#' @export
#'
#' @examples
#' translocation_bed <- data.frame(Chr_1 = paste0("chr", c(7, 8, 9, 9, 13, 15)),
#' Pos_1 = c(8663231, 70602300, 131457166, 33130549, 21746650, 42750778),
#' Chr_2 = paste0("chr", c(5, 1, 2, 6, 11, 2)),
#' Pos_2 = c(37709720, 91853200, 116376668, 43655549, 108585748, 214996194))
#'
#' sortTranslocationBed(translocation_bed)
#'

sortTranslocationBed <- function(translocation_bed){
    ## convert the class of columns
    translocation_bed$Event <- paste0("event_", c(1:nrow(translocation_bed)))
    translocation_bed$Chr_1 <- factor(translocation_bed$Chr_1, levels = paste0("chr", c(1:22, "X", "Y")))
    translocation_bed$Chr_2 <- factor(translocation_bed$Chr_2, levels = paste0("chr", c(1:22, "X", "Y")))

    ## Sort by column
    for (i in 1:nrow(translocation_bed)) {
        if(as.numeric(translocation_bed[i, "Chr_2"]) <= as.numeric(translocation_bed[i, "Chr_1"])){

            Chr_2 <- translocation_bed[i, "Chr_2"]
            Pos_2 <- translocation_bed[i, "Pos_2"]

            translocation_bed[i, "Chr_2"] <- translocation_bed[i, "Chr_1"]
            translocation_bed[i, "Pos_2"] <- translocation_bed[i, "Pos_1"]

            translocation_bed[i, "Chr_1"] <- Chr_2
            translocation_bed[i, "Pos_1"] <- Pos_2
        }else{}
    }

    ## Sort by row
    translocation_bed <- translocation_bed[order(as.numeric(translocation_bed$Chr_1)),]

    return(translocation_bed)
}

#' clusterOnePairofChromosomes
#' @description
#' Clustering on one pair of chromosomes
#'
#' @param translocation_bed data.frame, paired translocation bed data.frame
#' Columns:
#' -Chr_1 and Pos_1: First break point chromosome ID and location.
#' -Chr_2 and Pos_2: Second breal point chromosome ID and location.
#'
#' @return clustered translocation bed data.frame and clustering plots
# @export
#'
#' @examples
#' translocation_bed <- data.frame(Chr_1 = paste0("chr", c(1, 1, 1, 1, 1, 1)),
#' Pos_1 = c(8663000, 70602500, 131457200, 33130600, 21746650, 42750778),
#' Chr_2 = paste0("chr", c(2, 2, 2, 2, 2, 2)),
#' Pos_2 = c(37709720, 91853200, 116376668, 43655549, 108585748, 214996194))
#' translocation_bed_example <- data.frame()
#' for(i in 1:10){
#'    translocation_bed_new <- translocation_bed %>% dplyr::mutate(Pos_1 = Pos_1 + sample.int(10000, 6)) %>% dplyr::mutate(Pos_2 = Pos_2 + sample.int(10000, 6))
#'    translocation_bed_example <- rbind(translocation_bed_example, translocation_bed_new)
#'}
#' clusterOnePairofChromosomes(translocation_bed_example)
#' @noRd
clusterOnePairofChromosomes <- function(translocation_bed){

    sorted_translocation_bed <- sortTranslocationBed(translocation_bed)
    ## Check if the input bed dataframe has only one pair of chromosomes
    stopifnot(nrow(unique(sorted_translocation_bed %>% select(Chr_1, Chr_2))) == 1)

    pos_df <- translocation_bed %>% dplyr::select(Pos_1, Pos_2)
    scaled_data <- as.matrix(scale(pos_df))
    ## Decide the best k number
    gap_stat <- cluster::clusGap(scaled_data, FUN = kmeans, nstart = 25,
                        K.max = 10, B = 50)
    cluster_k_selection <- factoextra::fviz_gap_stat(gap_stat)
    gap <- gap_stat$Tab[, "gap"]
    se <- gap_stat$Tab[, "SE.sim"]
    decr <- diff(gap) <= 0
    k_number <- cluster::maxSE(gap, se, method = "firstSEmax", SE.factor = 1)
    distance <- factoextra::get_dist(pos_df)
    k_optimal <- stats::kmeans(pos_df, centers = k_number, nstart = 25)
    clustered_plot <- factoextra::fviz_cluster(k_optimal,
                 geom = "point",
                 data = pos_df,
                 ggtheme = ggpubr::theme_classic2()) +
        ggtitle(paste0("k = ", k_number)) + xlab(translocation_bed[1, "Chr_1"]) + ylab(translocation_bed[1, "Chr_2"])+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"))
    clustered_means <- k_optimal$centers %>% as.data.frame()
    clustered_means$Chr_1 <- translocation_bed[1, "Chr_1"]
    clustered_means$Chr_2 <- translocation_bed[1, "Chr_2"]
    clustered_means <- clustered_means[, c("Chr_1", "Pos_1", "Chr_2", "Pos_2")]
    out_figures <- ggpubr::ggarrange(cluster_k_selection, clustered_plot, ncol = 2, nrow = 1)
    out_clustered_means <- clustered_means
    out_list <- list("out_clustered_means" = out_clustered_means, "out_figures" = out_figures)
    return(out_list)
}

#' Project translocation to Cartesian coordinate system and clustering
#'
#' @param translocation_bed dataframe, paired translocation bed dataframe
#' Columns:
#' -Chr_1 and Pos_1: First break point chromosome ID and location.
#' -Chr_2 and Pos_2: Second breal point chromosome ID and location.
#'
#' @param method clustering method, now only support "K-means", will add other methods into the function
#'
#' @param least_translocation_number_of_clustering numeric,
#' when the Chr_1 to Chr_2 translocation number > least_translocation_number_of_clustering, perform clustering,
#' otherwise no need to clustering
#'
#' @return clustered result
#' @export
#'
#' @examples
#' \dontrun{
#' translocation_bed_1 <- data.frame(Chr_1 = paste0("chr", c(7, 8, 9, 9, 13, 15)),
#' Pos_1 = c(8663000, 70602500, 131457200, 33130600, 21746650, 42750778),
#' Chr_2 = paste0("chr", c(5, 1, 2, 6, 11, 2)),
#' Pos_2 = c(37709720, 91853200, 116376668, 43655549, 108585748, 214996194))

#' translocation_bed_2 <- data.frame(Chr_1 = paste0("chr", c(7, 8, 9, 9, 13, 15)),
#'                                   Pos_1 = c(8663000, 70602500, 131457200, 33130600, 21746650, 42750778)/100,
#'                                   Chr_2 = paste0("chr", c(5, 1, 2, 6, 11, 2)),
#'                                   Pos_2 = c(37709720, 91853200, 116376668, 43655549, 108585748, 214996194)/100)
#'
#' translocation_bed_1_2 <- rbind(translocation_bed_1, translocation_bed_2)
#' translocation_bed <- data.frame()

#' for(i in 1:10){
#'    translocation_bed_new <- translocation_bed_1_2 %>% dplyr::mutate(Pos_1 = Pos_1 + sample.int(1000, 12)) %>% dplyr::mutate(Pos_2 = Pos_2 + sample.int(1000, 12))
#'    translocation_bed <- rbind(translocation_bed, translocation_bed_new)
#'}
#'
#'projectTranslocationClustering(translocation_bed)
#'
projectTranslocationClustering <- function(translocation_bed, method = "K-Means", least_translocation_number_of_clustering = 10){
    sorted_translocation_bed <- sortTranslocationBed(translocation_bed)

    count_chr_1_chr_2 <- sorted_translocation_bed %>% dplyr::select(Chr_1, Chr_2) %>%
        dplyr::group_by_all() %>% dplyr::summarise(total_count = n(), .groups = 'drop') %>%
        as.data.frame()

    selected_chr_1_chr_2_pair <- count_chr_1_chr_2 %>% dplyr::filter(total_count >= least_translocation_number_of_clustering) %>%
        dplyr::select(Chr_1, Chr_2)

    out_table_list <- list()
    out_figure_list <- list()

    for (i in 1:nrow(selected_chr_1_chr_2_pair)) {
        input_one_pair_chromosome <- sorted_translocation_bed %>% dplyr::filter(Chr_1 == selected_chr_1_chr_2_pair[i,1], Chr_2 == selected_chr_1_chr_2_pair[i,2])
        cluster_one_pair_of_chromosomes_out <- clusterOnePairofChromosomes(input_one_pair_chromosome)
        out_table_list[[i]] <- cluster_one_pair_of_chromosomes_out[["out_clustered_means"]]
        out_figure_list[[i]] <- cluster_one_pair_of_chromosomes_out[["out_figures"]]

    }

    ## Make table into one
    one_table <- do.call("rbind", out_table_list)

    ## Make figures into one
    one_figure <- do.call(gridExtra::grid.arrange, c(out_figure_list, ncol = 2))

    result_list <- list("merged_trans_table" = one_table, "clustering_plot" = one_figure)
    return(result_list)
}
