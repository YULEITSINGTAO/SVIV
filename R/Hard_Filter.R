#' Hard Filter the VCF
#' @description Filter the VCF list by filter tag
#' @param VCF_list VCF_list, the list of input vcf files.
#' @param filter_standard character, choose from "PASS", "Precise", or "Both".
#' @return Filtered VCF_list
#' @export
#'
#' @examples
#'\dontrun{
#'Read_VCFs(sample_mapping_df)
#'}
#'
hardFilter <- function(VCF_list, filter_standard){

    Chromosome <- paste0("chr", c(1:22, "X","Y"))
    Filtered_VCF_list <- list()

    if(filter_standard == "PASS"){
        for (i in names(VCF_list)){
            temp_list <- list()
            for (j in names(VCF_list[[i]])){
                temp_list[[j]] <- as.data.frame(VCF_list[[i]][[j]]@fix) %>% filter(FILTER =="PASS", CHROM %in% Chromosome)
            }
            Filtered_VCF_list[[i]] <- temp_list
        }
    }

    if(filter_standard == "Precise"){
        for (i in names(VCF_list)){
            temp_list <- list()
            for (j in names(VCF_list[[i]])){
                temp_list[[j]] <- as.data.frame(VCF_list[[i]][[j]]@fix) %>% filter(grepl("PRECISE", INFO, ignore.case = TRUE), CHROM %in% Chromosome)
            }
            Filtered_VCF_list[[i]] <- temp_list
        }
    }

    if(filter_standard == "Both"){
        for (i in names(VCF_list)){
            temp_list <- list()
            for (j in names(VCF_list[[i]])){
                temp_list[[j]] <- as.data.frame(VCF_list[[i]][[j]]@fix) %>% filter(FILTER =="PASS", grepl("PRECISE", INFO, ignore.case = TRUE), CHROM %in% Chromosome)
            }
            Filtered_VCF_list[[i]] <- temp_list
        }
    }

    return(Filtered_VCF_list)

}
