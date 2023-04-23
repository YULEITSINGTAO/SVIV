#' Hard Filter the VCF
#' @description Filter the VCFList by filter tag
#'
#' @param VCFlist list, the list of input vcf files.
#' @param filterStandard character, choose from "PASS", "Precise", or "Both".
#' @return Filtered VCF_List
#' @export
#'
#' @examples
#'
#' ## Read the VCFlist example
#'
#'
#'
#'
#'
#'
#'


Hard_Filter <- function(VCF_List, Filter_Standard){

    Chromosome <- paste0("chr", c(1:22, "X","Y"))
    Filtered_VCF_List <- list()

    if(Filter_Standard == "PASS"){
        for (i in names(VCF_List)){
            temp_list <- list()
            for (j in names(VCF_List[[i]])){
                temp_list[[j]] <- as.data.frame(VCF_List[[i]][[j]]@fix) %>% filter(FILTER =="PASS", CHROM %in% Chromosome)
            }
            Filtered_VCF_List[[i]] <- temp_list
        }
    }

    if(Filter_Standard == "Precise"){
        for (i in names(VCF_List)){
            temp_list <- list()
            for (j in names(VCF_List[[i]])){
                temp_list[[j]] <- as.data.frame(VCF_List[[i]][[j]]@fix) %>% filter(grepl("PRECISE", INFO, ignore.case = TRUE), CHROM %in% Chromosome)
            }
            Filtered_VCF_List[[i]] <- temp_list
        }
    }

    if(Filter_Standard == "Both"){
        for (i in names(VCF_List)){
            temp_list <- list()
            for (j in names(VCF_List[[i]])){
                temp_list[[j]] <- as.data.frame(VCF_List[[i]][[j]]@fix) %>% filter(FILTER =="PASS", grepl("PRECISE", INFO, ignore.case = TRUE), CHROM %in% Chromosome)
            }
            Filtered_VCF_List[[i]] <- temp_list
        }
    }

    return(Filtered_VCF_List)

}
