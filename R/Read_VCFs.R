#### This is the function to read the VCF files of a cohort into R workspace ##
#' Read_VCFs
#'
#' @param sample_mapping_df
#'
#' @return VCF_list
#' @export
#'
#' @examples
#'
#' Read_VCFs(sample_mapping_df)
Read_VCFs <- function(sample_mapping_df){
    VCF_list <- list()
    for (i in 2:ncol(sample_mapping_df)){
        files_to_read <- sample_mapping_df[,i]
        VCF_list[[colnames(sample_mapping_df)[i]]] <- lapply(files_to_read, read.vcfR)
        names(VCF_list[[colnames(sample_mapping_df)[i]]]) <- sample_mapping_df[,1]

    }
    return(VCF_list)
}
