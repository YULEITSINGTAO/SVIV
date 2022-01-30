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

    ## check the if all files exist in the sample map
    print("Checking if all files exist")
    files_status <- 0
    for (n in colnames(sample_mapping_df)[-1]) {

        exist_vector <- file.exists(sample_mapping_df[, get("n")])

        if (sum(exist_vector) == nrow(sample_mapping_df)){
            print(paste("All files from caller:", get("n"), "exist"))
        }else{
            print(paste("File:", sample_mapping_df[,get("n")][!exist_vector], "from caller:", n, "is missing"))
            files_status <- files_status + 1
        }

    }
    if(files_status > 0){stop("Please check all files are exist")}


    VCF_list <- list()
    for (i in 2:ncol(sample_mapping_df)){
        files_to_read <- sample_mapping_df[,i]
        VCF_list[[colnames(sample_mapping_df)[i]]] <- lapply(files_to_read, read.vcfR)
        names(VCF_list[[colnames(sample_mapping_df)[i]]]) <- sample_mapping_df[,1]

    }
    return(VCF_list)
    print("Reading Finished")
}
