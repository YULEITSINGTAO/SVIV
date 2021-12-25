### The function is used to hard filter the raw VCF_list ###

### This function is used to filter the VCF_List

#' Hard_Filter the VCF_List
#'
#' @param VCF_List
#' @param Filter_Standard
#'
#' @return
#' @export
#'
#' @examples
Hard_Filter <- function(VCF_List, Filter_Standard){

    Chromosome <- paste0("chr", c(1:22, "X","Y"))


    if(Filter_Standard == "PASS"){
        for (i in 1:VCF_list){


        }
        for (i in 1:length(DELLY_SV_list)){
            DELLY_VCF_list[[i]] <- as.data.frame(DELLY_SV_list[[i]]@fix)
        }


    }

    if(Filter_Standard == "Precise"){


    }

    if(Filter_Standard == "Both"){


    }

    return ()

}
