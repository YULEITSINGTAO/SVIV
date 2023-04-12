#' SV type extraction
#' This function is used to extract specific SV types from the SV list
#' @param VCF_list
#' @param SV_type
#'
#' @return
#' @export
#'
#' @examples
#' SV_type_List(SVF_list, SV_type)
SV_type_List <- function(VCF_list, SV_type){
    SV_type_list <- list()
    for (Caller in names(VCF_list)) {
        if (SV_type == "DEL"){
            for (j in names(VCF_list[[Caller]])){

                Deletion_table <- VCF_list[[Caller]][[j]] %>% filter(grepl("SVTYPE=DEL",INFO)|grepl("SIMPLE_TYPE=DEL",INFO))
                bed_Deletion_table <- cbind(Deletion_table[,1:2],str_extract(string = Deletion_table[,8], pattern = "(?<=END=)\\d+"))
                names(bed_Deletion_table) <- c("Chr", "Start", "End")
                bed_Deletion_table$Start <- as.numeric(bed_Deletion_table$Start)
                bed_Deletion_table$End <- as.numeric(bed_Deletion_table$End)
                SV_type_list[[Caller]][[j]] <- union_bed_df(bed_Deletion_table)
                }
            }else if (SV_type == "DUP") {
                for (j in names(VCF_list[[Caller]])){

                ##########################################################################
                Duplication_table <- VCF_list[[Caller]][[j]] %>% filter(grepl("SVTYPE=DUP",INFO)|grepl("SIMPLE_TYPE=DUP",INFO))
                bed_Duplication_table <- cbind(Duplication_table[,1:2],str_extract(string = Duplication_table[,8], pattern = "(?<=END=)\\d+"))

                names(bed_Duplication_table) <- c("Chr", "Start", "End")
                bed_Duplication_table$Start <- as.numeric(bed_Duplication_table$Start)
                bed_Duplication_table$End <- as.numeric(bed_Duplication_table$End)
                SV_type_list[[Caller]][[j]] <- union_bed_df(bed_Duplication_table)
                ############################################################################
                }
            }else if (SV_type == "INV"){
                for (j in names(VCF_list[[Caller]])){
                    Inversion_table <- VCF_list[[Caller]][[j]] %>% filter(grepl("SVTYPE=INV",INFO)|grepl("SIMPLE_TYPE=INV",INFO))
                    bed_Inversion_table <- cbind(Inversion_table[,1:2],str_extract(string = Inversion_table[,8], pattern = "(?<=END=)\\d+"))

                    names(bed_Inversion_table) <- c("Chr", "Start", "End")
                    bed_Inversion_table$Start <- as.numeric(bed_Inversion_table$Start)
                    bed_Inversion_table$End <- as.numeric(bed_Inversion_table$End)
                    SV_type_list[[Caller]][[j]] <- union_bed_df(bed_Inversion_table)
                }

            }else if (SV_type == "INS"){
                for (j in names(VCF_list[[Caller]])){

                    Insertion_table <- VCF_list[[Caller]][[j]] %>% filter(grepl("SVTYPE=INS",INFO)|grepl("SIMPLE_TYPE=INS",INFO))
                    bed_Insertion_table <- cbind(Insertion_table[,1:2],str_extract(string = Insertion_table[,8], pattern = "(?<=END=)\\d+"))

                    names(bed_Insertion_table) <- c("Chr", "Start", "End")
                    bed_Insertion_table$Start <- as.numeric(bed_Insertion_table$Start)
                    bed_Insertion_table$End <- as.numeric(bed_Insertion_table$End)
                    SV_type_list[[Caller]][[j]] <- union_bed_df(bed_Insertion_table)
                }

            }else if (SV_type == "TRANS"){
                for (j in names(VCF_list[[Caller]])){
                    Translocation_table <- VCF_list[[Caller]][[j]] %>% filter(grepl("SVTYPE=BND",INFO))
                    chr_end <- as.data.frame(str_extract(string = Translocation_table[,5], pattern = "(?<=\\[).*(?=\\[)|(?<=\\]).*(?=\\])")
                                             %>%str_split_fixed(":", 2))

                    chr_start <- Translocation_table[,1:2]

                    chr_start[,2] <- as.numeric(chr_start[,2])
                    chr_end[,2] <- as.numeric(chr_end[,2])
                    Translocation_bed <- cbind(chr_start, chr_end)
                    colnames(Translocation_bed) <- c("Chr1", "Pos1", "Chr2","Pos2")
                    SV_type_list[[Caller]][[j]] <- Translocation_bed %>% filter(grepl("chr[0-9|X|Y]{1,}$", Chr1) & grepl("chr[0-9|X|Y]{1,}$", Chr2))
                    }
                } else{
                    print("Please input the required SV_type: One of DEL, DUP, INS, INV, TRANS")
                }
        }
    return(SV_type_list)
    }












