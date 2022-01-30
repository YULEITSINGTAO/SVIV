## Stats of SVs

SV_stats <- function(SV_list, Type = "ALL", Patient = "ALL", Caller = "ALL", Figure = True){

    if(Type = "ALL"){

        Type = c("DEL", "DUP", "INS", "INV", "TRANS")

    }
    if (Patient = "ALL"){

        Patient = names(VCF_list[[1]])


    }

    if(Caller = "ALL"){

        Caller = names(VCF_list)
    }



}
