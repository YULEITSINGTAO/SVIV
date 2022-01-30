sample_map <- list.files("~/shared/CLL_project/phs000879_SV/DELLY_Result", "*.vcf", full.names = TRUE)

sample_map <- data.frame(Patient = tools::file_path_sans_ext(basename(sample_map)), DELLY = sample_map)

VCF_list <- Read_VCFs(sample_map)

## Collect all gt into one dataframe
SV_list <- list()
for (i in names(VCF_list[["DELLY"]])){

    SV_list[[i]] <- as.data.frame(VCF_list[["DELLY"]][[i]]@fix) %>% mutate(patient = i)


}

SV <- do.call("rbind",SV_list)
filtered_SV <- SV %>% filter(FILTER == "PASS")
