# VCFComparison

## Introduction

Whole genome sequencing (WGS) is a powerful tool for detecting structure variations (SV) in genomes. Typically, multiple SV callers are applied to a single patient's WGS data to increase the comprehensiveness of the results. However, there is currently a lack of tools to compare and merge the results from different callers. To address this issue, we have developed an R package called VCFComparison. This package enables users to manipulate SV results from different callers and identify the most confident structure mutations. Furthermore, we propose a projection method for complex translocations, which projects each translocation to a point on Cartesian coordinate system. This is followed by a clustering method to characterize the key translocation information and summarize the translocation mutations for the user. Overall, our approach provides a useful framework for analyzing SVs in WGS data.

## Install Package

```
library(devtools)
install_github("YULEITSINGTAO/VCFComparison")
```

## Sample map

VCFComparison start from reading VCF files by sample map dataframe. Users can define a dataframe in the format as a sample map.  

| Sample ID  | Caller1 |Caller2|
| ------------- | ------------- | ------------- | 
| Sample_1  | Directory  | Directory |
| Sample_2  | Directory  | Directory |

## Read the vcf files into R
Via sample map, function `Read_VCFs()` is applied to read the VCF files into the working environment. 

```
VCF_list <- Read_VCFs(sample_map)
```
In the VCF_list, the VCF files are stored by sample ID

![Data struce of VCF_list](./vignettes/Figures/VCF_Structure_List.png)


## Filter and Extract VCF_list 
```
VCF_list <- Read_VCFs(sample_map)
```
