# VCFComparison

## Introduction

Whole genome sequencing (WGS) is a powerful tool for detecting structure variations (SV) in genomes. Typically, multiple SV callers are applied to a single patient's WGS data to increase the comprehensiveness of the results. However, there is currently a lack of tools to compare and merge the results from different callers. To address this issue, we have developed an R package called VCFComparison. This package enables users to manipulate SV results from different callers and identify the most confident structure mutations. Furthermore, we propose a projection method for complex translocations, which projects each translocation to a point on Cartesian coordinate system. This is followed by a clustering method to characterize the key translocation information and summarize the translocation mutations for the user. Overall, our approach provides a useful framework for analyzing SVs in WGS data.

## Install Package

```
library(devtools)
install_github("YULEITSINGTAO/VCFComparison")
```

## Sample map

| Sample ID  | Caller1 |Caller2|
| ------------- | ------------- | ------------- |
| Sample_1  | Directory  | Directory |
| Sample_2  | Directory  | Directory |

## Progress
### Data Process
- [x] 1. Import VCFs
- [ ] 2. Chop intervals (need to improve)
- [x] 3. Define operators to combine the results from different callers
- [ ] 4. Collect SV features 
- [ ] 5. Cluster translocation SVs
### Visualization
- [x] 1. x-y SV overview all chromosomes
- [x] 2. x-y SV overview specific chromosomes
- [ ] 3. Circos plot of SVs histogram

