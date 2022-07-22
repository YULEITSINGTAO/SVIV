# VCFComparison

## Introduction
For structure variation analysis based on whole genome sequencing (WGS), usually we apply different callers to get more comprehensive results (i.e. Multiple of files of a sample). However, we do not have specific tools to compare different results. Therefore,  we start a project to build a R package that focus on dealing with VCF files from different callers in order to get more accurate structure variations and visualize the results. 

This package is still under development, hope we can finish it soon!

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

