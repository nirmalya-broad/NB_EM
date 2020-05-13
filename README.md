## Objective

This project is part of the UMI normalization described in the README of [UmiNormalize](https://github.com/nirmalya-broad/UMINormalize). The basic algorithm of UMI normalization assumes that the reads from bacteria generated using [scDual-seq](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1340-x) represents underlying structure of bacterial transcripts. Since a bacterial transcript can span multiple genes, reads generated from that transcript can span multiple genes. Inference of the boundary of the original transcript is not straight forward and according to our knowledge, there is no available method to solve this problem. The first version of the UMI normalization collapsed those reads by making a relatively simple assumption:

Reads from two different transcripts with the same attached UMI are separated by at least a predefined number of bases.

The [first version of UMI](https://www.nature.com/articles/s41598-019-55633-6), we used a fixed value of 500 for this gap between reads to filter out the transcripts, which we estimated using basic histogram analysis of the read distribution. It provided reasonably accurate UMI normalization as can be seen from the following figure where are 

![alt text](https://github.com/nirmalya-broad/NB_EM/blob/master/assets/figures/UmiNorm1.png "Umi Normalization")



