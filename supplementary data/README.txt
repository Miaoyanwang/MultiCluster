###############################################################################
Supplementary data:

Three-way clustering of multi-tissue multi-individual gene expression data using constrained tensor decomposition

Miaoyan Wang, Jonathan Fischer, and Yun S. Song
University of California, Berkeley
###############################################################################

The supplementary data provide the results from the global and local analyses based on our tensor decomposition method. Each file is named as tensor_*_module_top10.txt, where * specifies the tissue group names. The tensor_global_module_top10.txt provides the results for global analysis of all 44 somatic tissues. 

###
For each tissue group, we run MultiCluster to obtain the top 10 expression modules (i.e., rank-1 tensor components). The number of genes in each module is decided based on a permutation procedure described in Materials and Methods. For each gene in the identified expression module, we apply tensor projection to assess the covariate (i.e. age-, sex-, and race-) effects. All covariates are pre-normalized to have mean 0 and sample s.e. 1.  

###
Each file has the following format:

ncomp: the index of the expression module. The modules are ranked by their singular values in the tensor decomposition. 

label: the sign of gene loading. "pos_end" indicates top gene with positive loading, whereas "neg_end" indicates top genes with negative loading. 

gene_loading_ranked: the value of gene loading in the corresponding tensor component.

gene_name: the gene names. The gene annotation is based on GENCODE version 19. 

age_effect: the estimated age effect obtained using tensor projection method.

p_for_age: uncorrected p-value for the age effect obtained using tensor projection method.

sex_effect: the estimated sex effect obtained using tensor projection method.

p_for_sex: uncorrected p-value for the sex effect obtained using tensor projection method.

race_effect: the estimated race effect obtained using tensor projection method.

p_for_race: uncorrected p-value for the race effect obtained using tensor projection method.



