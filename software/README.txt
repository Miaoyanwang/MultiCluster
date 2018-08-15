###############################################################################
Software accompaniment to:

Three-way clustering of multi-tissue multi-individual gene expression data using constrained tensor decomposition

Miaoyan Wang, Jonathan Fischer, and Yun S. Song
University of California, Berkeley

###############################################################################

Tensor decomposition with semi-nonnegative constraint:

This software tool takes an order-3 tensor as input, and outputs the triplets of singular vectors and singular values. In the context of multi-tissue, multi-individual gene expression, the input tensor should be in the form of gene-by-individual-by-tissue array. The decomposition algorithm imposes nonnegativity on the tissue-mode (i.e. Z-mode). 
 
The package contains the main function MultiCluster.m and auxiliary function positive.m.  See demo.m for an example of using the main function. 


function[output_vector_X,output_vector_Y,output_vector_Z,output_value]=MultiCluster(T,Ncomp)

###
Input:
T: an order-3 tensor 
Ncomp: number of components to extract 

###
Output:
output_vector_X: a matrix with Ncomp columns, where each column is an estimated singular vector in the X-mode. 
output_vector_Y: a matrix with Ncomp columns, where each column is an estimated singular vector in the Y-mode.
output_vector_Z: a matrix with Ncomp columns, where each column is an estimated singular vector in the Z-mode. All entries in the Z-mode are non-negative. 
output_value: a vector of length Ncomp, where each value is an estimated singular value.
