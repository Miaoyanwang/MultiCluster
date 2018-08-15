% creat an order-3 tensor
d1=10;
d2=20;
d3=100;

v1 = rand(1,d1);
v2 = rand(1,d2);
v3 = rand(1,d3);

v1=v1/norm(v1);
v2=v2/norm(v2);
v3=v3/norm(v3);

noise=0.001;
Ncomp=2;
%
Tensor=reshape(kron(kron(v3,v2),v1),[d1,d2,d3])+noise*reshape(rand(1,d1*d2,d3),[d1,d2,d3]);  

% perform tensor decomposition with semi-nonnegative constraint in the Z-mode

[output_vector_X,output_vector_Y,output_vector_Z,output_value]=MultiCluster(Tensor,Ncomp);

% output
output_vector_X(:,1)
output_vector_Y(:,1)
output_vector_Z(:,1)
output_value
