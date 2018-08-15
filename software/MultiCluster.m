% Two-Mode HOSVD with semi-nonnegative constraint
% Input:
% Tensor: an order-3 tensor 
% Ncomp: number of components to extract 
% Output:
% output_vector_X: a matrix with Ncomp columns, where each column is an estimated singular vector in the X-mode. 
% output_vector_Y: a matrix with Ncomp columns, where each column is an estimated singular vector in the Y-mode.
% output_vector_Z: a matrix with Ncomp columns, where each column is an estimated singular vector in the Z-mode. All entries in the Z-mode are forced to be non-negative. 
% output_value: a lenghth-Ncomp vector, where each value is an estimated singular value.


function[output_vector_X,output_vector_Y,output_vector_Z,output_value]=MultiCluster(Tensor,Ncomp)
d=size(Tensor);
d1=d(1);
d2=d(2);
d3=d(3);


output_vector_X=zeros(d1,Ncomp);
output_vector_Y=zeros(d2,Ncomp);
output_vector_Z=zeros(d3,Ncomp);
output_value=zeros(1,Ncomp);



for index =1:Ncomp

M=reshape(double(Tensor),[d1*d2,d3]);
[X,~,Z]=svds(M,1);
Matrix.Basis=reshape(X,[d1,d2,1]);
[X,~,Y]=svds(Matrix.Basis,1);
res0=Z'*tenmat(Tensor,3)*reshape(kron(X,Y'),[d1*d2,1]);
res0=res0.data;

%res=cp_als(Tensor,1);
%X=res.U{1}(:,1);
%Y=res.U{2}(:,1);
%Z=res.U{3}(:,1);
%res0=res.lambda;

thresh=0.1;
res=res0+thresh+1;

while((res-res0)>thresh)

res0=res;

M=kr(Y,X);
response=tenmat(Tensor,3)*M;
Z=response.data;
Z_option2=Z;
res=0;
res_option2=0; 
 
if(min(Z)<0)
Z(Z>0)=0;
Z=Z/norm(Z);
res=Z'*tenmat(Tensor,3)*reshape(kron(X,Y'),[d1*d2,1]);
res=res.data;
end
                                   
if(max(Z_option2)>0) 
Z_option2(Z_option2<0)=0;
Z_option2=Z_option2/norm(Z_option2);
res_option2=Z_option2'*tenmat(Tensor,3)*reshape(kron(X,Y'),[d1*d2,1]);
res_option2=res_option2.data;
end 

if(abs(res_option2)>abs(res))
Z=Z_option2;
res=res_option2;
end


M=kr(Z,Y);
X=tenmat(Tensor,1)*M;
X=X.data/norm(X.data);

M=kr(Z,X);
Y=tenmat(Tensor,2)*M;
Y=Y.data/norm(Y.data);

end
                                                                 
output_vector_X(:,index)=X;
output_vector_Y(:,index)=Y;
output_vector_Z(:,index)=Z; 
output_value(index)=res;

Tensor=Tensor-output_value(index)*reshape(kron(output_vector_Z(:,index),kron(output_vector_Y(:,index),output_vector_X(:,index))),[d1,d2,d3]);

end
end
                                         
