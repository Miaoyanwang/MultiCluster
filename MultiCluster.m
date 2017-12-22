% Two-Mode HOSVD with semi-nonnegative constraint
% Input:
% T: an order-3 tensor 
% Ncomp: number of components to extract 
% Output:
% output_vector_X: a matrix with Ncomp columns, where each column is an estimated singular vector in the X-mode. 
% output_vector_Y: a matrix with Ncomp columns, where each column is an estimated singular vector in the Y-mode.
% output_vector_Z: a matrix with Ncomp columns, where each column is an estimated singular vector in the Z-mode. All entries in the Z-mode are forced to be non-negative. 
% output_value: a lenghth-Ncomp vector, where each value is an estimated singular value.

function[output_vector_X,output_vector_Y,output_vector_Z,output_value]=MultiCluster(T,Ncomp)

d=size(T);
d1=d(1);
d2=d(2);
d3=d(3);
 output_vector_X=zeros(d1,Ncomp);
 output_vector_Y=zeros(d2,Ncomp);
 output_vector_Z=zeros(d3,Ncomp);
 output_value=zeros(1,Ncomp);

% Two-Mode HOSVD

for index =1:Ncomp

%option 1
 M=reshape(double(T),[d1*d2,d3]);
 [Lspace,~,~]=svds(M,1);
 Matrix.Basis=reshape(Lspace,[d1,d2,1]);
 [Lvector,~,Rvector]=svds(Matrix.Basis,1);
output_vector_X(:,index)=Lvector;
output_vector_Y(:,index)=Rvector;
[B,I]=sort(abs(Rvector));
vector=reshape(kron(Lvector,Rvector'),[1,d1*d2])*M;
                    
[Lvector,Rvector,vector]=positive(double(T),Lvector,Rvector,vector);
                    
output_vector_X(:,index)=Lvector;
output_vector_Y(:,index)=Rvector;
output_vector_Z(:,index)=vector;                  
                    
output_value(index)=reshape(kron(Lvector,Rvector'),[1,d1*d2])*M*vector';

                    

%option 2
M=reshape(double(T),[d1,d2*d3]);
[~,~,Rspace]=svds(M,1);
Matrix.Basis=reshape(Rspace,[d2,d3,1]);
[Lvector,~,Rvector]=svds(Matrix.Basis,1);
vector=M*reshape(kron(Lvector,Rvector'),[d2*d3,1]);
                      
[vector,Lvector,Rvector]=positive(double(T),vector,Lvector,Rvector');
output_value_new=vector'*M*reshape(kron(Lvector,Rvector),[d2*d3,1]);


if  output_value_new>output_value(index)
output_value(index)= output_value_new;
output_vector_X(:,index)=vector;
output_vector_Y(:,index)=Lvector;
output_vector_Z(:,index)=Rvector;
end

%option 3
T_perm=permute(T,[2 1 3]);
M=reshape(double(T_perm),[d2,d1*d3]);
[~,~,Rspace]=svds(M,1);
Matrix.Basis=reshape(Rspace,[d1,d3,1]);
[Lvector,~,Rvector]=svds(Matrix.Basis,1);
vector=M*reshape(kron(Lvector,Rvector'),[d1*d3,1]);
                      
[Lvector,vector,Rvector]=positive(double(T),Lvector,vector,Rvector');
output_value_new=vector'*M*reshape(kron(Lvector,Rvector),[d1*d3,1]); 


if  output_value_new>output_value(index)
output_value(index)= output_value_new;
output_vector_X(:,index)=Lvector;
output_vector_Y(:,index)=vector;
output_vector_Z(:,index)=Rvector;
end
                                
                                            
                                             
T=T-output_value(index)*reshape(kron(output_vector_Z(:,index),kron(output_vector_Y(:,index),output_vector_X(:,index))),[d1,d2,d3]);
 
end
                        
    
                                 
end


