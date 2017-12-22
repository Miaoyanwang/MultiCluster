% Semi-nonnegative constraint
function[Lvector,Rvector,vector]=positive(T,vector1,vector2,vector3)

d=size(T);
d1=d(1);
d2=d(2);
d3=d(3);

%if sumsqr(vector3(vector3<0))>=sumsqr(vector3(vector3>0));
if -min(vector3)>max(vector3)
vector3=-vector3;
end

vector3(vector3<0)=0;
M=reshape(double(T),[d1*d2,d3]);
M_hat=reshape(M*vector3',[d1,d2]);
[Lvector,~,Rvector]=svds(M_hat,1);
vector=vector3/norm(vector3);   
              
end


