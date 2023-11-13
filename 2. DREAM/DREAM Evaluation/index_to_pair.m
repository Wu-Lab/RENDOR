function [i,j]=index_to_pair(ind,n_tf)

j=ceil(ind/n_tf);
k=mod(ind,n_tf);
if k==0
    i=n_tf;
else
   i=k; 
end