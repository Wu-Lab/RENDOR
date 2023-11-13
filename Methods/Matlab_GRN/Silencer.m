function [output_network]=Silencer(mat)

[n_tf,n]=size(mat);
for i=1:n_tf
    mat(i,i)=0;
end

%% ********************** input matrix imputation ***************** 
mat(1:n_tf,1:n_tf)=(mat(1:n_tf,1:n_tf)+mat(1:n_tf,1:n_tf)')/2;
mat1=[mat;[zeros(n-n_tf,n_tf),eye(n-n_tf,n-n_tf)]]; 
mat1=(mat1+mat1')/2;
mat1=(mat1-min(mat1(:)))/(max(mat1(:))-min(mat1(:)));
mat1=mat1+min(mat1(mat1>0));
% mat1=corrcoef(mat1);
% mat1= mat1+eye(n)+min(mat1(:));

% mat1=[mat;[mat(:,(n_tf+1):end)',eye(n-n_tf,n-n_tf)]]; 
% mat1=(mat1+mat1')/2;
% mat1=(mat1-min(mat1(:)))./(max(mat1(:))-min(mat1(:)));
% mat1=mat1+min(mat1(mat1>0));

%% ********************** Silencer *********************************
temp = diag(diag((mat1-eye(n))*mat1));
net_new = (mat1-eye(n)+temp)*inv(mat1);
net_new = abs(net_new);
% net_new=(net_new-min(net_new(:)))./(max(net_new(:))-min(net_new(:)));

%% ****************************************************************
output_network=net_new(1:n_tf,:);

