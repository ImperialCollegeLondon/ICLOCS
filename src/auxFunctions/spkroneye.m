function [ A_new ] = spkroneye( N,A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[idx_i,idx_j,mat_val] = find(A);

N_row=size(A,1);
N_col=size(A,2);
M=length(idx_i);

idx_i_new=zeros(M*N,1);
idx_j_new=zeros(M*N,1);
mat_val_new=zeros(M*N,1);

for i=1:N
    idx_i_new((i-1)*M+1:i*M)=idx_i+(i-1)*N_row;
    idx_j_new((i-1)*M+1:i*M)=idx_j+(i-1)*N_col;
    mat_val_new((i-1)*M+1:i*M)=mat_val;
end

A_new=sparse(idx_i_new,idx_j_new,mat_val_new,N_row*N,N_col*N);

end
