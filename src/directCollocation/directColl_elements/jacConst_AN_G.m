function [ gz ] = jacConst_AN_G( dg, gz, M, nt, np, ng, nz, T, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    idx=data.FD.index.g;nfd=size(idx,2);
    dgz=[dg.dx dg.du];
    
    if nt
        Tj=kron(T,ones(1,ng)); 
        gz=gz+sparse(1:M*ng,idx(:,1),reshape((dg.dt{1}.*Tj)',M*ng,1),M*ng,nz);
    end
    for i=nt+1:nfd
      if np&&(nt<i)&&(i<nt+np+1)
        gz=gz+sparse(1:M*ng,idx(:,i),reshape(dg.dp{i-nt}',M*ng,1),M*ng,nz);  
      else   
       gz=gz+sparse(1:M*ng,idx(:,i),reshape(dgz{i-nt-np}',M*ng,1),M*ng,nz);
      end
    end


end

