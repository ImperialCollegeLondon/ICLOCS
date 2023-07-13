function [ Lzz ] = hessian_AN_wL( dL, Lzz, HL, M, nt, nz, T, DT, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 idx=data.FD.index.Ly; 
 nfd=size(idx,2);  
 if nfd&&nt==1
   Lt=(2*dL{1}.*T+DT*HL{1,1}.*(T.^2)).*data.map.w;  
   Lzz=Lzz+sparse(idx(:,1),idx(:,1),reshape(Lt',M,1),nz,nz);  
 elseif nfd&&nt==2
   Lt0t0=-(2*dL{1}.*T+DT*HL{1,1}.*(T.^2)).*data.map.w;  
   Ltftf=(2*dL{2}.*T+DT*HL{2,2}.*(T.^2)).*data.map.w;  
   Lt0tf=-(2*dL{1}.*T+DT*HL{1,2}.*(T.^2)).*data.map.w;  
   Lzz=Lzz+sparse(idx(:,1),idx(:,1),reshape(Lt0t0',M,1),nz,nz)+sparse(idx(:,2),idx(:,2),reshape(Ltftf',M,1),nz,nz)+sparse(idx(:,2),idx(:,1),reshape(Lt0tf',M,1),nz,nz);  
 end    
 
  for i=1+nt:nfd
   for j=1:i
       if any(HL{j,i},'all')
          if (nt==1&&(j==1))
            Lt=(dL{i}+DT*HL{j,i}.*T).*data.map.w;
          elseif (nt==2&&(j<=nt))
            if j==1
                Lt=-(dL{i}+DT*HL{j,i}.*T).*data.map.w;
            elseif j==2
                Lt=(dL{i}+DT*HL{j,i}.*T).*data.map.w;
            end
          else   
            Lt=DT*HL{j,i}.*data.map.w;
          end
         Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);
       end
    end
  end
    
    


end

