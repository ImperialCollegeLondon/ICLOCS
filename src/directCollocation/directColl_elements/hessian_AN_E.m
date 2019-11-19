function [ Ezz ] = hessian_AN_E( Ezz, HE, nt, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if iscell(HE)
      idx=data.FD.index.Ey; 
      idxv=data.FD.vindex.E;
      nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
      if nfd&& (nt==1 || length(idx)==1)
       Ezz(idx(1),idx(1))=HE{1,1};
     elseif nfd&&nt==2
       Ezz(idx(1),idx(1))=-HE{1,1};
       Ezz(idx(2),idx(2))=HE{2,2};
       Ezz(idx(2),idx(1))=-HE{1,2};
      end    
    for i=1+nt:nfd
        for j=1:i
            if (nt==1&&(j==1))
                Ezz(idx(i),idx(j))=HE{idxv(j),idxv(i)};
            elseif (nt==2&&(j<=nt))
                if j==1
                    Ezz(idx(i),idx(j))=-HE{idxv(j),idxv(i)};
                elseif j==2
                    Ezz(idx(i),idx(j))=HE{idxv(j),idxv(i)};
                end
            else   
                Ezz(idx(i),idx(j))=HE{idxv(j),idxv(i)};
            end
        end
    end
else
      idx=data.FD.index.Ey; 
      idxv=data.FD.vindex.E;
      nfd=size(idx,2);  % It works when the structure is determined for the worst case (not with random variables)
      if nfd&& (nt==1 || length(idx)==1)
       Ezz(idx(1),idx(1))=HE(1,1);
     elseif nfd&&nt==2
       Ezz(idx(1),idx(1))=-HE(1,1);
       Ezz(idx(2),idx(2))=HE(2,2);
       Ezz(idx(2),idx(1))=-HE(1,2);
      end    
    for i=1+nt:nfd
        for j=1:i
            if (nt==1&&(j==1))
                Ezz(idx(i),idx(j))=HE(idxv(j),idxv(i));
            elseif (nt==2&&(j<=nt))
                if j==1
                    Ezz(idx(i),idx(j))=-HE(idxv(j),idxv(i));
                elseif j==2
                    Ezz(idx(i),idx(j))=HE(idxv(j),idxv(i));
                end
            else   
                Ezz(idx(i),idx(j))=HE(idxv(j),idxv(i));
            end
        end
    end
end


end

