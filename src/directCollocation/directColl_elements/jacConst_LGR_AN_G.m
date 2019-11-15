function [ gz ] = jacConst_LGR_AN_G( dg, gz, M, nt, np, ng, nz, T, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    idx=data.FD.index.g;nfd=size(idx,2);
    dgz=[dg.dx dg.du];
    
      for i=1:(m+n)
        gz=gz+sparse(1:M*ng,idx(:,i),reshape(dgz{i},M*ng,1),M*ng,nz);
      end

      if np
          for i=1:np
            gz=gz+sparse(1:M*ng,idx(:,i+m+n),reshape(dg.dp{i},M*ng,1),M*ng,nz);  
          end
      end

      if nt>=2 && ~isempty(dg.dt)
         for i=1:nt
             if i==1
                 Tj=kron((1-T)/2,ones(1,ng)); 
             elseif i==nt
                 Tj=kron((T+1)/2,ones(1,ng)); 
             end
             
             gz=gz+sparse(1:M*ng,idx(:,i+m+n+np),reshape((dg.dt{i}.*Tj),M*ng,1),M*ng,nz);
         end
      elseif nt>=2 && isempty(dg.dt)
            idx=data.FD.index.g;nfd=size(idx,2);
            et0=e*data.FD.vector.g.et0;etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
            ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;
            for i=nfd-nt+1:nfd
                gp=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
                gm=g(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
                gz=gz+sparse(1:M*ng,idx(:,i),reshape((gp-gm)/(2*e),M*ng,1),M*ng,nz);
            end
      end

end

