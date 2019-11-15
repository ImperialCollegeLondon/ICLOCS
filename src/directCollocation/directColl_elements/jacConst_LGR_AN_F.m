function [ fz,Jf ] = jacConst_LGR_AN_F( df, fz, M, n, m, np, nt, nz, f, X, U, P, t0, tf, T, e, DT_ratio_diff, DTLP, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

  idx=data.FD.index.f;    
  % Compute and store the derivate of fz with respect to tf

  for i=1:n
    fz=fz-sparse(1:M*n,idx(:,i),reshape(repmat(data.t_segment_end,1,n).*df.dx{i},M*n,1),M*n,nz); 
  end
  for i=1:m
   fz=fz-sparse(1:M*n,idx(:,n+i),reshape(repmat(data.t_segment_end,1,n).*df.du{i},M*n,1),M*n,nz); 
  end
  Jaf=[df.dx, df.du];
  
  if np
      for i=1:np
        fz=fz-sparse(1:M*n,idx(:,n+m+i),reshape(repmat(data.t_segment_end,1,n).*df.dp{i},M*n,1),M*n,nz);
      end
  Jaf=[Jaf, df.dp];
  end
  
  if nt>=2 && ~isempty(df.dt) && ~isempty(cell2mat(df.dt))
     for i=1:nt
         if i==1 
            fz_t=f(X,U,P,(tf-t0)/2*T+(tf+t0)/2,vdat)*.5-repmat(data.t_segment_end,1,n).*(df.dt{1}).*(kron((1-T)/2,ones(1,n)));
         elseif i==nt 
            fz_t=-f(X,U,P,(tf-t0)/2*T+(tf+t0)/2,vdat)*.5-repmat(data.t_segment_end,1,n).*(df.dt{1}).*(kron((1+T)/2,ones(1,n)));
         end
         fz=fz+sparse(1:M*n,idx(:,n+m+np+i),reshape(fz_t.',M*n,1),M*n,nz);
     end
     Jaf=[Jaf, df.dt];
  elseif nt>=2 && (isempty(df.dt) || ~isempty(cell2mat(df.dt)))
    idx=data.FD.index.f;nfd=size(idx,2);
    et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;et=data.FD.vector.f.et;
    ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

    if data.options.adaptseg==1 
        for i=nfd-nt+1:nfd
            fp=(repmat(data.t_segment_end+data.t_segment_mat_m*(et{i}'*e),1,n)).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
            fm=(repmat(data.t_segment_end-data.t_segment_mat_m*(et{i}'*e),1,n)).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
            fz=fz+sparse(1:M*n,idx(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
        end
    else
        for i=nfd-nt+1:nfd
            fp=(DTLP+etf(i)/2*DT_ratio_diff-et0(i)/2*DT_ratio_diff).*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i))/2.*T+(tf+etf(i)+t0+et0(i))/2,vdat);
            fm=(DTLP-etf(i)/2*DT_ratio_diff+et0(i)/2*DT_ratio_diff).*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i))/2.*T+(tf-etf(i)+t0-et0(i))/2,vdat);
            ft=sparse(1:M*n,idx(:,i),reshape(-(fp-fm)/(2*e),M*n,1),M*n,nz);
            fz=fz+ft;
        end
    end
     Jaf=[Jaf, cell(1,1)];
  end
  

  

  Jf=vertcat(Jaf);    
    


end

