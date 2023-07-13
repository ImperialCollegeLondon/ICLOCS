function [ fz,Jf ] = jacConst_AN_F( df, fz, M, n, m, np, nt, nz, f, X, U, P, t0, tf, T, e, vdat, data )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

  idx=data.FD.index.f;    
  % Compute and store the derivate of fz with respect to tf

   if nt==1 && ~isempty(df.dt) && ~isempty(cell2mat(df.dt))
     fz_tf=f(X,U,P,(tf-t0)*T+t0,vdat)+(tf-t0)*(df.dt{1}).*(kron(T,ones(1,n)));
     fz=sparse(1:M*n,idx(:,1),reshape(fz_tf.',M*n,1),M*n,nz);
     Jaf=[df.dt];
   elseif nt==2 && ~isempty(df.dt) && ~isempty(cell2mat(df.dt))
     fz_t0=-(f(X,U,P,(tf-t0)*T+t0,vdat)+(tf-t0)*(df.dt{1}).*(kron(T,ones(1,n))));
     fz_tf=f(X,U,P,(tf-t0)*T+t0,vdat)+(tf-t0)*(df.dt{1}).*(kron(T,ones(1,n)));
     fz=fz+sparse(1:M*n,idx(:,1),reshape(fz_t0.',M*n,1),M*n,nz)+sparse(1:M*n,idx(:,2),reshape(fz_tf.',M*n,1),M*n,nz);
     Jaf=[df.dt, df.dt];
   else
         idx=data.FD.index.f;
         et0=e*data.FD.vector.f.et0;etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
         ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;

        for i=1:nt
           fp=f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0-et0(i)).*T+t0+et0(i),vdat);
           fm=f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0+et0(i)).*T+t0-et0(i),vdat);
           ft=((tf+etf(i)-t0-et0(i))*fp-(tf-etf(i)-t0+et0(i))*fm)/(2*e);
           fz=fz+sparse(1:M*n,idx(:,i),reshape(ft',M*n,1),M*n,nz);
        end
        Jaf=[cell(nt,nt)];
   end
  if np 
  for i=1:np
      if ~isempty(df.dp{i})
            fz=fz+sparse(1:M*n,idx(:,nt+i),reshape((tf-t0)*df.dp{i},M*n,1),M*n,nz);
      end
  end
   Jaf=[Jaf, df.dp];
  end
  for i=1:n
     
  fz=fz+sparse(1:M*n,idx(:,nt+np+i),reshape((tf-t0)*df.dx{i}',M*n,1),M*n,nz); 
  
  end
  for i=1:m
   fz=fz+sparse(1:M*n,idx(:,nt+np+n+i),reshape((tf-t0)*df.du{i}',M*n,1),M*n,nz); 
  end
  Jaf=[Jaf, df.dx, df.du];
  Jf=vertcat(Jaf);
    
    


end

