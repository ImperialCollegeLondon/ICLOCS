function [ Lz, JL ] = gradientCost_LGR_AN_L( dL, Lz,L,nt,np,n,m,nz,X,Xr,U,Ur,P,t,T,t0,tf,vdat,data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    idx=data.FD.index.Ly;
    for  i=1:n
       dLx=data.map.w.*(data.t_segment_end.*dL.dx(:,i));
       Lz=Lz+sparse(1,idx(:,i),dLx,1,nz);
       JaL{i}=dL.dx(:,i);
    end  
    

    for i=1:m
       dLu=data.map.w.*(data.t_segment_end.*dL.du(:,i));  
       Lz=Lz+sparse(1,idx(:,n+i),dLu,1,nz);
       JaL{i+n}=dL.du(:,i);
    end
    

    if np
        for i=1:np
            dLp=data.map.w*(data.t_segment_end.*dL.dp(:,i));
            Lz=Lz+sparse(1,idx(:,n+m+i),dLp,1,nz);
            JaL{i+n+m}=dL.dp(:,i);
        end
    end
    

    if nt==2 && ~isempty(dL.dt)
        for i=1:nt
            if i==1 
                dLt=-0.5*data.map.w'*L(X,Xr,U,Ur,P,t,vdat)+data.map.w'*(data.t_segment_end.*(dL.dt.*(1-T)/2)); %t0
            elseif i==2 
                dLt=0.5*data.map.w'*L(X,Xr,U,Ur,P,t,vdat)+data.map.w'*(data.t_segment_end.*(dL.dt.*(1+T)/2)); %tf
            end
            Lz=Lz+sparse(1,idx(:,n+m+np+i),dLt,1,nz);
            JaL{i+n+m+np}=dL.dt(:);
        end
    end
   
    
    JL=vertcat(JaL);
    
end

