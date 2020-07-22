function c=g(x,u,p,t,vdat)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% Warp function
%------------- BEGIN CODE --------------
g_unscaled=vdat.InternalDynamics;
ng_group=nargout(g_unscaled);
if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale_back, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale_back, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale_back, vdat.Pshift );
    end
end

if ng_group==1
    c=[];
elseif ng_group==2
    if vdat.resmin && vdat.ng_eq 
        c=[];
    else
        [~,c] = g_unscaled(x,u,p,t,vdat);
    end
else
    [~,ceq,cneq] = g_unscaled(x,u,p,t,vdat);
    if vdat.resmin && vdat.ng_eq 
        c=cneq;
    else
        c=[ceq cneq];
    end
end

if isfield(vdat,'gFilter')
    c(:,vdat.gFilter)=[];
end

if strcmp(vdat.mode.currentMode,'Feasibility')
    c=[c-p(:,end-vdat.mode.np*2+1:end-vdat.mode.np) c+p(:,end-vdat.mode.np+1:end)];
end
