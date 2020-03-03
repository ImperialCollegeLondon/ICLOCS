function c=g_eq(x,u,p,t,vdat)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% Warp function
%------------- BEGIN CODE --------------
g_unscaled=vdat.InternalDynamics;
ng_group=nargout(g_unscaled);
if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
end

if ng_group==1
    c=[];
elseif ng_group==2
    if vdat.ng_eq 
        [~,c] = g_unscaled(x,u,p,t,vdat);
    end
else
    [~,c,~] = g_unscaled(x,u,p,t,vdat);
end


