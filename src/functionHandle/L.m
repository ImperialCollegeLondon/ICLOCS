function stageCost=L(x,xr,u,ur,p,t,vdat)

% L - Returns the stage cost.
% Warp function
%------------- BEGIN CODE --------------
L_unscaled=vdat.functions_unscaled{1};
if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale_back, vdat.Xshift );
    if ~isempty(xr)
        xr=scale_variables_back( xr, vdat.Xscale_back, vdat.Xshift );
    end
    u=scale_variables_back( u, vdat.Uscale_back, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale_back, vdat.Pshift );
    end
    if ~isempty(ur)
        ur=scale_variables_back( ur, vdat.Uscale_back, vdat.Ushift );
    end
    if strcmp(vdat.mode.currentMode,'Feasibility')
        stageCost=0*t;
    else
        stageCost=L_unscaled(x,xr,u,ur,p,t,vdat);
    end
else
    if strcmp(vdat.mode.currentMode,'Feasibility')
        stageCost=0*t;
    else
        stageCost=L_unscaled(x,xr,u,ur,p,t,vdat);
    end
end

%------------- END OF CODE --------------